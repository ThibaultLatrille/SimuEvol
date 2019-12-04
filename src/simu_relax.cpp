#include "argparse.hpp"
#include "random.hpp"
#include "statistic.hpp"

using namespace TCLAP;
using namespace std;

normal_distribution<double> normal_distr(0.0, 1.0);
uniform_int_distribution<u_long> chr_distr(0, 1);

double distance(vector<double> const &v) {
    return sqrt(
        accumulate(v.begin(), v.end(), 0.0, [](double a, double const &b) { return a + b * b; }));
}

vector<double> spherical_coord(u_long n, double radius) {
    vector<double> coord(n, 0);

    for (u_long i = 0; i < n; ++i) { coord[i] = normal_distr(generator); }

    double norm = radius / distance(coord);

    for (u_long i = 0; i < n; ++i) { coord[i] *= norm; }
    return coord;
}

class Chromosome {
  public:
    vector<double> phenotype;

    explicit Chromosome() = default;

    explicit Chromosome(u_long n, double radius) { phenotype = spherical_coord(n, radius); };

    void add_mutation(vector<double> &mut) {
        assert(mut.size() == phenotype.size());
        for (u_long i = 0; i < mut.size(); ++i) { phenotype[i] += mut[i]; }
    };
};

class Being {
  public:
    vector<Chromosome> genome;
    double fitness{1.0};
    double d{0.0};

    static u_long k;
    static double mu;
    static double r;
    static u_long complexity;
    static u_long pleiotropy;
    static double a;
    static double q;

    explicit Being() {
        genome.reserve(2 * k);
        auto chr = Chromosome(Being::complexity, 1.0 / (2 * k));
        for (u_long i = 0; i < k; ++i) {
            genome.push_back(chr);
            genome.push_back(chr);
        }
        update_fitness();
    };

    Being(Being const &mum, Being const &dad) {
        genome.reserve(2 * k);
        for (u_long i = 0; i < k; ++i) {
            genome.push_back(mum.genome[2 * i + chr_distr(generator)]);
            genome.push_back(dad.genome[2 * i + chr_distr(generator)]);
        }
        update_fitness();
    };

    void update_fitness() {
        vector<double> p(complexity, 0.0);
        for (auto const &chr : genome) {
            for (u_long i = 0; i < complexity; ++i) { p[i] += chr.phenotype[i]; }
        }
        d = distance(p);
        fitness = exp(-a * pow(d, q));
    };

    void mutation() {
        poisson_distribution<u_long> poisson_distr(mu);
        u_long poisson_draw = poisson_distr(generator);

        if (poisson_draw > 0) {
            uniform_int_distribution<u_long> chr_distr(0, 2 * k - 1);

            for (u_long mutation{0}; mutation < poisson_draw; mutation++) {
                auto mutation_coord = spherical_coord(pleiotropy, r);
                if (pleiotropy < complexity) {
                    mutation_coord.resize(complexity, 0.0);
                    shuffle(mutation_coord.begin(), mutation_coord.end(), generator);
                }
                genome[chr_distr(generator)].add_mutation(mutation_coord);
            };

            update_fitness();
        }
    };
};

class Population {
  public:
    // Parameters
    u_long population_size{100};
    u_long elapsed{0};

    // Beings composing the population
    vector<Being> parents{};
    double delta_f{0};
    double delta_f_5pct_low{0};
    double delta_f_5pct_high{0};
    double delta_f_50pct_low{0};
    double delta_f_50pct_high{0};

    explicit Population(u_long population_size) : population_size{population_size} {
        parents.resize(population_size);
    };

    void mutation() {
        for (auto &being : parents) { being.mutation(); }
    };

    void selection() {
        vector<double> fitness_beings(parents.size(), 0.0);
        transform(parents.begin(), parents.end(), fitness_beings.begin(),
            [](Being &b) { return b.fitness; });
        discrete_distribution<u_long> parent_distr(fitness_beings.begin(), fitness_beings.end());

        vector<Being> children;
        children.reserve(parents.size());
        for (u_long i_being{0}; i_being < parents.size(); i_being++) {
            u_long mum = parent_distr(generator);
            u_long dad = mum;
            while (dad == mum) { dad = parent_distr(generator); }
            children.emplace_back(parents[mum], parents[dad]);
        }
        mean_delta_f(parents, children);
        parents = move(children);
    };

    void random_mating() {
        uniform_int_distribution<u_long> parent_distr(0, parents.size() - 1);

        vector<Being> children;
        children.reserve(parents.size());
        for (u_long i_being{0}; i_being < parents.size(); i_being++) {
            u_long mum = parent_distr(generator);
            u_long dad = mum;
            while (dad == mum) { dad = parent_distr(generator); }
            children.emplace_back(parents[mum], parents[dad]);
        }
        mean_delta_f(parents, children);
        parents = move(children);
    };

    void run(
        string &output_filename, u_long nbr_intervals, u_long interval_length, u_long relax_time) {
        // Run under selection
        cout << "Running under selection" << endl;
        output_state(output_filename);
        u_long last_pct = 0;
        for (u_long sample{1}; sample <= nbr_intervals; sample++) {
            for (u_long gen{1}; gen <= interval_length; gen++) {
                mutation();
                selection();
                elapsed++;
            }
            output_state(output_filename);
            u_long pct = 100 * sample / nbr_intervals;
            if (pct != last_pct) {
                cout << pct << "% (" << interval_length * sample << " generations)" << endl;
                last_pct = pct;
            }
        }
        cout << "Running under relaxation" << endl;
        // Relax selection
        for (u_long relax{1}; relax <= relax_time; relax++) {
            mutation();
            random_mating();
            elapsed++;
            output_state(output_filename);
        }
    };

    void mean_delta_f(vector<Being> const &pop, vector<Being> const &next_pop) {
        vector<double> delta_fitnesses(population_size, 0);
        for (size_t i = 0; i < delta_fitnesses.size(); ++i) {
            delta_fitnesses[i] = next_pop[i].fitness - pop[i].fitness;
        }
        sort(delta_fitnesses.begin(), delta_fitnesses.end());
        delta_f = avg(delta_fitnesses);
        auto low_bound = static_cast<u_long>(0.05 * (population_size - 1));
        auto fifty_bound = static_cast<u_long>(0.5 * (population_size - 1));
        auto up_bound = static_cast<u_long>(0.95 * (population_size - 1));
        delta_f_5pct_low =
            accumulate(delta_fitnesses.begin(), delta_fitnesses.begin() + low_bound, 0.0) /
            low_bound;
        delta_f_5pct_high =
            accumulate(delta_fitnesses.begin() + up_bound, delta_fitnesses.end(), 0.0) /
            (population_size - up_bound);

        delta_f_50pct_low =
            accumulate(delta_fitnesses.begin(), delta_fitnesses.begin() + fifty_bound, 0.0) /
            fifty_bound;
        delta_f_50pct_high =
            accumulate(delta_fitnesses.begin() + fifty_bound, delta_fitnesses.end(), 0.0) /
            (population_size - fifty_bound);
    }

    void output_state(string &output_filename) const {
        auto low_bound = static_cast<u_long>(0.05 * (population_size - 1));
        auto up_bound = static_cast<u_long>(0.95 * (population_size - 1));

        vector<double> fitnesses(population_size, 0);
        transform(parents.begin(), parents.end(), fitnesses.begin(),
            [](Being const &b) { return b.fitness; });

        vector<double> distances(population_size, 0);
        transform(
            parents.begin(), parents.end(), distances.begin(), [](Being const &b) { return b.d; });

        sort(fitnesses.begin(), fitnesses.end());
        sort(distances.begin(), distances.end());

        double f_tot = sum(fitnesses);
        double entropy = 0;
        for (auto const &f : fitnesses) { entropy += f * log(f / f_tot); }
        entropy = exp(-entropy / f_tot) / fitnesses.size();

        vector<double> p(Being::complexity, 0.0);
        for (auto const &being : parents) {
            for (auto const &chr : being.genome) {
                for (u_long i = 0; i < Being::complexity; ++i) { p[i] += chr.phenotype[i]; }
            }
        }
        double d = distance(p);
        double fitness = exp(-Being::a * pow(d, Being::q));

        ofstream tsv;
        tsv.open(output_filename + ".tsv", ios_base::app);
        tsv << elapsed << "\t" << 1 << "\tFitness lower bound\t" << fitnesses[low_bound] << endl;
        tsv << elapsed << "\t" << 1 << "\tFitness upper bound\t" << fitnesses[up_bound] << endl;
        tsv << elapsed << "\t" << 1 << "\tMean fitness\t" << avg(fitnesses) << endl;
        tsv << elapsed << "\t" << 1 << "\tFitness of mean phenotype\t" << fitness << endl;
        tsv << elapsed << "\t" << 1 << "\tFitness entropy\t" << entropy << endl;
        tsv << elapsed << "\t" << 2 << "\tDistance lower bound\t" << distances[low_bound] << endl;
        tsv << elapsed << "\t" << 2 << "\tDistance upper bound\t" << distances[up_bound] << endl;
        tsv << elapsed << "\t" << 2 << "\tMean distance\t" << avg(distances) << endl;
        tsv << elapsed << "\t" << 3 << "\tDistance of mean phenotype\t" << d << endl;
        tsv << elapsed << "\t" << 4 << "\tΔF mean\t" << delta_f << endl;
        tsv << elapsed << "\t" << 4 << "\tΔF 5% low\t" << delta_f_5pct_low << endl;
        tsv << elapsed << "\t" << 4 << "\tΔF 5% high\t" << delta_f_5pct_high << endl;
        tsv << elapsed << "\t" << 4 << "\tΔF 50% low\t" << delta_f_50pct_low << endl;
        tsv << elapsed << "\t" << 4 << "\tΔF 50% high\t" << delta_f_50pct_high << endl;
        tsv << elapsed << "\t" << 0 << "\tVar fitness\t" << variance(fitnesses) << endl;
        tsv << elapsed << "\t" << 0 << "\tVar distance\t" << variance(distances) << endl;
        tsv.close();
    };
};

u_long Being::k = 1;
double Being::mu = 0;
double Being::r = 0;
u_long Being::complexity = 1;
u_long Being::pleiotropy = 1;
double Being::a = 0;
double Being::q = 0;

class SimuEvolArgParse : public OutputArgParse {
  public:
    explicit SimuEvolArgParse(CmdLine &cmd) : OutputArgParse(cmd) {}

    TCLAP::ValueArg<u_long> pop_size{"n", "pop_size", "population size", false, 100, "u_long", cmd};
    TCLAP::ValueArg<u_long> chromosomes{
        "k", "chromosomes", "number of chromosomes per individual", false, 23, "u_long", cmd};
    TCLAP::ValueArg<double> mu{"m", "mutation_rate_per_generation",
        "mean number of mutations per individual per generation", false, 10, "double", cmd};
    TCLAP::ValueArg<double> radius{"r", "radius",
        "effect of a mutation (radius of the move in the phenotypic space)", false, 1e-2, "double",
        cmd};
    TCLAP::ValueArg<u_long> complexity{"c", "complexity",
        "complexity (dimension) of the phenotypic space", false, 10, "u_long", cmd};
    TCLAP::ValueArg<u_long> pleiotropy{"p", "pleiotropy",
        "pleiotropy of mutations (p<=c, and p=c if not specified)", false, 0, "u_long", cmd};
    TCLAP::ValueArg<double> peakness{"a", "peakness",
        "'a' parameter (peakness) of the fitness function (exp(-a*(d^q))", false, 0.5, "double",
        cmd};
    TCLAP::ValueArg<double> epistasis{"q", "epistasis",
        "'q' parameter (epistasis) of fitness function (exp(-a*(d^q))", false, 3.0, "double", cmd};
};

int main(int argc, char *argv[]) {
    CmdLine cmd{"SimuRelax", ' ', "0.1"};
    SimuEvolArgParse args(cmd);
    cmd.parse(argc, argv);

    u_long arg_seed = args.seed.getValue();
    cout << "Random generator seed: " << arg_seed << endl;
    generator.seed(arg_seed);

    string output_path{args.output_path.getValue()};

    u_long pop_size{args.pop_size.getValue()};
    cout << "Population of " << pop_size << " individuals." << endl;

    Being::k = args.chromosomes.getValue();
    cout << "Each individual has " << Being::k << " pairs of chromosomes." << endl;

    Being::mu = args.mu.getValue();
    cout << "Each individual has on average " << Being::mu << " mutations per generations." << endl;

    Being::r = args.radius.getValue();
    cout << "Each mutation move the phenotype at a distance " << Being::r << " from the parent."
         << endl;

    Being::complexity = args.complexity.getValue();
    cout << "The phenotypic space is of dimension " << Being::complexity << "." << endl;

    Being::pleiotropy = args.pleiotropy.getValue();
    if (Being::pleiotropy == 0) {
        Being::pleiotropy = Being::complexity;
    } else {
        assert(Being::pleiotropy <= Being::complexity);
    }
    cout << "Each mutation has a pleiotropic effect on " << Being::pleiotropy << " dimensions."
         << endl;

    Being::a = args.peakness.getValue();
    cout << "The peakness of the fitness function is " << Being::a << "." << endl;

    Being::q = args.epistasis.getValue();
    cout << "The epistasis of the fitness function is " << Being::q << "." << endl;

    u_long nbr_intervals = 500;
    u_long interval_length = 1;
    u_long relax_length = 500;

    output_path += "/" + to_string(pop_size) + "_" + to_string(Being::k) + "_" +
                   to_string(Being::mu) + "_" + to_string(Being::r) + "_" +
                   to_string(Being::complexity) + "_" + to_string(Being::pleiotropy) + "_" +
                   to_string(Being::a) + "_" + to_string(Being::q);
    cout << "The data will be written in " << output_path << endl;


    ofstream tsv;
    tsv.open(output_path + ".tsv");
    tsv << "seed\t" << arg_seed << endl;
    tsv << "pop_size\t" << pop_size << endl;
    tsv << "chromosomes\t" << Being::k << endl;
    tsv << "mutation_rate_per_generation\t" << Being::mu << endl;
    tsv << "radius\t" << Being::r << endl;
    tsv << "complexity\t" << Being::complexity << endl;
    tsv << "pleiotropy\t" << Being::pleiotropy << endl;
    tsv << "peakness\t" << Being::a << endl;
    tsv << "epistasis\t" << Being::q << endl;
    tsv << "t\t" << nbr_intervals * interval_length << endl;
    tsv.close();

    assert(pop_size > 1);
    assert(Being::complexity > 0);
    assert(Being::pleiotropy > 0);
    assert(Being::k > 0);
    Population population(pop_size);
    population.run(output_path, nbr_intervals, interval_length, relax_length);

    return 0;
}