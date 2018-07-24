#include <iostream>
#include <random>
#include <fstream>
#include <set>
#include <iomanip>

using namespace std;

#define DOCOPT_HEADER_ONLY

#include "docopt.cpp/docopt.h"

double seed{0};
default_random_engine generator(seed);
normal_distribution<double> normal_distr(0.0, 1.0);
uniform_int_distribution<u_long> chr_distr(0, 1);

double distance(vector<double> const &v) {
    return sqrt(accumulate(v.begin(), v.end(), 0.0, [](double a, double const &b) { return a + b * b; }));
}

vector<double> spherical_coord(u_long n, double radius) {
    vector<double> coord(n, 0);

    for (int i = 0; i < n; ++i) {
        coord[i] = normal_distr(generator);
    }

    double norm = radius / distance(coord);

    for (int i = 0; i < n; ++i) {
        coord[i] *= norm;
    }
    return coord;
}

// Mean of a vector
double mean(vector<double> const &v) {
    return accumulate(v.begin(), v.end(), 0.0) / v.size();
}

// Variance of a vector
double var(vector<double> const &v) {
    double s2 = accumulate(v.begin(), v.end(), 0.0, [](double a, double const &b) { return a + b * b; }) / v.size();
    double s = mean(v);
    return s2 - s * s;
}

class Chromosome {
public:
    vector<double> phenotype;

    explicit Chromosome() = default;

    explicit Chromosome(u_long n, double radius) {
        phenotype = spherical_coord(n, radius);
    };

    void add_mutation(vector<double> &mut) {
        assert(mut.size() == phenotype.size());
        for (int i = 0; i < mut.size(); ++i) {
            phenotype[i] += mut[i];
        }
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
    static u_long n;
    static u_long m;
    static double a;
    static double q;

    explicit Being() {
        genome.reserve(2 * k);
        auto chr = Chromosome(Being::n, 1.0 / (2 * k));
        for (int i = 0; i < k; ++i) {
            genome.push_back(chr);
            genome.push_back(chr);
        }
        update_fitness();
    };

    Being(Being const &mum, Being const &dad) {
        genome.reserve(2 * k);
        for (int i = 0; i < k; ++i) {
            genome.emplace_back(mum.genome[2 * i + chr_distr(generator)]);
            genome.emplace_back(dad.genome[2 * i + chr_distr(generator)]);
        }
        update_fitness();
    };

    void update_fitness() {
        vector<double> p(n, 0.0);
        for (auto const &chr: genome) {
            for (int i = 0; i < n; ++i) {
                p[i] += chr.phenotype[i];
            }
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
                auto mutation_coord = spherical_coord(m, r);
                if (m < n) {
                    mutation_coord.resize(n, 0.0);
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
    vector<Being> beings{};

    explicit Population(u_long population_size) : population_size{population_size} {
        beings.resize(population_size);
    };

    void mutation() {
        for (auto &being: beings) {
            being.mutation();
        }
    };

    void selection() {
        vector<double> fitness_beings(beings.size(), 0.0);
        transform(beings.begin(), beings.end(), fitness_beings.begin(), [](Being &b) {
            return b.fitness;
        });
        discrete_distribution<u_long> parent_distr(fitness_beings.begin(), fitness_beings.end());

        vector<Being> next_beings;
        next_beings.reserve(beings.size());
        for (u_long i_being{0}; i_being < beings.size(); i_being++) {
            u_long mum = parent_distr(generator);
            u_long dad = mum;
            while (dad == mum) {
                dad = parent_distr(generator);
            }
            next_beings.emplace_back(Being(beings[mum], beings[dad]));
        }
        beings = move(next_beings);
    };

    void random_mating() {
        uniform_int_distribution<u_long> parent_distr(0, beings.size() - 1);

        vector<Being> next_beings;
        next_beings.reserve(beings.size());
        for (u_long i_being{0}; i_being < beings.size(); i_being++) {
            u_long mum = parent_distr(generator);
            u_long dad = mum;
            while (dad == mum) {
                dad = parent_distr(generator);
            }
            next_beings.emplace_back(Being(beings[mum], beings[dad]));
        }
        beings = move(next_beings);
    };

    void run(string &output_filename, u_long nbr_intervals, u_long interval_length, u_long relax_time) {
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

    void output_state(string &output_filename) const {
        auto low_bound = static_cast<u_long>(0.05 * (population_size - 1));
        auto up_bound = static_cast<u_long>(0.95 * (population_size - 1));

        vector<double> fitnesses(population_size, 0);
        vector<double> distances(population_size, 0);

        transform(beings.begin(), beings.end(), fitnesses.begin(), [](Being const &b) {
            return b.fitness;
        });
        transform(beings.begin(), beings.end(), distances.begin(), [](Being const &b) {
            return b.d;
        });

        sort(fitnesses.begin(), fitnesses.end());
        sort(distances.begin(), distances.end());

        double f_tot = accumulate(fitnesses.begin(), fitnesses.end(), 0.0);
        double entropy = 0;
        for (auto const &f: fitnesses) {
            entropy += f * log(f / f_tot);
        }
        entropy = exp(-entropy / f_tot) / fitnesses.size();

        auto min_over_max = max_element(fitnesses.begin(), fitnesses.end());
        
        vector<double> p(Being::n, 0.0);
        for (auto const &being: beings) {
            for (auto const &chr: being.genome) {
                for (int i = 0; i < Being::n; ++i) {
                    p[i] += chr.phenotype[i];
                }
            }
        }
        double d = distance(p);
        double fitness = exp(-Being::a * pow(d, Being::q));

        ofstream tsv;
        tsv.open(output_filename + ".tsv", ios_base::app);
        tsv << elapsed << "\t" << 1 << "\tFitness lower bound\t" << fitnesses[low_bound] << endl;
        tsv << elapsed << "\t" << 1 << "\tFitness upper bound\t" << fitnesses[up_bound] << endl;
        tsv << elapsed << "\t" << 1 << "\tMean fitness\t" << mean(fitnesses) << endl;
        tsv << elapsed << "\t" << 1 << "\tFitness of mean phenotype\t" << fitness << endl;
        tsv << elapsed << "\t" << 1 << "\tFitness entropy\t" << entropy << endl;

        tsv << elapsed << "\t" << 2 << "\tDistance lower bound\t" << distances[low_bound] << endl;
        tsv << elapsed << "\t" << 2 << "\tDistance upper bound\t" << distances[up_bound] << endl;
        tsv << elapsed << "\t" << 2 << "\tMean distance\t" << mean(distances) << endl;

        tsv << elapsed << "\t" << 3 << "\tDistance of mean phenotype\t" << d << endl;

        tsv << elapsed << "\t" << 0 << "\tVar fitness\t" << var(fitnesses) << endl;
        tsv << elapsed << "\t" << 0 << "\tVar distance\t" << var(distances) << endl;
        tsv.close();
    };
};

u_long Being::k = 1;
double Being::mu = 0;
double Being::r = 0;
u_long Being::n = 1;
u_long Being::m = 1;
double Being::a = 0;
double Being::q = 0;

static char const USAGE[] =
        R"(
Usage:
      SimuRelax [--pop_size=<100>] [--k=<23>] [--mu=<1e-3>] [--r=<1e-3>] [--n=<10>] [--m=<10>] [--a=<0.5>] [--q=<2.0>] [--dir=<path>]
      SimuRelax --help
      SimuRelax --version

Options:
-h --help                    show this help message and exit
--version                    show version and exit
--pop_size=<100>             specify the population size [default: 100]
--k=<23>                     specify the number of chromosomes per individual [default: 23]
--mu=<10.0>                  specify the mean number of mutations per individual per generation [default: 10.0]
--r=<0.01>                   specify the effect of a mutation (radius of the move in the phenotypic space) [default: 1e-2]
--n=<10>                     specify the complexity (dimension) of the phenotypic space [default: 10]
--m=<10>                     specify the pleiotropy of mutations (m<=n, and m=n if not specified) [default: 10]
--a=<1.0>                    specify the a parameter (peakness) of the fitness function (exp(-a*(d^q)) [default: 0.5]
--q=<2.0>                    specify the q parameter (epistasis) of fitness function (exp(-a*(d^q)) [default: 2.0]
--dir=<path>                 specify the output data folder [default: ../data_relax]

)";

int main(int argc, char *argv[]) {

    auto args = docopt::docopt(USAGE,
                               {argv + 1, argv + argc},
                               true,              // show help if requested
                               "SimuRelax 0.1");  // version string

    u_long pop_size = 100;
    if (args["--pop_size"]) {
        pop_size = static_cast<u_long>(args["--pop_size"].asLong());
    }
    cout << "Population of " << pop_size << " individuals." << endl;

    Being::k = 23;
    if (args["--k"]) {
        Being::k = static_cast<u_long>(args["--k"].asLong());
    }
    cout << "Each individual has " << Being::k << " pairs of chromosomes." << endl;

    Being::mu = 10.0;
    if (args["--mu"]) {
        Being::mu = stod(args["--mu"].asString());
    }
    cout << "Each individual has on average " << Being::mu << " mutations per generations." << endl;

    Being::r = 1e-2;
    if (args["--r"]) {
        Being::r = stod(args["--r"].asString());
    }
    cout << "Each mutation move the phenotype at a distance " << Being::r << " from the parent." << endl;

    Being::n = 10;
    if (args["--n"]) {
        Being::n = static_cast<u_long>(args["--n"].asLong());
    }
    cout << "The phenotypic space is of dimension " << Being::n << "." << endl;

    Being::m = Being::n;
    if (args["--m"]) {
        Being::m = static_cast<u_long>(args["--m"].asLong());
        assert(Being::m <= Being::n);
    }
    cout << "Each mutation has a pleiotropic effect on " << Being::m << " dimensions." << endl;

    Being::a = 0.5;
    if (args["--a"]) {
        Being::a = stod(args["--a"].asString());
    }
    cout << "The flatness of the fitness function is " << Being::a << "." << endl;

    Being::q = 3.0;
    if (args["--q"]) {
        Being::q = stod(args["--q"].asString());
    }
    cout << "The epistasis of the fitness function is " << Being::q << "." << endl;

    u_long nbr_intervals = 1000;
    u_long interval_length = 10;
    u_long relax_length = 1000;

    string path{"../data_relax"};
    if (args["--dir"]) {
        path = args["--dir"].asString();
    }
    path += "/" + to_string(pop_size) + "_" + to_string(Being::k) + "_" + to_string(Being::mu) + "_" +
            to_string(Being::r) + "_" + to_string(Being::n) + "_" + to_string(Being::m) + "_" + to_string(Being::a) +
            "_" + to_string(Being::q);
    cout << "The data will be written in " << path << endl;


    ofstream tsv;
    tsv.open(path + ".tsv");
    tsv << "ne\t" << pop_size << endl;
    tsv << "k\t" << Being::k << endl;
    tsv << "mu\t" << Being::mu << endl;
    tsv << "r\t" << Being::r << endl;
    tsv << "n\t" << Being::n << endl;
    tsv << "m\t" << Being::m << endl;
    tsv << "a\t" << Being::a << endl;
    tsv << "q\t" << Being::q << endl;
    tsv << "t\t" << nbr_intervals * interval_length << endl;
    tsv.close();

    assert(pop_size > 1);
    assert(Being::n > 0);
    assert(Being::m > 0);
    assert(Being::k > 0);
    Population population(pop_size);
    population.run(path, nbr_intervals, interval_length, relax_length);

    return 0;
}