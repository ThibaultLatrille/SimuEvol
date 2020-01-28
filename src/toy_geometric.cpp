#include "argparse.hpp"
#include "random.hpp"
#include "statistic.hpp"

using namespace TCLAP;
using namespace std;

class Sequence {
    unsigned nbr_sites;
    unsigned nbr_right;
    double pop_size;
    double distance;
    double radius;
    double peakness;
    double epistasis;
    std::unordered_map<std::string, SummaryStatistic> summary_stats;

    double computeFitness(double const &d) { return exp(-peakness * pow(abs(d), epistasis)); }

    double selection_coef(double const &d, double const &dMutant) {
        double fm = computeFitness(dMutant), f = computeFitness(d);
        return (fm - f) / f;
    }

    double Pfix(double const &beta, double const &d, double const &delta_d) {
        double S = 4 * selection_coef(d, d + delta_d) * beta;
        if ((abs(S)) < 1e-4) {
            return 1 + S / 2;
        } else {
            return S / (1.0 - exp(-S));
        }
    }

    void next_substitution(bool burn_in) {
        double right_pfix = Pfix(pop_size, distance, radius);
        double right_flow = right_pfix * (nbr_sites - nbr_right) / nbr_sites;

        double left_pfix = Pfix(pop_size, distance, -radius);
        double left_flow = left_pfix * nbr_right / nbr_sites;

        double sub_flow = right_flow + left_flow;
        double draw = uniform_real_distribution<double>(0.0, sub_flow)(generator);

        if (draw < right_flow) {
            nbr_right++;
            distance += radius;
        } else {
            nbr_right--;
            distance -= radius;
        }

        if (!burn_in) {
            summary_stats["nbr_right"].add(nbr_right);
            summary_stats["distance"].add(distance);
            summary_stats["RatioSubFlow"].add(right_flow / left_flow);
            summary_stats["SubFlow"].add(sub_flow);
        }
    }

  public:
    explicit Sequence(unsigned nbr_sites, double const &pop_size, double const &radius,
        double const &peakness, double const &epistasis)
        : nbr_sites{nbr_sites},
          nbr_right{nbr_sites / 2},
          pop_size{pop_size},
          distance{0},
          radius{radius},
          peakness{peakness},
          epistasis{epistasis} {
        distance = radius * (nbr_right - (static_cast<double>(nbr_sites) / 2));
        summary_stats["nbr_right"] = SummaryStatistic();
        summary_stats["distance"] = SummaryStatistic();
        summary_stats["RatioSubFlow"] = SummaryStatistic();
        summary_stats["SubFlow"] = SummaryStatistic();
    };

    void run_substitutions(unsigned last, bool burn_in) {
        for (unsigned i = 0; i < last; ++i) { next_substitution(burn_in); }
    }

    void trace(string const &output) const {
        Trace trace;
        for (auto const &ss : summary_stats) {
            cout << ss.first + "-mean=" << ss.second.mean() << endl;
            cout << ss.first + "-std=" << ss.second.std() << endl;
            trace.add(ss.first + "-mean", ss.second.mean());
            trace.add(ss.first + "-std", ss.second.std());
            trace.add(ss.first + "-abs-mean", ss.second.abs_mean());
            trace.add(ss.first + "-abs-std", ss.second.abs_std());
        }
        trace.write_tsv(output);
    }
};

class CombinatorialStability : public OutputArgParse {
  public:
    explicit CombinatorialStability(CmdLine &cmd) : OutputArgParse(cmd) {}

    TCLAP::ValueArg<double> population_size{
        "", "population_size", "population size", false, 1e4, "double", cmd};
    TCLAP::ValueArg<unsigned> nbr_sites{
        "", "nbr_sites", "Number of sites", false, 300, "unsigned", cmd};
    TCLAP::ValueArg<double> radius{
        "", "radius", "radius of mutations", false, 0.1, "double", cmd};
    TCLAP::ValueArg<double> peakness{"", "peakness",
        "'alpha' parameter (peakness) of the fitness function (exp(-alpha*(d^beta))", false, 0.5,
        "double", cmd};
    TCLAP::ValueArg<double> epistasis{"", "epistasis",
        "'beta' parameter (epistasis) of fitness function (exp(-alpha*(d^beta))", false, 2.0,
        "double", cmd};
};

int main(int argc, char *argv[]) {
    CmdLine cmd{"BinaryStateGeo", ' ', "0.1"};
    CombinatorialStability args(cmd);
    cmd.parse(argc, argv);

    unsigned arg_seed = args.seed.getValue();
    string output_path{args.output_path.getValue()};
    cout << "Random generator seed: " << arg_seed << endl;
    generator.seed(arg_seed);

    unsigned nbr_sites{args.nbr_sites.getValue()};
    cout << "Sequence with " << nbr_sites << " sites." << endl;

    double pop_size{args.population_size.getValue()};
    cout << "Population of " << pop_size << " individuals." << endl;

    ofstream tsv;
    tsv.open(output_path + ".tsv");
    tsv << "seed\t" << arg_seed << endl;
    tsv << "nbr_sites\t" << nbr_sites << endl;
    tsv << "pop_size\t" << pop_size << endl;
    tsv.close();

    assert(pop_size > 1);
    Sequence seq(nbr_sites, pop_size, args.radius.getValue(), args.peakness.getValue(),
        args.epistasis.getValue());
    seq.run_substitutions(nbr_sites * 100, true);
    seq.run_substitutions(nbr_sites * 100, false);
    seq.trace(output_path);

    return 0;
}