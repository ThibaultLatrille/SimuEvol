#include "argparse.hpp"
#include "random.hpp"
#include "statistic.hpp"

using namespace TCLAP;
using namespace std;

double computePFolded(double const &deltaG) {
    double factor = exp(-deltaG);
    return (factor / (1 + factor));
}

double selection_coefficient_ddG(double const &deltaG, double const &deltaGMutant) {
    double fm = computePFolded(deltaGMutant), f = computePFolded(deltaG);
    return (fm - f) / f;
}

double Pfix(double const &pop_size, double const &deltaG, double const &ddG) {
    double S = 4 * selection_coefficient_ddG(deltaG, deltaG + ddG) * pop_size;
    if ((abs(S)) < 1e-4) {
        return 1 + S / 2;
    } else {
        return S / (1.0 - exp(-S));
    }
}

class Sequence {
    unsigned nbr_sites;
    unsigned nbr_states;
    unsigned distance{0};
    double dG{-10};
    double ddG{0.1};
    double pop_size;
    std::unordered_map<std::string, SummaryStatistic> summary_stats;

    void next_substitution(bool burn_in) {
        double deleterious_pfix = Pfix(pop_size, dG, ddG);
        double deleterious_flow = deleterious_pfix * nbr_states * (nbr_sites - distance) / (nbr_sites * nbr_states);

        double advantageous_pfix = Pfix(pop_size, dG, -ddG);
        double advantageous_flow = advantageous_pfix * distance / (nbr_sites * nbr_states);

        double sub_flow = deleterious_flow + advantageous_flow;
        double draw = uniform_real_distribution<double>(0.0, sub_flow)(generator);

        if (draw < deleterious_flow) {
            distance++;
            dG += ddG;
        } else {
            distance--;
            dG -= ddG;
        }

        if (!burn_in) {
            summary_stats["distance"].add(distance);
            summary_stats["dG"].add(dG);
            summary_stats["SubFlow"].add(sub_flow);
        }
    }

  public:
    explicit Sequence(unsigned nbr_sites, unsigned nbr_states, double const &pop_size)
        : nbr_sites{nbr_sites}, nbr_states{nbr_states}, pop_size{pop_size} {
        summary_stats["dG"] = SummaryStatistic();
        summary_stats["distance"] = SummaryStatistic();
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

    TCLAP::ValueArg<double> pop_size{"", "population_size", "population size", false, 1e4, "double", cmd};
    TCLAP::ValueArg<unsigned> nbr_sites{"", "nbr_sites", "Number of sites", false, 300, "unsigned", cmd};
    TCLAP::ValueArg<unsigned> nbr_states{"", "nbr_states", "Number of states", false, 10, "unsigned", cmd};
};

int main(int argc, char *argv[]) {
    CmdLine cmd{"SimuRelax", ' ', "0.1"};
    CombinatorialStability args(cmd);
    cmd.parse(argc, argv);

    unsigned arg_seed = args.seed.getValue();
    string output_path{args.output_path.getValue()};
    cout << "Random generator seed: " << arg_seed << endl;
    generator.seed(arg_seed);

    unsigned nbr_sites{args.nbr_sites.getValue()};
    unsigned nbr_states{args.nbr_states.getValue()};
    cout << "Sequence of " << nbr_states << " states with " << nbr_sites << " sites." << endl;

    double pop_size{args.pop_size.getValue()};
    cout << "Population of " << pop_size << " individuals." << endl;

    ofstream tsv;
    tsv.open(output_path + ".tsv");
    tsv << "seed\t" << arg_seed << endl;
    tsv << "nbr_sites\t" << nbr_sites << endl;
    tsv << "pop_size\t" << pop_size << endl;
    tsv.close();

    assert(pop_size > 1);
    Sequence seq(nbr_sites, nbr_states, pop_size);
    seq.run_substitutions(nbr_sites * 100, true);
    seq.run_substitutions(nbr_sites * 100, false);
    seq.trace(output_path);

    return 0;
}