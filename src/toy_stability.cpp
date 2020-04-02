#include "argparse.hpp"
#include "random.hpp"
#include "statistic.hpp"

using namespace TCLAP;
using namespace std;
uniform_real_distribution<double> unif(0.0, 1.0);

class Sequence {
  private:
    u_long nbr_states;
    double delta_x;
    u_long exon_size;
    double x{0.3};
    double pop_size;
    double alpha;
    double gamma;
    double beta;
    std::unordered_map<std::string, SummaryStatistic> summary_stats = {};

  public:
    explicit Sequence(u_long exon_size, u_long nbr_states, double const &pop_size,
        double const &alpha, double const &gamma, double const &beta)
        : nbr_states{nbr_states},
          delta_x{1.0 / exon_size},
          exon_size{exon_size},
          pop_size{pop_size},
          alpha{alpha},
          gamma{gamma},
          beta{beta} {};

    void run_substitutions(u_long last, bool burn_in) {
        for (u_long i = 0; i < last; ++i) { next_substitution(burn_in); }
    }

    void trace(string const &output) const {
        Trace trace;
        for (auto const &ss : summary_stats) {
            cout << ss.first + "-mean=" << ss.second.mean() << endl;
            cout << ss.first + "-std=" << ss.second.std() << endl;
            trace.add(ss.first + "-mean", ss.second.mean());
            trace.add(ss.first + "-std", ss.second.std());
        }
        trace.write_tsv(output);
    }

  private:
    double computeFitness(double const &new_x) {
        double factor = exp(-beta * (alpha + gamma * exon_size * new_x));
        return (factor / (1 + factor));
    }

    double selection_coefficient(double const &deltax) {
        double fm = computeFitness(x + deltax), f = computeFitness(x);
        return (fm - f) / f;
    }

    double Pfix(double const &s) {
        double S = 4 * s * pop_size;
        if ((abs(S)) < 1e-4) {
            return 1 + S / 2;
        } else {
            return S / (1.0 - exp(-S));
        }
    }

    void next_substitution(bool burn_in) {
        double neutral_flow = (nbr_states - 2) * x / (nbr_states - 1);

        double del_s = selection_coefficient(delta_x);
        double del_pfix = Pfix(del_s);
        double del_flow = (x == 1.0 ? 0.0 : del_pfix * (1 - x));

        double adv_s = selection_coefficient(-delta_x);
        double adv_pfix = Pfix(adv_s);
        double adv_flow = (x == 0 ? 0.0 : adv_pfix * x / (nbr_states - 1));

        bool adv_drawn = unif(generator) <= adv_flow / (del_flow + adv_flow);
        if (adv_drawn) {
            x -= delta_x;
        } else {
            x += delta_x;
        }
        if (x < 0.0) { x = -x; }
        if (x > 1.0) { x = 2.0 - x; }
        assert(x > 0.0 and x < 1.0);

        if (burn_in) { return; }

        summary_stats["p_adv_sub_count"].add(adv_drawn);
        summary_stats["p_adv_sub_flow"].add(adv_flow / (del_flow + adv_flow));
        summary_stats["s_sub"].add(adv_s * adv_drawn + (1 - adv_drawn) * del_s);
        summary_stats["x"].add(x);
        summary_stats["dG"].add(alpha + gamma * x);
        summary_stats["dNdN0"].add(del_flow + neutral_flow + adv_flow);
    }
};

class CombinatorialStability : public OutputArgParse {
  public:
    explicit CombinatorialStability(CmdLine &cmd) : OutputArgParse(cmd) {}

    TCLAP::ValueArg<double> population_size{
        "", "population_size", "population size", false, 1e4, "double", cmd};
    TCLAP::ValueArg<u_long> exon_size{
        "", "exon_size", "Number of sites", false, 300, "u_long", cmd};
    TCLAP::ValueArg<u_long> nbr_states{
        "", "nbr_states", "Number of states", false, 2, "u_long", cmd};
    TCLAP::ValueArg<double> alpha{
        "", "alpha", "ΔG minimum for the optimal sequence", false, -118, "double", cmd};
    TCLAP::ValueArg<double> gamma{
        "", "gamma", "Maximum increase in ΔG for the pesimal sequence", false, 300, "double", cmd};
    TCLAP::ValueArg<double> beta{"", "beta", "beta=1/kT in mol/kcal", false, 1.686, "double", cmd};

    void add_to_trace(Trace &trace) {
        OutputArgParse::add_2_trace(trace);
        trace.add("population_size", population_size.getValue());
        trace.add("#exon_size", exon_size.getValue());
        trace.add("#nbr_states", nbr_states.getValue());
        trace.add("alpha", alpha.getValue());
        trace.add("gamma", gamma.getValue());
        trace.add("beta", beta.getValue());
    }
};

int main(int argc, char *argv[]) {
    CmdLine cmd{"BinaryStateStab", ' ', "0.1"};
    CombinatorialStability args(cmd);
    cmd.parse(argc, argv);

    string output_path{args.output_path.getValue()};
    generator.seed(args.seed.getValue());

    Trace parameters;
    args.add_to_trace(parameters);
    parameters.write_tsv(output_path + ".parameters");

    u_long exon_size{args.exon_size.getValue()};
    Sequence seq(exon_size, args.nbr_states.getValue(), args.population_size.getValue(),
        args.alpha.getValue(), args.gamma.getValue(), args.beta.getValue());
    seq.run_substitutions(exon_size * 1000, true);
    seq.run_substitutions(exon_size * 1000, false);
    seq.trace(output_path);

    return 0;
}