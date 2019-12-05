#include "argparse.hpp"
#include "fitness_additive.hpp"
#include "process.hpp"

using namespace TCLAP;
using namespace std;

class SimuEvolArgParse : public SimuArgParse {
  public:
    explicit SimuEvolArgParse(CmdLine &cmd) : SimuArgParse(cmd) {}

    TCLAP::ValueArg<std::string> preferences_path{
        "f", "preferences", "input site-specific preferences path", true, "", "string", cmd};
    TCLAP::ValueArg<double> beta{
        "b", "beta", "Stringency parameter of the fitness profiles", false, 1.0, "double", cmd};
    TCLAP::ValueArg<u_long> nbr_grid_step{"d", "nbr_grid_step",
        "Number of intervals in which discretize the brownian motion", false, 100, "u_long", cmd};
};

int main(int argc, char *argv[]) {
    CmdLine cmd{"SimuDiv", ' ', "0.1"};
    SimuEvolArgParse args(cmd);
    cmd.parse(argc, argv);

    generator.seed(args.seed.getValue());
    string output_path{args.output_path.getValue()};

    // Nucleotide matrix
    double mutation_rate_per_generation{args.mutation_rate_per_generation.getValue()};
    double generation_time{args.generation_time.getValue()};
    NucleotideRateMatrix nuc_matrix(
        args.nuc_matrix_path.getValue(), mutation_rate_per_generation / generation_time, true);

    // Tree
    double root_age{args.root_age.getValue()};
    bool branch_wise_correlation{args.branch_wise_correlation.getValue()};
    Tree tree(args.newick_path.getValue());
    tree.set_root_age(root_age);

    // Fitness model
    string preferences_path{args.preferences_path.getValue()};
    u_long exon_size{args.exons.getValue()};
    SequenceAdditiveModel seq_fit_profiles(preferences_path, 1.0 / 4, exon_size);
    u_long nbr_sites = seq_fit_profiles.nbr_sites();
    assert(0 <= exon_size and exon_size <= nbr_sites);

    // Process discretization
    double beta{args.beta.getValue()};
    assert(beta >= 0.0);
    u_long nbr_grid_step = args.nbr_grid_step.getValue();
    assert(nbr_grid_step > 0);
    time_grid_step = root_age / nbr_grid_step;
    LogMultivariate log_multivariate(beta, mutation_rate_per_generation, generation_time);
    CorrelationMatrix correlation_matrix(args.precision_path.getValue(),
        args.fix_pop_size.getValue(), args.fix_mut_rate.getValue(), args.fix_gen_time.getValue());
    Eigen::SelfAdjointEigenSolver<EMatrix> eigen_solver(correlation_matrix);
    EMatrix transform_matrix =
        eigen_solver.eigenvectors() * eigen_solver.eigenvalues().cwiseSqrt().asDiagonal();

    // Save the parameters of the simulation
    Trace parameters;
    args.add_to_trace(parameters);
    tree.add_to_trace(parameters);
    parameters.add("site_preferences_path", preferences_path);
    parameters.add("preferences_beta", beta);
    parameters.add("#exons", seq_fit_profiles.nbr_exons());
    parameters.add("#codonsites", nbr_sites);
    parameters.add("#nucleotidesites", nbr_sites * 3);
    nuc_matrix.add_to_trace(parameters);
    correlation_matrix.add_to_trace(parameters);
    parameters.write_tsv(output_path + ".parameters");

    // Initialisation of the root sequence
    init_alignments(output_path, tree.nb_leaves(), nbr_sites * 3);
    auto root_sequence = make_unique<Sequence>(
        seq_fit_profiles, log_multivariate, nuc_matrix, transform_matrix, branch_wise_correlation);
    root_sequence->at_equilibrium();

    // Simulation along the tree
    Process simu_process(tree, root_sequence);
    simu_process.run(output_path);

    // Save the result summary statistics
    double expected_subs =
        nbr_sites * 3 * (mutation_rate_per_generation / generation_time) * tree.total_length();
    simu_process.summary(output_path, expected_subs);

    return 0;
}