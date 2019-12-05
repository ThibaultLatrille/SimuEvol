#include "argparse.hpp"
#include "fitness_geometric.hpp"
#include "process.hpp"

using namespace TCLAP;
using namespace std;

class SimuEvolArgParse : public SimuArgParse {
  public:
    explicit SimuEvolArgParse(CmdLine &cmd) : SimuArgParse(cmd) {}
    TCLAP::ValueArg<double> pop_size{
        "b", "population_size", "Effective population size", false, 1.0, "double", cmd};
    TCLAP::ValueArg<u_long> nbr_grid_step{"d", "nbr_grid_step",
        "Number of intervals in which discretize the brownian motion", false, 100, "u_long", cmd};
    TCLAP::ValueArg<u_long> complexity{
        "", "complexity", "Complexity of the fitness landscape", false, 3, "u_long", cmd};
    TCLAP::ValueArg<u_long> nbr_exons{"", "nbr_exons", "Number of exons", false, 30, "u_long", cmd};
};

int main(int argc, char *argv[]) {
    CmdLine cmd{"SimuGeo", ' ', "0.1"};
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
    u_long complexity = args.complexity.getValue();
    assert(complexity > 0);
    u_long nbr_exons{args.nbr_exons.getValue()};
    u_long exon_size{args.exons.getValue()};
    GeometricModel seq_fitness(exon_size, nbr_exons, complexity);
    u_long nbr_sites = seq_fitness.nbr_sites();
    assert(exon_size <= nbr_sites);

    // Process discretization
    double pop_size{args.pop_size.getValue()};
    assert(pop_size >= 0.0);
    u_long nbr_grid_step = args.nbr_grid_step.getValue();
    assert(nbr_grid_step > 0);
    time_grid_step = root_age / nbr_grid_step;
    LogMultivariate log_multivariate(pop_size, mutation_rate_per_generation, generation_time);
    CorrelationMatrix correlation_matrix(args.precision_path.getValue(),
        args.fix_pop_size.getValue(), args.fix_mut_rate.getValue(), args.fix_gen_time.getValue());
    Eigen::SelfAdjointEigenSolver<EMatrix> eigen_solver(correlation_matrix);
    EMatrix transform_matrix =
        eigen_solver.eigenvectors() * eigen_solver.eigenvalues().cwiseSqrt().asDiagonal();

    // Save the parameters of the simulation
    Trace parameters;
    args.add_to_trace(parameters);
    tree.add_to_trace(parameters);
    parameters.add("complexity", complexity);
    parameters.add("pop_size", pop_size);
    parameters.add("#codonsites", nbr_sites);
    parameters.add("#nucleotidesites", nbr_sites * 3);
    nuc_matrix.add_to_trace(parameters);
    correlation_matrix.add_to_trace(parameters);
    parameters.write_tsv(output_path + ".parameters");

    // Initialisation of the root sequence
    init_alignments(output_path, tree.nb_leaves(), nbr_sites * 3);
    auto root_sequence = make_unique<Sequence>(
        seq_fitness, log_multivariate, nuc_matrix, transform_matrix, branch_wise_correlation);
    u_long burn_in_aa_changes = exon_size * 2;
    u_long equilibrium_nbr_rounds = 15;
    cout << "DNA Sequence optimizing site marginals for " << equilibrium_nbr_rounds
         << " rounds, and running for " << burn_in_aa_changes << " amino-acid changes (per exon)"
         << endl;
    root_sequence->at_equilibrium(equilibrium_nbr_rounds);
    root_sequence->burn_in(burn_in_aa_changes);

    // Simulation along the tree
    Process simu_process(tree, root_sequence);
    simu_process.run(output_path);

    // Save the result summary statistics
    double expected_subs =
        nbr_sites * 3 * (mutation_rate_per_generation / generation_time) * tree.total_length();
    simu_process.summary(output_path, expected_subs);
    return 0;
}