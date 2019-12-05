#include "argparse.hpp"
#include "fitness_folding.hpp"
#include "process.hpp"

using namespace TCLAP;
using namespace std;

class SimuEvolArgParse : public SimuArgParse {
  public:
    explicit SimuEvolArgParse(CmdLine &cmd) : SimuArgParse(cmd) {}

    TCLAP::ValueArg<double> pop_size{
        "n", "population_size", "Population size (at the root)", false, 500, "u_long", cmd};
    TCLAP::ValueArg<u_long> nbr_grid_step{"d", "nbr_grid_step",
        "Number of intervals in which to discretize the brownian motion", false, 100, "u_long",
        cmd};
    TCLAP::ValueArg<string> pdb_folder{"", "pdb_folder", "The folder containing the .pdb files",
        false, "data/pdbfiles", "string", cmd};
    TCLAP::ValueArg<u_long> nbr_exons{
        "", "nbr_exons", "Number of exons in the protein", false, 5000, "u_long", cmd};
    TCLAP::ValueArg<double> cut_off{"", "cut_off",
        "The distance (in angstrom) to determine if 2 sites are in contact", false, 7.0, "double",
        cmd};
};

int main(int argc, char *argv[]) {
    CmdLine cmd{"SimuFold", ' ', "0.1"};
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
    u_long nbr_exons{args.nbr_exons.getValue()};
    u_long exon_size{args.exons.getValue()};
    StabilityModel seq_fitness(
        args.pdb_folder.getValue(), exon_size, nbr_exons, args.cut_off.getValue());
    u_long nbr_sites = seq_fitness.nbr_sites();
    assert(exon_size <= 300);
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

    Trace parameters;
    args.add_to_trace(parameters);
    tree.add_to_trace(parameters);
    parameters.add("pop_size", pop_size);
    parameters.add("cut_off", args.cut_off.getValue());
    parameters.add("pdb_folder", args.pdb_folder.getValue());
    parameters.add("#codonsites", nbr_sites);
    parameters.add("#nucleotidesites", nbr_sites * 3);
    nuc_matrix.add_to_trace(parameters);
    correlation_matrix.add_to_trace(parameters);
    parameters.write_tsv(output_path + ".parameters");

    // Initialisation of the root sequence
    init_alignments(output_path, tree.nb_leaves(), nbr_sites * 3);
    auto root_sequence = make_unique<Sequence>(
        seq_fitness, log_multivariate, nuc_matrix, transform_matrix, branch_wise_correlation);
    string equilibrium_filename = output_path + ".equilibrium";
    ifstream input_seq(equilibrium_filename + ".fasta");
    if (input_seq) {
        string dna_string;
        getline(input_seq, dna_string);
        root_sequence->set_from_dna_seq(dna_string);
        cout << "DNA Sequence at equilibrium found and starting from it." << endl;
    } else {
        cout << "AA Sequence found (in pdb) and starting from it." << endl;
        root_sequence->set_from_exon_aa_seq(seq_fitness.seq());
    }
    u_long burn_in_aa_changes = 5 * exon_size;
    u_long equilibrium_nbr_rounds = 15;
    cout << "DNA Sequence optimizing site marginals for " << equilibrium_nbr_rounds
         << " rounds, and running for " << burn_in_aa_changes << " amino-acid changes (per exon)"
         << endl;
    root_sequence->at_equilibrium(equilibrium_nbr_rounds, 1.0e2);
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