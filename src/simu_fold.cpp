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
    TCLAP::SwitchArg initialisation{
        "", "initialisation", "Initialize the amino-acid sequence.", cmd, false};
    TCLAP::ValueArg<double> cut_off{"", "cut_off",
        "The distance (in angstrom) to determine if 2 sites are in contact", false, 7.0, "double",
        cmd};
};

int main(int argc, char *argv[]) {
    CmdLine cmd{"SimuDiv", ' ', "0.1"};
    SimuEvolArgParse args(cmd);
    cmd.parse(argc, argv);

    u_long arg_seed = args.seed.getValue();
    cout << "Random generator seed: " << arg_seed << endl;
    generator.seed(arg_seed);

    string pdb_folder{args.pdb_folder.getValue()};
    string newick_path{args.newick_path.getValue()};
    string nuc_matrix_path{args.nuc_matrix_path.getValue()};
    string output_path{args.output_path.getValue()};
    string precision_path{args.precision_path.getValue()};
    bool fix_pop_size{args.fix_pop_size.getValue()};
    bool fix_mut_rate{args.fix_mut_rate.getValue()};
    bool fix_gen_time{args.fix_gen_time.getValue()};
    double mutation_rate_per_generation{args.mutation_rate_per_generation.getValue()};
    assert(mutation_rate_per_generation > 0.0);
    double root_age{args.root_age.getValue()};
    assert(root_age > 0.0);
    double generation_time{args.generation_time.getValue()};
    assert(generation_time > 0.0);
    assert(generation_time < root_age);
    u_long nbr_grid_step = args.nbr_grid_step.getValue();
    assert(nbr_grid_step > 0);
    time_grid_step = root_age / nbr_grid_step;
    double pop_size{args.pop_size.getValue()};
    assert(pop_size >= 0.0);
    bool branch_wise_correlation{args.branch_wise_correlation.getValue()};
    double cut_off{args.cut_off.getValue()};

    u_long nbr_exons{args.nbr_exons.getValue()};
    u_long exon_size{args.exons.getValue()};
    StabilityModel seq_fitness(pdb_folder, exon_size, nbr_exons, cut_off);
    u_long nbr_sites = seq_fitness.nbr_sites();
    assert(exon_size <= 300);
    assert(exon_size <= nbr_sites);

    Tree tree(newick_path);
    tree.set_root_age(root_age);

    NucleotideRateMatrix nuc_matrix(
        nuc_matrix_path, mutation_rate_per_generation / generation_time, true);

    LogMultivariate log_multivariate(pop_size, mutation_rate_per_generation, generation_time);
    CorrelationMatrix correlation_matrix(precision_path, fix_pop_size, fix_mut_rate, fix_gen_time);

    Trace parameters;
    parameters.add("seed", arg_seed);
    parameters.add("output_path", output_path);
    parameters.add("pdb_folder", pdb_folder);
    parameters.add("tree_path", newick_path);
    parameters.add("#tree_nodes", tree.nb_nodes());
    parameters.add("#tree_branches", tree.nb_branches());
    parameters.add("#tree_leaves", tree.nb_leaves());
    parameters.add("tree_ultrametric", tree.is_ultrametric());
    parameters.add("tree_min_distance_to_root_in_year", tree.min_distance_to_root());
    parameters.add("tree_max_distance_to_root_in_year", tree.max_distance_to_root());
    parameters.add("#codonsites", nbr_sites);
    parameters.add("#nucleotidesites", nbr_sites * 3);
    parameters.add("pop_size", pop_size);
    parameters.add("nucleotide_matrix_path", output_path);
    parameters.add("mutation_rate_per_generation", mutation_rate_per_generation);
    nuc_matrix.add_to_trace(parameters);
    parameters.add("generation_time", generation_time);
    parameters.add("exon_size", exon_size);
    parameters.add("fix_pop_size", fix_pop_size);
    parameters.add("fix_mut_rate", fix_mut_rate);
    parameters.add("fix_gen_time", fix_gen_time);
    correlation_matrix.add_to_trace(parameters);
    parameters.write_tsv(output_path + ".parameters");

    init_alignments(output_path, tree.nb_leaves(), nbr_sites * 3);
    Eigen::SelfAdjointEigenSolver<EMatrix> eigen_solver(correlation_matrix);
    EMatrix transform_matrix =
        eigen_solver.eigenvectors() * eigen_solver.eigenvalues().cwiseSqrt().asDiagonal();

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

    int burn_in_aa_changes = 5 * exon_size;
    int equilibrium_nbr_rounds = 15;
    cout << "DNA Sequence optimizing site marginals for " << equilibrium_nbr_rounds
         << " rounds, and running for " << burn_in_aa_changes << " amino-acid changes (per exon)"
         << endl;
    root_sequence->at_equilibrium(equilibrium_nbr_rounds, 1.0e2);
    root_sequence->burn_in(burn_in_aa_changes);


    // root_sequence.write_matrices(output_path);
    Process simu_process(tree, root_sequence);
    simu_process.run(output_path);

    long nbr_non_synonymous, nbr_synonymous;
    tie(nbr_non_synonymous, nbr_synonymous) = simu_process.nbr_substitutions();

    dnd0_event_tot /= tree.total_length();
    dnd0_count_tot /= tree.total_length();
    dnds_event_tot /= tree.total_length();
    dnds_count_tot /= tree.total_length();

    cout << "dnd0_event_tot is :" << dnd0_event_tot << endl;
    cout << "dnd0_count_tot is :" << dnd0_count_tot << endl;
    cout << "dnds_event_tot is :" << dnds_event_tot << endl;
    cout << "dnds_count_tot is :" << dnds_count_tot << endl;

    // .txt output
    Trace trace;
    trace.add("#substitutions", nbr_synonymous + nbr_non_synonymous);
    trace.add("#substitutions_per_site",
        static_cast<double>(nbr_synonymous + nbr_non_synonymous) / nbr_sites);
    trace.add("#synonymous_substitutions", nbr_synonymous);
    trace.add("#non_synonymous_substitutions", nbr_non_synonymous);
    trace.add("dnd0_event_tot", dnd0_event_tot);
    trace.add("dnd0_count_tot", dnd0_count_tot);
    trace.add("dnds_event_tot", dnds_event_tot);
    trace.add("dnds_count_tot", dnds_count_tot);
    trace.add("<Pfix>", pfix_tot / nbr_non_synonymous);
    trace.add("<s>", s_tot / nbr_non_synonymous);
    trace.add("<S=4Nes>", S_tot / nbr_non_synonymous);
    trace.add("<|s|>", s_abs_tot / nbr_non_synonymous);
    trace.add("<|S=4Nes|>", S_abs_tot / nbr_non_synonymous);
    trace.write_tsv(output_path);

    tracer_traits.write_tsv(output_path + ".traits");
    tracer_fossils.write_tsv(output_path + ".fossils");
    // tracer_substitutions.write_tsv(output_path + ".substitutions");
    tracer_sequences.write_tsv(output_path + ".sequences");

    cout << "Simulation computed." << endl;
    cout << nbr_sites * 3 * (mutation_rate_per_generation / generation_time) * tree.total_length()
         << " expected substitutions." << endl;
    cout << nbr_synonymous + nbr_non_synonymous << " simulated substitutions." << endl;
    cout << nbr_non_synonymous << " simulated non-synonymous substitutions." << endl;
    cout << "Statistics summarized in: " << output_path + ".tsv" << endl;
    cout << "Fasta file in: " << output_path + ".fasta" << endl;
    cout << "Alignment (.ali) file in: " << output_path + ".ali" << endl;
    return 0;
}