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

    u_long arg_seed = args.seed.getValue();
    cout << "Random generator seed: " << arg_seed << endl;
    generator.seed(arg_seed);

    string preferences_path{args.preferences_path.getValue()};
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
    double beta{args.beta.getValue()};
    assert(beta >= 0.0);
    bool branch_wise_correlation{args.branch_wise_correlation.getValue()};

    u_long exon_size{args.exons.getValue()};
    SequenceAdditiveModel seq_fit_profiles(preferences_path, 1.0 / 4, exon_size);
    u_long nbr_sites = seq_fit_profiles.nbr_sites();
    assert(0 <= exon_size and exon_size <= nbr_sites);

    Tree tree(newick_path);
    tree.set_root_age(root_age);

    NucleotideRateMatrix nuc_matrix(
        nuc_matrix_path, mutation_rate_per_generation / generation_time, true);

    LogMultivariate log_multivariate(beta, mutation_rate_per_generation, generation_time);
    CorrelationMatrix correlation_matrix(precision_path, fix_pop_size, fix_mut_rate, fix_gen_time);

    Trace parameters;
    parameters.add("seed", arg_seed);
    parameters.add("output_path", output_path);
    parameters.add("tree_path", newick_path);
    parameters.add("#tree_nodes", tree.nb_nodes());
    parameters.add("#tree_branches", tree.nb_branches());
    parameters.add("#tree_leaves", tree.nb_leaves());
    parameters.add("tree_ultrametric", tree.is_ultrametric());
    parameters.add("tree_min_distance_to_root_in_year", tree.min_distance_to_root());
    parameters.add("tree_max_distance_to_root_in_year", tree.max_distance_to_root());
    parameters.add("site_preferences_path", preferences_path);
    parameters.add("#codonsites", seq_fit_profiles.nbr_sites());
    parameters.add("#nucleotidesites", seq_fit_profiles.nbr_sites() * 3);
    parameters.add("preferences_beta", beta);
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
        seq_fit_profiles, log_multivariate, nuc_matrix, transform_matrix, branch_wise_correlation);
    root_sequence->at_equilibrium();
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
        static_cast<double>(nbr_synonymous + nbr_non_synonymous) / seq_fit_profiles.nbr_sites());
    trace.add("#synonymous_substitutions", nbr_synonymous);
    trace.add("#non_synonymous_substitutions", nbr_non_synonymous);
    trace.add("dnd0_event_tot", dnd0_event_tot);
    trace.add("dnd0_count_tot", dnd0_count_tot);
    trace.add("dnds_event_tot", dnds_event_tot);
    trace.add("dnds_count_tot", dnds_count_tot);
    for (auto const &ss : FitnessState::summary_stats) {
        trace.add(ss.first + "-mean", ss.second.mean());
        trace.add(ss.first + "-var", ss.second.variance());
        trace.add(ss.first + "-abs-mean", ss.second.abs_mean());
        trace.add(ss.first + "-var", ss.second.abs_variance());
    }
    trace.write_tsv(output_path);

    tracer_traits.write_tsv(output_path + ".traits");
    tracer_fossils.write_tsv(output_path + ".fossils");
    tracer_substitutions.write_tsv(output_path + ".substitutions");
    tracer_sequences.write_tsv(output_path + ".sequences");

    cout << "Simulation computed." << endl;
    cout << nbr_sites * 3 * (mutation_rate_per_generation / generation_time) * tree.total_length()
         << " expected substitutions." << endl;
    cout << nbr_synonymous + nbr_non_synonymous << " simulated substitutions." << endl;
    cout << "Statistics summarized in: " << output_path + ".tsv" << endl;
    cout << "Fasta file in: " << output_path + ".fasta" << endl;
    cout << "Alignment (.ali) file in: " << output_path + ".ali" << endl;
    return 0;
}