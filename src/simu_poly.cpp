#include "argparse.hpp"
#include "fitness_additive.hpp"
#include "wright_fisher.hpp"

using namespace TCLAP;
using namespace std;

class SimuPolyArgParse : public SimuArgParse {
  public:
    explicit SimuPolyArgParse(CmdLine &cmd) : SimuArgParse(cmd) {}
    TCLAP::ValueArg<std::string> preferences_path{
        "f", "preferences", "input site-specific preferences path", true, "", "string", cmd};
    TCLAP::ValueArg<double> beta{
        "b", "beta", "Stringency parameter of the fitness profiles", false, 1.0, "double", cmd};
    TCLAP::ValueArg<u_long> pop_size{
        "n", "population_size", "Population size (at the root)", false, 500, "u_long", cmd};
    TCLAP::ValueArg<u_long> sample_size{
        "p", "sample_size", "Sample size (at the leaves)", false, 20, "u_long", cmd};
    TCLAP::ValueArg<double> noise_sigma{"", "noise_sigma",
        "The Ornstein–Uhlenbeck sigma (0<sigma) applied to Ne at each generation", false, 0.0,
        "double", cmd};
    TCLAP::ValueArg<double> noise_theta{"", "noise_theta",
        "The Ornstein–Uhlenbeck theta (0<=theta<1) applied to Ne at each generation", false, 0.9,
        "double", cmd};
};

int main(int argc, char *argv[]) {
    CmdLine cmd{"SimuPoly", ' ', "0.1"};
    SimuPolyArgParse args(cmd);
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
    double beta{args.beta.getValue()};
    assert(beta >= 0.0);
    bool branch_wise_correlation{args.branch_wise_correlation.getValue()};
    u_long pop_size{args.pop_size.getValue()};
    u_long sample_size{args.sample_size.getValue()};
    assert(sample_size <= pop_size);
    assert(sample_size > 0);
    double noise_sigma{args.noise_sigma.getValue()};
    assert(noise_sigma >= 0.0);
    double noise_theta{args.noise_theta.getValue()};
    assert(noise_theta < 1.0);

    u_long exon_size{args.exons.getValue()};
    SequenceAdditiveModel seq_fit_profiles(preferences_path, beta / (4 * pop_size), exon_size);
    u_long nbr_sites = seq_fit_profiles.nbr_sites();
    assert(0 <= exon_size and exon_size <= nbr_sites);

    Tree tree(newick_path);
    tree.set_root_age(root_age);

    u_long burn_in = 100 * pop_size;
    NucleotideRateMatrix nuc_matrix(nuc_matrix_path, mutation_rate_per_generation, true);

    LogMultivariate log_multivariate(pop_size, mutation_rate_per_generation, generation_time);
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
    parameters.add("#codonsites", nbr_sites);
    parameters.add("#nucleotidesites", nbr_sites * 3);
    parameters.add("preferences_beta", beta);
    parameters.add("nucleotide_matrix_path", output_path);
    parameters.add("mutation_rate_per_generation", mutation_rate_per_generation);
    nuc_matrix.add_to_trace(parameters);
    parameters.add("generation_time", generation_time);
    parameters.add("#generations_burn_in", burn_in);
    parameters.add("population_size", pop_size);
    parameters.add("sample_size", sample_size);
    parameters.add("exon_size", exon_size);
    parameters.add("fix_pop_size", fix_pop_size);
    parameters.add("fix_mut_rate", fix_mut_rate);
    parameters.add("fix_gen_time", fix_gen_time);
    parameters.add("noise_sigma", noise_sigma);
    parameters.add("noise_theta", noise_theta);
    correlation_matrix.add_to_trace(parameters);
    parameters.write_tsv(output_path + ".parameters");


    init_alignments(output_path, tree.nb_leaves(), nbr_sites * 3);
    Eigen::SelfAdjointEigenSolver<EMatrix> eigen_solver(correlation_matrix);
    EMatrix transform_matrix =
        eigen_solver.eigenvectors() * eigen_solver.eigenvalues().cwiseSqrt().asDiagonal();

    auto root_population = make_unique<Population>(seq_fit_profiles, sample_size, log_multivariate,
        nuc_matrix, transform_matrix, branch_wise_correlation, noise_sigma, noise_theta);
    root_population->burn_in(burn_in);

    Process simu_process(tree, root_population);
    simu_process.run(output_path);

    tracer_leaves.write_tsv(output_path + ".leaves");
    tracer_nodes.write_tsv(output_path + ".nodes");
    tracer_substitutions.write_tsv(output_path + ".substitutions");
    tracer_traits.write_tsv(output_path + ".traits");
    tracer_fossils.write_tsv(output_path + ".fossils");

    Population::timer_cout();

    cout << "Simulation computed." << endl;
    cout << nbr_sites * 3 * (mutation_rate_per_generation / generation_time) * tree.total_length()
         << " expected substitutions." << endl;
    cout << Exon::nbr_fixations << " simulated substitutions." << endl;
    cout << "Statistics summarized in: " << output_path + ".tsv" << endl;
    cout << "Fasta file in: " << output_path + ".fasta" << endl;
    cout << "Alignment (.ali) file in: " << output_path + ".ali" << endl;
    return 0;
}