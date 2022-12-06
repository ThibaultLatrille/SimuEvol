#include "argparse.hpp"
#include "fitness_stability.hpp"
#include "wright_fisher.hpp"

using namespace TCLAP;
using namespace std;

int main(int argc, char *argv[]) {
    CmdLine cmd{"PolyStab", ' ', "0.1"};
    SimuPolyArgParse args(cmd);
    StabilityArgParse args_fitness(cmd);
    TreeArgParse args_tree(cmd);
    cmd.parse(argc, argv);

    u_long arg_seed = args.seed.getValue();
    cout << "Random generator seed: " << arg_seed << endl;
    generator.seed(arg_seed);

    string newick_path{args_tree.newick_path.getValue()};
    string nuc_matrix_path{args.nuc_matrix_path.getValue()};
    string output_path{args.output_path.getValue()};
    string precision_path{args_tree.precision_path.getValue()};
    bool fix_pop_size{args_tree.fix_pop_size.getValue()};
    bool fix_mut_rate{args_tree.fix_mut_rate.getValue()};
    bool fix_gen_time{args_tree.fix_gen_time.getValue()};
    bool fix_gBGC{args_tree.fix_gBGC.getValue()};
    double mutation_rate_per_generation{args.mutation_rate_per_generation.getValue()};
    assert(mutation_rate_per_generation > 0.0);
    double root_age{args_tree.root_age.getValue()};
    assert(root_age > 0.0);
    double generation_time{args.generation_time.getValue()};
    assert(generation_time > 0.0);
    assert(generation_time < root_age);
    double gBGC{args.gBGC.getValue()};
    bool branch_wise_correlation{args_tree.branch_wise_correlation.getValue()};
    u_long pop_size{static_cast<u_long>(args.pop_size.getValue())};
    u_long sample_size{args.sample_size.getValue()};
    assert(sample_size <= pop_size);
    assert(sample_size > 0);
    double noise_sigma{args.noise_sigma.getValue()};
    assert(noise_sigma >= 0.0);
    double noise_theta{args.noise_theta.getValue()};
    assert(noise_theta < 1.0);

    u_long exon_size{args.exons.getValue()};
    StabilityModel fitness_model(exon_size, args_fitness);
    u_long nbr_sites = fitness_model.nbr_sites();
    assert(exon_size <= nbr_sites);

    Tree tree(newick_path);
    tree.set_root_age(root_age);

    u_long burn_in = 100 * pop_size;
    NucleotideRateMatrix nuc_matrix(nuc_matrix_path, mutation_rate_per_generation, true);

    LogMultivariate log_multivariate(pop_size, mutation_rate_per_generation, generation_time, gBGC);
    CorrelationMatrix correlation_matrix(precision_path, fix_pop_size, fix_mut_rate, fix_gen_time, fix_gBGC);

    Trace parameters;
    args.add_to_trace(parameters);
    args_fitness.add_to_trace(parameters);
    args_tree.add_to_trace(parameters);
    tree.add_to_trace(parameters);
    parameters.add("#codonsites", nbr_sites);
    parameters.add("#nucleotidesites", nbr_sites * 3);
    parameters.add("#generations_burn_in", burn_in);
    nuc_matrix.add_to_trace(parameters);
    correlation_matrix.add_to_trace(parameters);
    parameters.write_tsv(output_path + ".parameters");


    init_alignments(output_path, tree.nb_leaves(), nbr_sites * 3);
    Eigen::SelfAdjointEigenSolver<EMatrix> eigen_solver(correlation_matrix);
    EMatrix transform_matrix =
        eigen_solver.eigenvectors() * eigen_solver.eigenvalues().cwiseSqrt().asDiagonal();

    auto root_population = make_unique<Population>(fitness_model, sample_size, log_multivariate,
        nuc_matrix, transform_matrix, branch_wise_correlation, noise_sigma, noise_theta);
    root_population->set_from_aa_seq(fitness_model.aa_seq());
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
    double expected_subs = static_cast<double>(nbr_sites) * 3 *
                           (mutation_rate_per_generation / generation_time) * tree.total_length();
    cout << expected_subs << " expected substitutions." << endl;
    cout << Exon::nbr_fixations << " simulated substitutions." << endl;
    cout << "Statistics summarized in: " << output_path + ".tsv" << endl;
    cout << "Fasta file in: " << output_path + ".fasta" << endl;
    cout << "Alignment (.ali) file in: " << output_path + ".ali" << endl;
    return 0;
}