#include "argparse.hpp"
#include "codon.hpp"
#include "io.hpp"
#include "matrices.hpp"
#include "random.hpp"
#include "statistic.hpp"
#include "tree.hpp"

using namespace TCLAP;
using namespace std;

struct Substitution {
    u_long site{0};
    char codon_from{0};
    char codon_to{0};
    char nuc_from{0};
    char nuc_to{0};
    Matrix4x4 non_syn_opp_matrix;
    Matrix4x4 syn_opp_matrix;
};

bool is_synonymous(Substitution const &s) {
    return (Codon::codon_to_aa_array[s.codon_from] == Codon::codon_to_aa_array[s.codon_to]);
}

// Class representing DNA sequence.
class Exon {
  public:
    // The number of sites in the exon (each position is a codon, thus the DNA sequence is 3
    // times greater).
    u_long nbr_sites;
    u_long position;

    // The sequence of codons.
    vector<char> codon_seq;

    // The fitness profiles of amino-acids.
    vector<array<double, 20>> aa_fitness_profiles;

    // The selection coefficient against the current amino-acid
    double sel_coef;
    // The probability to randomize the fitness landscape.
    double proba_permutation;
    // Shuffle all sites or only the current one
    bool all_sites;

    // Constructor of Exon.
    // size: the size of the DNA sequence.
    explicit Exon(vector<array<double, 20>> const &fitness_profiles, const u_long &position,
        double const s, double const proba_permutation, bool const all_sites)
        : nbr_sites{static_cast<u_long>(fitness_profiles.size())},
          position{position},
          codon_seq(nbr_sites, 0),
          aa_fitness_profiles{fitness_profiles},
          sel_coef{s},
          proba_permutation{proba_permutation},
          all_sites{all_sites} {}

    // Method computing the next substitution event to occur, and the time for it to happen.
    // This method takes a time as input. If there is enough time given, the substitution event is
    // computed. This method also returns the time (given as input) decremented by the time needed
    // for the substitution event to occur. If there is not enough time given, no substitution event
    // is computed and the method returns 0. time_left: The time available for a substitution event
    // to occur.
    double next_substitution(NucleotideRateMatrix const &nuc_matrix, double &time_left,
        vector<Substitution> &substitutions_vec) {
        // Number of possible substitutions is 9 times the number of sites (3 substitutions for each
        // 3 possible positions).
        u_long nbr_substitutions{9 * nbr_sites};

        // Vector of substitution rates.
        vector<double> substitution_rates(nbr_substitutions, 0);

        // Sum of substitution rates.
        double total_substitution_rates{0.};

        Matrix4x4 non_syn_opp_matrix = Matrix4x4::Zero();
        Matrix4x4 syn_opp_matrix = Matrix4x4::Zero();

        // For all site of the sequence.
        for (u_long site{0}; site < nbr_sites; site++) {
            // Codon original before substitution.

            char codon_from = codon_seq[site];

            // Array of neighbors of the original codon (codons differing by only 1 mutation).
            array<tuple<char, char, char>, 9> neighbors =
                Codon::codon_to_neighbors_array[codon_from];

            // For all possible neighbors.
            for (char neighbor{0}; neighbor < 9; neighbor++) {
                // Codon after mutation, Nucleotide original and Nucleotide after mutation.
                char codon_to{0}, n_from{0}, n_to{0};
                tie(codon_to, n_from, n_to) = neighbors[neighbor];

                // Assign the substitution rate given by the method substitution rate.
                // Rate of substitution initialized to 0 (deleterious mutation)
                double rate_substitution{0.};

                // If the mutated amino-acid is a stop codon, the rate of fixation is 0.
                // Else, if the mutated and original amino-acids are non-synonymous, we compute the
                // rate of fixation. Note that, if the mutated and original amino-acids are
                // synonymous, the rate of fixation is 1.
                if (Codon::codon_to_aa_array[codon_to] == 20) {
                    rate_substitution = 0.;
                } else if (Codon::codon_to_aa_array[codon_from] !=
                           Codon::codon_to_aa_array[codon_to]) {
                    rate_substitution =
                        nuc_matrix(n_from, n_to) *
                        rate_fixation(aa_fitness_profiles[site], codon_from, codon_to);
                    non_syn_opp_matrix(n_from, n_to) += nuc_matrix(n_from, n_to);

                } else {
                    rate_substitution = nuc_matrix(n_from, n_to);
                    syn_opp_matrix(n_from, n_to) += nuc_matrix(n_from, n_to);
                }

                substitution_rates[9 * site + neighbor] = rate_substitution;
                // Increment the sum of substitution rates
                total_substitution_rates += rate_substitution;
            }
        }

        // Decrement the time by drawing from an exponential distribution (mean equal to the inverse
        // sum of substitution rates).
        time_left -= exponential_distribution<double>(total_substitution_rates)(generator);

        // Substitute the sequence if the time is positive, else there is no substitution but the
        // time left is set to 0.
        if (time_left > 0. and total_substitution_rates != 0.) {
            discrete_distribution<u_long> substitution_distr(
                substitution_rates.begin(), substitution_rates.end());

            u_long index = substitution_distr(generator);
            u_long site = index / 9;
            char codom_from = codon_seq[site];

            // Array of neighbors of the original codon (codons differing by only 1 mutation).
            array<tuple<char, char, char>, 9> neighbors =
                Codon::codon_to_neighbors_array[codom_from];

            // Codon after mutation, Nucleotide original and Nucleotide after mutation.
            char codon_to{0}, n_from{0}, n_to{0};

            tie(codon_to, n_from, n_to) = neighbors[index % 9];

            Substitution sub = {
                site, codom_from, codon_to, n_from, n_to, non_syn_opp_matrix, syn_opp_matrix};

            assert(non_syn_opp_matrix.sum() < 10e10);
            assert(syn_opp_matrix.sum() < 10e10);

            substitutions_vec.push_back(sub);

            if (sel_coef != 0.0) {
                char aa_from = Codon::codon_to_aa_array[codon_to];
                char aa_to = Codon::codon_to_aa_array[codon_to];
                aa_fitness_profiles[site][aa_from] += sel_coef;
                aa_fitness_profiles[site][aa_to] -= sel_coef;
            }

            if (proba_permutation != 0.0) {
                // Random shuffle of the fitness landscape
                uniform_real_distribution<double> unif_rand_proba(0, 1);
                double rand_uni = unif_rand_proba(generator);
                if (rand_uni < proba_permutation) {
                    if (all_sites) {
                        for (u_long site_shuffle{0}; site_shuffle < nbr_sites; site_shuffle++) {
                            shuffle(aa_fitness_profiles[site_shuffle].begin(),
                                aa_fitness_profiles[site_shuffle].end(), generator);
                        }
                    } else {
                        shuffle(aa_fitness_profiles[site].begin(), aa_fitness_profiles[site].end(),
                            generator);
                    }
                }
            }

            codon_seq[site] = codon_to;

        } else if (time_left < 0.) {
            time_left = 0.;
        }
        return time_left;
    }

    // Method computing all substitution event occurring during a given time-frame.
    // t: time during which substitution events occur (typically branch length).
    // This method is o(n²) where n is the number of sites, but can take into account epistatic
    // effects
    void run_substitutions(
        NucleotideRateMatrix const &nuc_matrix, double t, vector<Substitution> &subs) {
        while (t > 0) { t = next_substitution(nuc_matrix, t, subs); }
    }
};

static double time_grid_step;

class Sequence {
  public:
    // TimeElapsed
    double time_from_root{0};

    // Blocks
    vector<Exon> exons;

    LogMultivariate log_multivariate;
    double beta{0};
    double generation_time{0};
    NucleotideRateMatrix nuc_matrix;

    CorrelationMatrix const &cor_matrix;
    vector<Substitution> substitutions;

    Sequence(vector<array<double, 20>> const &fitness_profiles, LogMultivariate &log_multi,
        u_long exon_size, NucleotideRateMatrix nucleotide_matrix,
        CorrelationMatrix const &cor_matrix, double const s, double const proba_permutation,
        bool const all_sites)
        : exons{},
          log_multivariate{log_multi},
          beta{log_multivariate.beta()},
          generation_time{log_multivariate.generation_time()},
          nuc_matrix{move(nucleotide_matrix)},
          cor_matrix{cor_matrix} {
        auto dv = std::div(static_cast<long>(fitness_profiles.size()), exon_size);
        if (dv.rem != 0) { dv.quot++; }
        exons.reserve(static_cast<long>(dv.quot));
        for (int exon{0}; exon < dv.quot; exon++) {
            size_t begin_exon = exon * exon_size;
            size_t end_exon = min(begin_exon + exon_size, fitness_profiles.size());

            std::vector<array<double, 20>> exon_profiles(
                fitness_profiles.begin() + begin_exon, fitness_profiles.begin() + end_exon);

            exons.emplace_back(Exon(exon_profiles, begin_exon, s, proba_permutation, all_sites));
        }
        assert(nbr_sites() == fitness_profiles.size());
        cout << exons.size() << " exons created." << endl;
        if (dv.rem != 0) { cout << "Last exon is " << dv.rem << " sites long." << endl; }
    }

    u_long nbr_sites() const {
        u_long sites = 0;
        for (auto const &exon : exons) { sites += exon.nbr_sites; }
        return sites;
    }

    u_long nbr_nucleotides() const { return 3 * nbr_sites(); }

    // Set the the DNA sequence to the mutation-selection equilibrium.
    void at_equilibrium() {
        // For all site of the sequence.
        for (auto &exon : exons) {
            for (u_long site{0}; site < exon.nbr_sites; site++) {
                array<double, 64> codon_freqs =
                    codon_frequencies(exon.aa_fitness_profiles[site], nuc_matrix);
                discrete_distribution<char> freq_codon_distr(
                    codon_freqs.begin(), codon_freqs.end());
                char chosen_codon = freq_codon_distr(generator);
                exon.codon_seq[site] = chosen_codon;
                if (exon.sel_coef != 0.0) {
                    exon.aa_fitness_profiles[site][Codon::codon_to_aa_array[chosen_codon]] -=
                        exon.sel_coef;
                }
            }
        }
    }

    double predicted_omega() const {
        double sub_flow{0.}, mut_flow{0.};

        for (auto const &exon : exons) {
            double dn{0}, d0{0};
            tie(dn, d0) = predicted_dn_d0(exon.aa_fitness_profiles, nuc_matrix, 1.0);
            sub_flow += dn;
            mut_flow += d0;
        }
        return sub_flow / mut_flow;
    };

    void run_forward(double t_max, Tree const &tree) {
        double time = 0.0;
        Eigen::SelfAdjointEigenSolver<EMatrix> eigen_solver(cor_matrix);
        EMatrix transform =
            eigen_solver.eigenvectors() * eigen_solver.eigenvalues().cwiseSqrt().asDiagonal();
        EVector sampled_vector = EVector::Zero(cor_matrix.dimensions);

        double step_in_year = time_grid_step;

        while (time < t_max) {
            if (time + step_in_year >= t_max) { step_in_year = t_max - time; }
            double ratio = sqrt(step_in_year / tree.max_distance_to_root());
            for (int dim = 0; dim < cor_matrix.dimensions; dim++) {
                sampled_vector(dim) = normal_distrib(generator);
            }
            log_multivariate += ratio * (transform * sampled_vector);

            beta = log_multivariate.beta();
            generation_time = log_multivariate.generation_time();
            nuc_matrix.set_mutation_rate(log_multivariate.mu());

            double generations_per_step = step_in_year / generation_time;
            for (auto &exon : exons) {
                exon.run_substitutions(nuc_matrix, generations_per_step, substitutions);
            }
            time += step_in_year;
            time_from_root += step_in_year;
        }
    }

    double theoretical_dnds(double beta) const {
        double sub_flow{0.}, mut_flow{0.};

        for (auto const &exon : exons) {
            double dn{0}, d0{0};
            tie(dn, d0) = predicted_dn_d0(exon.aa_fitness_profiles, nuc_matrix, beta);
            sub_flow += dn;
            mut_flow += d0;
        }
        return sub_flow / mut_flow;
    };

    double flow_dnds(double beta) const {
        double sub_flow{0.}, mut_flow{0.};

        for (auto const &exon : exons) {
            double dn{0}, d0{0};
            tie(dn, d0) = flow_dn_d0(exon.aa_fitness_profiles, exon.codon_seq, nuc_matrix, beta);
            sub_flow += dn;
            mut_flow += d0;
        }
        return sub_flow / mut_flow;
    };

    void node_trace(string const &output_filename, Tree::NodeIndex node, Tree &tree,
        Sequence const *parent) const {
        string node_name = tree.node_name(node);

        tree.set_tag(node, "population_size", to_string(beta));
        tree.set_tag(node, "generation_time_in_year", d_to_string(generation_time));
        tree.set_tag(node, "mutation_rate_per_generation", d_to_string(nuc_matrix.mutation_rate));

        if (!tree.is_root(node)) {
            assert(parent != nullptr);
            tree.set_tag(node, "Branch_population_size", to_string((beta + parent->beta) / 2));
            tree.set_tag(node, "Branch_generation_time_in_year",
                d_to_string((generation_time + parent->generation_time) / 2));
            tree.set_tag(node, "Branch_mutation_rate_per_generation",
                d_to_string((nuc_matrix.mutation_rate + parent->nuc_matrix.mutation_rate) / 2));
            tree.set_tag(
                node, "Branch_dNdS_pred", d_to_string(theoretical_dnds((beta + parent->beta) / 2)));
            tree.set_tag(
                node, "Branch_dNdS_flow", d_to_string(flow_dnds((beta + parent->beta) / 2)));
        }

        if (tree.is_leaf(node)) {
            // If the node is a leaf, output the DNA sequence and name.
            write_sequence(output_filename, node_name, this->get_dna_str());
        }
    }

    // Method returning the DNA string corresponding to the codon sequence.
    string get_dna_str() const {
        string dna_str{};
        dna_str.reserve(nbr_nucleotides());

        // For each site of the sequence.
        for (auto const &exon : exons) {
            for (u_long site{0}; site < exon.nbr_sites; site++) {
                // Assert there is no stop in the sequence.
                assert(Codon::codon_to_aa_array[exon.codon_seq[site]] != 20);

                // Translate the site to a triplet of DNA nucleotides
                array<char, 3> triplet = Codon::codon_to_triplet_array[exon.codon_seq[site]];
                for (char position{0}; position < 3; position++) {
                    dna_str += Codon::nucleotides[triplet[position]];
                }
            }
        }
        return dna_str;  // return the DNA sequence as a string.
    }

    void write_matrices(string const &output_filename) {
        // For all site of the sequence.
        Trace trace;
        for (auto const &exon : exons) {
            for (u_long site{0}; site < exon.nbr_sites; site++) {
                // Codon original before substitution.
                trace.add("site", exon.position + site + 1);
                array<double, 64> codon_freqs =
                    codon_frequencies(exon.aa_fitness_profiles[site], nuc_matrix);
                for (char codon_from{0}; codon_from < 64; codon_from++) {
                    // For all possible neighbors.
                    for (auto &neighbor : Codon::codon_to_neighbors_array[codon_from]) {
                        // Codon after mutation, Nucleotide original and Nucleotide after mutation.
                        char codon_to{0}, n_from{0}, n_to{0};
                        tie(codon_to, n_from, n_to) = neighbor;

                        if (Codon::codon_to_aa_array[codon_to] != 20) {
                            string key = "q_" + Codon::codon_string(codon_from) + "_" +
                                         Codon::codon_string(codon_to);
                            double fix =
                                codon_freqs[codon_from] * nuc_matrix(n_from, n_to) *
                                rate_fixation(exon.aa_fitness_profiles[site], codon_from, codon_to);
                            trace.add(key, fix);
                        }
                    }
                }
            }
        }
        trace.write_tsv(output_filename + ".matrices");
    }
};

class Process {
  private:
    static double years_computed;
    Tree &tree;
    vector<Sequence *> sequences;  // Vector of sequence DNA.

  public:
    // Constructor
    Process(Tree &intree, Sequence &root_seq) : tree{intree}, sequences() {
        sequences.resize(tree.nb_nodes());
        sequences[tree.root()] = &root_seq;
    }

    void run(string const &output_filename) {
        sequences[tree.root()]->node_trace(output_filename, tree.root(), tree, nullptr);
        run_recursive(tree.root(), output_filename);
        ofstream nhx;
        nhx.open(output_filename + ".nhx");
        nhx << tree.as_string() << endl;
        nhx.close();
    }

    // Recursively iterate through the subtree.
    void run_recursive(Tree::NodeIndex node, string const &output_filename) {
        // Substitutions of the DNA sequence is generated.

        if (!tree.is_root(node)) {
            sequences[node]->run_forward(tree.node_length(node), tree);

            years_computed += tree.node_length(node);
            cout << years_computed << " years computed ("
                 << static_cast<int>(100 * years_computed / tree.total_length()) << "%)." << endl;

            sequences[node]->node_trace(output_filename, node, tree, sequences[tree.parent(node)]);
        }

        if (tree.is_leaf(node)) {
            // If the node is a leaf, output the DNA sequence and name.
            write_sequence(output_filename, tree.node_name(node), sequences[node]->get_dna_str());
        } else {
            // If the node is internal, iterate through the direct children.
            for (auto &child : tree.children(node)) {
                sequences[child] = new Sequence(*sequences[node]);
                run_recursive(child, output_filename);
            }
        }
    }

    // Recursively iterate through the subtree and count the number of substitutions.
    tuple<long, long> nbr_substitutions() {
        long nbr_non_syn{0}, nbr_syn{0};

        for (auto const seq : sequences) {
            long syn =
                count_if(seq->substitutions.begin(), seq->substitutions.end(), is_synonymous);
            nbr_syn += syn;
            nbr_non_syn += seq->substitutions.size() - syn;
        }

        return make_tuple(nbr_non_syn, nbr_syn);
    }

    // Recursively iterate through the subtree and count the number of substitutions.
    tuple<double, double> evolutionary_rates(string const &source, string const &target) {
        double dn{0}, ds{0};

        for (auto const seq : sequences) {
            for (auto const &substitution : seq->substitutions) {
                double non_syn_opp{0}, syn_opp{0};
                for (auto &dna_source : source) {
                    char nuc_source{Codon::nuc_to_index.at(dna_source)};
                    for (auto &dna_target : target) {
                        char nuc_target{Codon::nuc_to_index.at(dna_target)};
                        non_syn_opp += substitution.non_syn_opp_matrix(nuc_source, nuc_target);
                        syn_opp += substitution.syn_opp_matrix(nuc_source, nuc_target);
                    }
                }

                size_t source_find = source.find(Codon::nucleotides[substitution.nuc_from]);
                size_t target_find = target.find(Codon::nucleotides[substitution.nuc_to]);
                if (source_find != string::npos and target_find != string::npos) {
                    if (is_synonymous(substitution)) {
                        ds += 1. / syn_opp;
                    } else {
                        dn += 1. / non_syn_opp;
                    }
                }
            }
        }
        return make_tuple(dn, ds);
    }

    // Simulated omega from the substitutions
    double simulated_omega(string const &source, string const &target) {
        double dn{0}, ds{0};
        tie(dn, ds) = evolutionary_rates(source, target);
        if (ds == .0) {
            cerr << "There is no synonymous substitutions generated by the simulation, dN/dS can't "
                    "be computed!"
                 << endl;
            return .0;
        } else {
            return dn / ds;
        }
    }
};

// Initialize static variables
double Process::years_computed = 0.0;

class SimuEvolArgParse : public SimuArgParse {
  public:
    explicit SimuEvolArgParse(CmdLine &cmd) : SimuArgParse(cmd) {}

    TCLAP::ValueArg<double> selection{"", "selection",
        "Selection coefficient given the current amino-acid", false, 0.0, "double", cmd};
    TCLAP::ValueArg<double> shuffle_proba{"", "shuffle_proba",
        "Probability to randomize the fitness landscape (once a substitution occured)", false, 0.0,
        "double", cmd};
    SwitchArg shuffle_all{"", "shuffle_all",
        "All sites are affected by the random shuffling (instead of just the one where the "
        "substitution occured)",
        cmd, false};
};

int main(int argc, char *argv[]) {
    CmdLine cmd{"SimuDiv", ' ', "0.1"};
    SimuEvolArgParse args(cmd);
    cmd.parse(argc, argv);

    string preferences_path{args.preferences_path.getValue()};
    string newick_path{args.newick_path.getValue()};
    string nuc_matrix_path{args.nuc_matrix_path.getValue()};
    string output_path{args.output_path.getValue()};
    string correlation_path{args.correlation_path.getValue()};
    double mu{args.mu.getValue()};
    assert(mu > 0.0);
    double root_age{args.root_age.getValue()};
    assert(root_age > 0.0);
    double generation_time{args.generation_time.getValue()};
    assert(generation_time > 0.0);
    assert(generation_time < root_age);
    u_long nbr_grid_step = 100;
    time_grid_step = root_age / nbr_grid_step;
    double beta{args.beta.getValue()};
    assert(beta > 0.0);
    double s{args.selection.getValue()};
    double p{args.shuffle_proba.getValue()};
    assert(p >= 0.0);
    assert(p <= 1.0);
    bool all_sites{args.shuffle_all.getValue()};

    vector<array<double, 20>> fitness_profiles = open_preferences(preferences_path, beta);
    u_long nbr_sites = fitness_profiles.size();
    u_long exon_size{args.exons.getValue()};
    if (exon_size == 0) { exon_size = nbr_sites; }
    assert(0 <= exon_size and exon_size <= nbr_sites);

    Tree tree(newick_path);
    tree.set_root_age(root_age);

    NucleotideRateMatrix nuc_matrix(nuc_matrix_path, mu, true);

    LogMultivariate log_multivariate(beta, generation_time, mu);
    CorrelationMatrix correlation_matrix(correlation_path);

    Trace parameters;
    parameters.add("output_path", output_path);
    parameters.add("tree_path", newick_path);
    parameters.add("#tree_nodes", tree.nb_nodes());
    parameters.add("#tree_branches", tree.nb_branches());
    parameters.add("#tree_leaves", tree.nb_leaves());
    parameters.add("tree_ultrametric", tree.is_ultrametric());
    parameters.add("tree_min_distance_to_root", tree.min_distance_to_root());
    parameters.add("tree_max_distance_to_root", tree.max_distance_to_root());
    parameters.add("site_preferences_path", preferences_path);
    parameters.add("#codonsites", fitness_profiles.size());
    parameters.add("#nucleotidesites", fitness_profiles.size() * 3);
    parameters.add("preferences_beta", beta);
    parameters.add("nucleotide_matrix_path", output_path);
    parameters.add("mutation_rate_per_generation", mu);
    nuc_matrix.add_to_trace(parameters);
    parameters.add("generation_time_in_year", generation_time);
    parameters.add("exon_size", exon_size);
    correlation_matrix.add_to_trace(parameters);
    parameters.add("selection_coefficient_current_amino_acid", s);
    parameters.add("fitness_landscape_shuffle_probability", p);
    parameters.add("fitness_landscape_shuffle_all_sites", all_sites);
    parameters.write_tsv(output_path + ".parameters");

    init_alignments(output_path, tree.nb_leaves(), nbr_sites * 3);
    Sequence root_sequence(fitness_profiles, log_multivariate, exon_size, nuc_matrix,
        correlation_matrix, s, p, all_sites);
    root_sequence.at_equilibrium();
    root_sequence.write_matrices(output_path);

    Process simu_process(tree, root_sequence);
    simu_process.run(output_path);

    long nbr_non_synonymous, nbr_synonymous;
    tie(nbr_non_synonymous, nbr_synonymous) = simu_process.nbr_substitutions();

    // .txt output
    Trace trace;
    trace.add("#substitutions", nbr_synonymous + nbr_non_synonymous);
    trace.add("#substitutions_per_site",
        static_cast<double>(nbr_synonymous + nbr_non_synonymous) / fitness_profiles.size());
    trace.add("#synonymous_substitutions", nbr_synonymous);
    trace.add("#non_synonymous_substitutions", nbr_non_synonymous);
    trace.add("omega_predicted", root_sequence.predicted_omega());
    trace.add(
        "omega_simulated", simu_process.simulated_omega(Codon::nucleotides, Codon::nucleotides));
    trace.write_tsv(output_path);

    cout << "Simulation computed." << endl;
    cout << nbr_synonymous + nbr_non_synonymous << " substitutions." << endl;
    cout << "Statistics summarized in: " << output_path + ".tsv" << endl;
    cout << "Fasta file in: " << output_path + ".fasta" << endl;
    cout << "Alignment (.ali) file in: " << output_path + ".ali" << endl;
    return 0;
}