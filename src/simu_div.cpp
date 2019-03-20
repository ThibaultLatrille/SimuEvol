#include "argparse.hpp"
#include "codon.hpp"
#include "io.hpp"
#include "matrices.hpp"
#include "random.hpp"
#include "statistic.hpp"
#include "tree.hpp"

using namespace TCLAP;
using namespace std;

class Substitution {
  private:
    u_long site{0};
    char codon_from{0};
    char codon_to{0};

  public:
    double time_between_event;
    double non_syn_sub_flow;
    double non_syn_mut_flow;
    double syn_mut_flow;

    explicit Substitution(u_long _site, char codon_from, char codon_to, double time_between_event,
        double non_syn_sub_flow, double non_syn_mut_flow, double syn_mut_flow)
        : site{_site},
          codon_from{codon_from},
          codon_to{codon_to},
          time_between_event{time_between_event},
          non_syn_sub_flow{non_syn_sub_flow},
          non_syn_mut_flow{non_syn_mut_flow},
          syn_mut_flow{syn_mut_flow} {
        assert(site >= 0);
    }

    bool is_synonymous() const {
        return (codon_from != codon_to and
                Codon::codon_to_aa_array[codon_from] == Codon::codon_to_aa_array[codon_to]);
    }

    bool is_non_synonymous() const {
        return (codon_from != codon_to and
                Codon::codon_to_aa_array[codon_from] != Codon::codon_to_aa_array[codon_to]);
    }
};


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
    double next_substitution(NucleotideRateMatrix const &nuc_matrix, double beta, double time_start,
        double time_end, vector<Substitution> &substitutions_vec) {
        // Number of possible substitutions is 9 times the number of sites (3 substitutions for each
        // 3 possible positions).
        u_long nbr_substitutions{9 * nbr_sites};

        // Vector of substitution rates.
        vector<double> substitution_rates(nbr_substitutions, 0);

        // Sum of substitution rates.
        double total_substitution_rates{0.};

        double non_syn_sub_flow{0.0}, non_syn_mut_flow{0.0}, syn_mut_flow{0.0};

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
                if (Codon::codon_to_aa_array[codon_to] != 20) {
                    rate_substitution = nuc_matrix(n_from, n_to);
                    if (Codon::codon_to_aa_array[codon_from] !=
                        Codon::codon_to_aa_array[codon_to]) {
                        non_syn_mut_flow += rate_substitution;
                        rate_substitution *=
                            rate_fixation(aa_fitness_profiles[site], codon_from, codon_to, beta);
                        non_syn_sub_flow += rate_substitution;
                    } else {
                        syn_mut_flow += rate_substitution;
                    }
                }

                substitution_rates[9 * site + neighbor] = rate_substitution;
                // Increment the sum of substitution rates
                total_substitution_rates += rate_substitution;
            }
        }

        // Increment the time by drawing from an exponential distribution (mean equal to the inverse
        // sum of substitution rates).
        double time_draw = exponential_distribution<double>(total_substitution_rates)(generator);

        // Substitute the sequence if the time is positive, else there is no substitution but the
        // time left is set to 0.
        assert(total_substitution_rates != 0.0);
        if (time_start + time_draw <= time_end) {
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

            assert(non_syn_mut_flow < 10e10);
            assert(non_syn_sub_flow < 10e10);
            assert(syn_mut_flow < 10e10);

            substitutions_vec.emplace_back(Substitution(site, codom_from, codon_to, time_draw,
                non_syn_sub_flow, non_syn_mut_flow, syn_mut_flow));

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

        } else {
            substitutions_vec.emplace_back(Substitution(0, -1, -1, time_end - time_start,
                non_syn_sub_flow, non_syn_mut_flow, syn_mut_flow));

            time_start = time_end;
        }
        return time_start;
    }

    // Method computing all substitution event occurring during a given time-frame.
    // t: time during which substitution events occur (typically branch length).
    // This method is o(n²) where n is the number of sites, but can take into account epistatic
    // effects
    void run_substitutions(NucleotideRateMatrix const &nuc_matrix, double beta, double t_start,
        double t_end, vector<Substitution> &subs) {
        while (t_start < t_end) {
            t_start = next_substitution(nuc_matrix, beta, t_start, t_end, subs);
        }
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
        exons.reserve(dv.quot);
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
                    codon_frequencies(exon.aa_fitness_profiles[site], nuc_matrix, beta);
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

    void run_forward(double t_max, Tree const &tree) {
        double time_current = 0.0;
        Eigen::SelfAdjointEigenSolver<EMatrix> eigen_solver(cor_matrix);
        EMatrix transform =
            eigen_solver.eigenvectors() * eigen_solver.eigenvalues().cwiseSqrt().asDiagonal();
        EVector sampled_vector = EVector::Zero(cor_matrix.dimensions);

        double step_in_year = time_grid_step;

        while (time_current < t_max) {
            if (time_current + step_in_year >= t_max) { step_in_year = t_max - time_current; }
            double ratio = sqrt(step_in_year / tree.max_distance_to_root());
            for (int dim = 0; dim < cor_matrix.dimensions; dim++) {
                sampled_vector(dim) = normal_distrib(generator);
            }
            log_multivariate += ratio * (transform * sampled_vector);

            beta = log_multivariate.beta();
            generation_time = log_multivariate.generation_time();
            nuc_matrix.set_mutation_rate(log_multivariate.mu() / generation_time);

            for (auto &exon : exons) {
                exon.run_substitutions(nuc_matrix, beta, time_current, time_current + step_in_year, substitutions);
            }
            time_current += step_in_year;
            time_from_root += step_in_year;
        }
    }

    double predicted_dn_dn0(NucleotideRateMatrix const &rates, double const &relative_pop) const {
        double dn{0.}, dn0{0.};
        for (auto const &exon : exons) {
            double exon_dn{0}, exon_d0{0};
            tie(exon_dn, exon_d0) =
                ::predicted_dn_dn0(exon.aa_fitness_profiles, rates, relative_pop);
            dn += exon_dn;
            dn0 += exon_d0;
        }
        return dn / dn0;
    };

    double sequence_wise_predicted_dn_dn0(Sequence const &parent, NucleotideRateMatrix const &rates,
        double const &relative_pop) const {
        double dn{0.}, dn0{0.};

        for (size_t i = 0; i < exons.size(); i++) {
            double exon_dn{0}, exon_dn0{0};
            tie(exon_dn, exon_dn0) = ::flow_dn_dn0(
                exons[i].aa_fitness_profiles, parent.exons[i].codon_seq, rates, relative_pop);
            dn += exon_dn;
            dn0 += exon_dn0;
        }
        return dn / dn0;
    };

    double count_based_dn_dn0() const {
        double dn{0}, dn0{0};
        for (auto const &substitution : substitutions) {
            if (substitution.is_non_synonymous()) { dn++; }
            dn0 += substitution.non_syn_mut_flow * substitution.time_between_event;
        }
        return dn / dn0;
    }

    double flow_based_dn_dn0() const {
        double dn{0}, dn0{0};
        for (auto const &substitution : substitutions) {
            dn0 += substitution.non_syn_mut_flow * substitution.time_between_event;
            dn += substitution.non_syn_sub_flow * substitution.time_between_event;
        }
        return dn / dn0;
    }

    // Simulated omega from the substitutions
    double event_based_dn_ds() const {
        double dn{0}, ds{0};
        for (auto const &substitution : substitutions) {
            if (substitution.is_synonymous()) {
                ds += 1.0 / substitution.syn_mut_flow;
            } else if (substitution.is_non_synonymous()) {
                dn += 1.0 / substitution.non_syn_mut_flow;
            }
        }
        if (ds == .0) {
            cerr << "There is no synonymous substitutions, dN/dS can't be computed!" << endl;
            return .0;
        } else {
            return dn / ds;
        }
    }

    // Simulated omega from the substitutions
    double count_based_dn_ds() const {
        double dn{0}, ds{0}, dn0{0}, ds0{0};
        for (auto const &substitution : substitutions) {
            dn0 += substitution.non_syn_mut_flow * substitution.time_between_event;
            ds0 += substitution.syn_mut_flow * substitution.time_between_event;
            if (substitution.is_synonymous()) {
                ds++;
            } else if (substitution.is_non_synonymous()) {
                dn++;
            }
        }
        if (ds == .0) {
            cerr << "There is no synonymous substitutions, dN/dS can't be computed!" << endl;
            return .0;
        } else {
            return (dn * ds0) / (dn0 * ds);
        }
    }

    void node_trace(string const &output_filename, Tree::NodeIndex node, Tree &tree,
        Sequence const *parent) const {
        string node_name = tree.node_name(node);

        tree.set_tag(node, "population_size", to_string(beta));
        tree.set_tag(node, "generation_time_in_year", d_to_string(generation_time));
        tree.set_tag(node, "mutation_rate", d_to_string(nuc_matrix.mutation_rate));

        if (!tree.is_root(node)) {
            assert(parent != nullptr);
            double avg_beta = (parent->beta + beta) / 2;
            NucleotideRateMatrix avg_nuc_matrix{nuc_matrix};
            avg_nuc_matrix += parent->nuc_matrix;
            avg_nuc_matrix /= 2;

            tree.set_tag(node, "Branch_population_size", to_string(avg_beta));
            tree.set_tag(node, "Branch_generation_time_in_year",
                d_to_string((generation_time + parent->generation_time) / 2));
            tree.set_tag(node, "Branch_mutation_rate",
                d_to_string((nuc_matrix.mutation_rate + parent->nuc_matrix.mutation_rate) / 2));
            tree.set_tag(node, "Branch_LogNe", to_string(log10(avg_beta)));
            tree.set_tag(node, "Branch_dNdN0_predicted",
                d_to_string(predicted_dn_dn0(avg_nuc_matrix, avg_beta)));
            tree.set_tag(node, "Branch_dNdN0_sequence_wise_predicted",
                d_to_string(sequence_wise_predicted_dn_dn0(*parent, avg_nuc_matrix, avg_beta)));
            tree.set_tag(node, "Branch_dNdN0_count_based", d_to_string(count_based_dn_dn0()));
            tree.set_tag(node, "Branch_dNdN0_flow_based", d_to_string(flow_based_dn_dn0()));
            tree.set_tag(node, "Branch_dNdS_event_based", d_to_string(event_based_dn_ds()));
            tree.set_tag(node, "Branch_dNdS_count_based", d_to_string(count_based_dn_ds()));
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
                    codon_frequencies(exon.aa_fitness_profiles[site], nuc_matrix, beta);
                for (char codon_from{0}; codon_from < 64; codon_from++) {
                    // For all possible neighbors.
                    for (auto &neighbor : Codon::codon_to_neighbors_array[codon_from]) {
                        // Codon after mutation, Nucleotide original and Nucleotide after mutation.
                        char codon_to{0}, n_from{0}, n_to{0};
                        tie(codon_to, n_from, n_to) = neighbor;

                        if (Codon::codon_to_aa_array[codon_to] != 20) {
                            string key = "q_" + Codon::codon_string(codon_from) + "_" +
                                         Codon::codon_string(codon_to);
                            double fix = codon_freqs[codon_from] * nuc_matrix(n_from, n_to) *
                                         rate_fixation(exon.aa_fitness_profiles[site], codon_from,
                                             codon_to, beta);
                            trace.add(key, fix);
                        }
                    }
                }
            }
        }
        trace.write_tsv(output_filename + ".matrices");
    }

    void clear() { substitutions.clear(); }
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
        sequences[node]->clear();

        if (!tree.is_root(node)) {
            sequences[node]->run_forward(tree.node_length(node), tree);

            years_computed += tree.node_length(node);
            cout << years_computed << " years computed ("
                 << static_cast<int>(100 * years_computed / tree.total_length()) << "%)." << endl;

            sequences[node]->node_trace(output_filename, node, tree, sequences[tree.parent(node)]);
        }

        // If the node is internal, iterate through the direct children.
        for (auto &child : tree.children(node)) {
            sequences[child] = new Sequence(*sequences[node]);
            run_recursive(child, output_filename);
        }
    }

    // Recursively iterate through the subtree and count the number of substitutions.
    tuple<long, long> nbr_substitutions() {
        long nbr_non_syn{0}, nbr_syn{0};

        for (auto const seq : sequences) {
            long syn = count_if(seq->substitutions.begin(), seq->substitutions.end(),
                [](Substitution const &s) { return s.is_synonymous(); });
            nbr_syn += syn;
            nbr_non_syn += seq->substitutions.size() - syn;
        }

        return make_tuple(nbr_non_syn, nbr_syn);
    }
};

// Initialize static variables
double Process::years_computed = 0.0;

class SimuEvolArgParse : public SimuArgParse {
  public:
    explicit SimuEvolArgParse(CmdLine &cmd) : SimuArgParse(cmd) {}

    TCLAP::ValueArg<u_long> nbr_grid_step{"d", "nbr_grid_step",
        "Number of intervals in which discretize the brownian motion", false, 100, "u_long", cmd};
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
    u_long nbr_grid_step = args.nbr_grid_step.getValue();
    assert(nbr_grid_step > 0);
    time_grid_step = root_age / nbr_grid_step;
    double beta{args.beta.getValue()};
    assert(beta >= 0.0);
    double s{args.selection.getValue()};
    double p{args.shuffle_proba.getValue()};
    assert(p >= 0.0);
    assert(p <= 1.0);
    bool all_sites{args.shuffle_all.getValue()};

    vector<array<double, 20>> fitness_profiles = open_preferences(preferences_path, 1.0);
    u_long nbr_sites = fitness_profiles.size();
    u_long exon_size{args.exons.getValue()};
    if (exon_size == 0) { exon_size = nbr_sites; }
    assert(0 <= exon_size and exon_size <= nbr_sites);

    Tree tree(newick_path);
    tree.set_root_age(root_age);

    NucleotideRateMatrix nuc_matrix(nuc_matrix_path, mu / generation_time, true);

    LogMultivariate log_multivariate(beta, generation_time, mu);
    CorrelationMatrix correlation_matrix(correlation_path);

    Trace parameters;
    parameters.add("output_path", output_path);
    parameters.add("tree_path", newick_path);
    parameters.add("#tree_nodes", tree.nb_nodes());
    parameters.add("#tree_branches", tree.nb_branches());
    parameters.add("#tree_leaves", tree.nb_leaves());
    parameters.add("tree_ultrametric", tree.is_ultrametric());
    parameters.add("tree_min_distance_to_root_in_year", tree.min_distance_to_root());
    parameters.add("tree_max_distance_to_root_in_year", tree.max_distance_to_root());
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
    trace.write_tsv(output_path);

    cout << "Simulation computed." << endl;
    cout << nbr_sites * 3 * (mu / generation_time) * tree.total_length()
         << " expected substitutions." << endl;
    cout << nbr_synonymous + nbr_non_synonymous << " simulated substitutions." << endl;
    cout << "Statistics summarized in: " << output_path + ".tsv" << endl;
    cout << "Fasta file in: " << output_path + ".fasta" << endl;
    cout << "Alignment (.ali) file in: " << output_path + ".ali" << endl;
    return 0;
}