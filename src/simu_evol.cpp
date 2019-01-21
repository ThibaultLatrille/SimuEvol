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
    unsigned site{0};
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

double rate_fixation(array<double, 20> const &preferences, char codon_from, char codon_to) {
    double rate_fix = 1.0;
    // Selective strength between the mutated and original amino-acids.
    double s{0.};
    s = preferences[Codon::codon_to_aa_array[codon_to]];
    s -= preferences[Codon::codon_to_aa_array[codon_from]];
    // If the selective strength is 0, the rate of fixation is neutral.
    // Else, the rate of fixation is computed using population genetic formulas
    // (Kimura).
    if (fabs(s) > Codon::epsilon) {
        // The substitution rate is the mutation rate multiplied by the rate of
        // fixation.
        rate_fix = s / (1 - exp(-s));
    }
    return rate_fix;
}

// Class representing DNA sequence.
class Sequence_dna {
  private:
    // The number of sites in the sequence (each position is a codon, thus the DNA sequence is 3
    // times greater).
    unsigned const nbr_sites;

    // The sequence of codons.
    vector<char> codon_seq;

    // The matrix of mutation rates between nucleotides.
    NucleotideRateMatrix const &mutation_rate_matrix;

    // The fitness profiles of amino-acids.
    vector<array<double, 20>> aa_fitness_profiles;

    // The selection coefficient against the current amino-acid
    double sel_coef;

    // The probability to randomize the fitness landscape.
    double proba_permutation;

    // Shuffle all sites or only the current one
    bool all_sites;

  public:
    // Constructor of Sequence_dna.
    // size: the size of the DNA sequence.
    explicit Sequence_dna(NucleotideRateMatrix const &mutation_rate,
        vector<array<double, 20>> const &fitness_profiles, double const s,
        double const proba_permutation, bool const all_sites)
        : nbr_sites{static_cast<unsigned>(fitness_profiles.size())},
          codon_seq(nbr_sites, 0),
          mutation_rate_matrix{mutation_rate},
          aa_fitness_profiles{fitness_profiles},
          sel_coef{s},
          proba_permutation{proba_permutation},
          all_sites{all_sites} {
        ;
    }

    // Method translating the codon sequence to amino-acid sequence.
    string translate() const {
        string protein(nbr_sites, ' ');

        // For each site
        for (unsigned site{0}; site < nbr_sites; site++) {
            // Use the codon_to_aa_array to translate site to amino-acid.
            protein[site] = Codon::amino_acids[Codon::codon_to_aa_array[codon_seq[site]]];
        }

        return protein;  // return the amino-acid sequence as a string.
    }

    // Method computing the next substitution event to occur, and the time for it to happen.
    // This method takes a time as input. If there is enough time given, the substitution event is
    // computed. This method also returns the time (given as input) decremented by the time needed
    // for the substitution event to occur. If there is not enough time given, no substitution event
    // is computed and the method returns 0. time_left: The time available for a substitution event
    // to occur.
    double next_substitution(double &time_left, vector<Substitution> &substitutions_vec) {
        // Number of possible substitutions is 9 times the number of sites (3 substitutions for each
        // 3 possible positions).
        unsigned nbr_substitutions{9 * nbr_sites};

        // Vector of substitution rates.
        vector<double> substitution_rates(nbr_substitutions, 0);

        // Sum of substitution rates.
        double total_substitution_rates{0.};

        Matrix4x4 non_syn_opp_matrix = Matrix4x4::Zero();
        Matrix4x4 syn_opp_matrix = Matrix4x4::Zero();

        // For all site of the sequence.
        for (unsigned site{0}; site < nbr_sites; site++) {
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
                        mutation_rate_matrix(n_from, n_to) *
                        rate_fixation(aa_fitness_profiles[site], codon_from, codon_to);
                    non_syn_opp_matrix(n_from, n_to) += mutation_rate_matrix(n_from, n_to);

                } else {
                    rate_substitution = mutation_rate_matrix(n_from, n_to);
                    syn_opp_matrix(n_from, n_to) += mutation_rate_matrix(n_from, n_to);
                }

                substitution_rates[9 * site + neighbor] = rate_substitution;
                // Increment the sum of substitution rates
                total_substitution_rates += rate_substitution;
            }
        }

        // Decrement the time by inverse of the sum of substitution rates.
        time_left -= 1. / total_substitution_rates;

        // Substitute the sequence if the time is positive, else there is no substitution but the
        // time left is set to 0.
        if (time_left > 0. and total_substitution_rates != 0.) {
            discrete_distribution<unsigned> substitution_distr(
                substitution_rates.begin(), substitution_rates.end());

            unsigned index = substitution_distr(generator);
            unsigned site = index / 9;
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
                        for (unsigned site_shuffle{0}; site_shuffle < nbr_sites; site_shuffle++) {
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
    void run_substitutions(double t, vector<Substitution> &subs) {
        while (t > 0) { t = next_substitution(t, subs); }
    }


    // Method computing the equilibrium frequencies for one site.
    // aa_fitness_profil: The amino-acid fitness profil of the given site.
    array<double, 64> codon_frequencies(array<double, 20> const &aa_fitness_profil) {
        Vector4x1 nuc_frequencies = mutation_rate_matrix.nuc_frequencies;

        array<double, 64> codon_frequencies{0};

        // Initialize the total sum of equilibrium frequency at 0.
        double total_frequencies{0.};

        // For each site of the vector of the site frequencies.
        for (char codon{0}; codon < 64; codon++) {
            double codon_freq{1.};

            // For all nucleotides in the codon
            for (auto &nuc : Codon::codon_to_triplet_array[codon]) {
                codon_freq *= nuc_frequencies(nuc);
            }

            if (Codon::codon_to_aa_array[codon] != 20) {
                codon_frequencies[codon] =
                    codon_freq * exp(aa_fitness_profil[Codon::codon_to_aa_array[codon]]);
            } else {
                codon_frequencies[codon] = 0.;
            }

            // Increment the total sum of equilibrium frequencies.
            total_frequencies += codon_frequencies[codon];
        }

        // Normalize the vector of equilibrium frequencies.
        for (char codon{0}; codon < 64; codon++) { codon_frequencies[codon] /= total_frequencies; }

        return codon_frequencies;
    };

    // Set the the DNA sequence to the mutation-selection equilibrium.
    void at_equilibrium() {
        // For all site of the sequence.
        for (unsigned site{0}; site < nbr_sites; site++) {
            array<double, 64> codon_freqs = codon_frequencies(aa_fitness_profiles[site]);
            discrete_distribution<char> freq_codon_distr(codon_freqs.begin(), codon_freqs.end());
            char chosen_codon = freq_codon_distr(generator);
            codon_seq[site] = chosen_codon;
            if (sel_coef != 0.0) {
                aa_fitness_profiles[site][Codon::codon_to_aa_array[chosen_codon]] -= sel_coef;
            }
        }
    }

    // Theoretical computation of the predicted omega
    double predicted_omega(string const &source, string const &target, bool averaged = true) {
        vector<double> omega_per_site(nbr_sites, 0.);
        vector<double> dn_per_site(nbr_sites, 0.);
        vector<double> ds_per_site(nbr_sites, 0.);

        // For all site of the sequence.
        for (unsigned site{0}; site < nbr_sites; site++) {
            // Codon original before substitution.

            array<double, 64> codon_freqs = codon_frequencies(aa_fitness_profiles[site]);
            double dn{0.}, d0{0.};
            for (char codon_from{0}; codon_from < 64; codon_from++) {
                if (Codon::codon_to_aa_array[codon_from] != 20) {
                    // For all possible neighbors.
                    for (auto &neighbor : Codon::codon_to_neighbors_array[codon_from]) {
                        // Codon after mutation, Nucleotide original and Nucleotide after mutation.
                        char codon_to{0}, n_from{0}, n_to{0};
                        tie(codon_to, n_from, n_to) = neighbor;

                        size_t source_find = source.find(Codon::nucleotides[n_from]);
                        size_t target_find = target.find(Codon::nucleotides[n_to]);
                        if (source_find != string::npos and target_find != string::npos) {
                            // If the mutated amino-acid is a stop codon, the rate of fixation is 0.
                            // Else, if the mutated and original amino-acids are non-synonymous, we
                            // compute the rate of fixation. Note that, if the mutated and original
                            // amino-acids are synonymous, the rate of fixation is 1.
                            if (Codon::codon_to_aa_array[codon_to] != 20 and
                                Codon::codon_to_aa_array[codon_from] !=
                                    Codon::codon_to_aa_array[codon_to]) {
                                dn +=
                                    codon_freqs[codon_from] * mutation_rate_matrix(n_from, n_to) *
                                    rate_fixation(aa_fitness_profiles[site], codon_from, codon_to);
                                d0 += codon_freqs[codon_from] * mutation_rate_matrix(n_from, n_to);
                            }
                        }
                    }
                }
            }
            dn_per_site[site] = dn;
            ds_per_site[site] = d0;
            omega_per_site[site] = dn / d0;
        }
        if (averaged) {
            return avg(omega_per_site);
        } else {
            return sum(dn_per_site) / sum(ds_per_site);
        }
    }

    // Method returning the DNA string corresponding to the codon sequence.
    string get_dna_str() const {
        // The DNA string is 3 times larger than the codon sequence.
        string dna_str(nbr_sites * 3, ' ');

        // For each site of the sequence.
        for (unsigned site{0}; site < nbr_sites; site++) {
            // Assert there is no stop in the sequence.
            assert(Codon::codon_to_aa_array[codon_seq[site]] != 20);

            // Translate the site to a triplet of DNA nucleotides
            array<char, 3> triplet = Codon::codon_to_triplet_array[codon_seq[site]];
            for (char position{0}; position < 3; position++) {
                dna_str[3 * site + position] = Codon::nucleotides[triplet[position]];
            }
        }
        return dna_str;  // return the DNA sequence as a string.
    }

    void write_matrices(string const &output_filename) {
        // For all site of the sequence.
        Trace trace;
        for (unsigned site{0}; site < nbr_sites; site++) {
            // Codon original before substitution.
            trace.add("site", site + 1);
            array<double, 64> codon_freqs = codon_frequencies(aa_fitness_profiles[site]);
            for (char codon_from{0}; codon_from < 64; codon_from++) {
                // For all possible neighbors.
                for (auto &neighbor : Codon::codon_to_neighbors_array[codon_from]) {
                    // Codon after mutation, Nucleotide original and Nucleotide after mutation.
                    char codon_to{0}, n_from{0}, n_to{0};
                    tie(codon_to, n_from, n_to) = neighbor;

                    if (Codon::codon_to_aa_array[codon_to] != 20) {
                        string key = "q_" + Codon::codon_string(codon_from) + "_" +
                                     Codon::codon_string(codon_to);
                        double fix = codon_freqs[codon_from] * mutation_rate_matrix(n_from, n_to) *
                                     rate_fixation(aa_fitness_profiles[site], codon_from, codon_to);
                        trace.add(key, fix);
                    }
                }
            }
        }
        trace.write_tsv(output_filename + ".matrices");
    }
};


class Process {
  private:
    const Tree &tree;
    vector<Sequence_dna *> sequences;    // Vector of sequence DNA.
    vector<Substitution> substitutions;  // The number of substitutions.

  public:
    // Constructor
    Process(const Tree &intree, Sequence_dna &root_seq) : tree{intree}, sequences() {
        sequences.resize(tree.nb_nodes());
        sequences[tree.root()] = &root_seq;
    }

    // Recursively iterate through the subtree and count the number of substitutions.
    tuple<long, long> nbr_substitutions() {
        long nbr_non_syn{0}, nbr_syn{0};

        nbr_syn = count_if(substitutions.begin(), substitutions.end(), is_synonymous);
        nbr_non_syn = substitutions.size() - nbr_syn;

        return make_tuple(nbr_non_syn, nbr_syn);
    }

    void run(string const &output_filename) { run_recursive(tree.root(), output_filename); }

    // Recursively iterate through the subtree.
    void run_recursive(Tree::NodeIndex node, string const &output_filename) {
        // Substitutions of the DNA sequence is generated.
        sequences[node]->run_substitutions(tree.node_length(node), substitutions);

        if (tree.is_leaf(node)) {
            // If the node is a leaf, output the DNA sequence and name.
            write_sequence(output_filename, tree.node_name(node), sequences[node]->get_dna_str());
        } else {
            // If the node is internal, iterate through the direct children.
            for (auto &child : tree.children(node)) {
                sequences[child] = new Sequence_dna(*sequences[node]);
                run_recursive(child, output_filename);
            }
        }
    }

    // Recursively iterate through the subtree and count the number of substitutions.
    tuple<double, double> evolutionary_rates(string const &source, string const &target) {
        double dn{0}, ds{0};

        for (auto &substitution : substitutions) {
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

class SimuEvolArgParse : public SimuArgParse {
  public:
    explicit SimuEvolArgParse(CmdLine &cmd) : SimuArgParse(cmd) {}

    TCLAP::ValueArg<double> mu{"m", "mu", "Mutation rate", false, 2.5, "double", cmd};
    TCLAP::ValueArg<double> beta{
        "b", "beta", "Effective population size (relative)", false, 1.0, "double", cmd};
    TCLAP::ValueArg<double> selection{"s", "selection",
        "Selection coefficient given the current amino-acid", false, 0.0, "double", cmd};
    TCLAP::ValueArg<double> shuffle_proba{"p", "shuffle_proba",
        "Probability to randomize the fitness landscape (once a substitution occured)", false, 0.0,
        "double", cmd};
    SwitchArg shuffle_all{"a", "shuffle_all",
        "All sites are affected by the random shuffling (instead of just the one where the "
        "substitution occured)",
        cmd, false};
};

int main(int argc, char *argv[]) {
    CmdLine cmd{"SimuEvol", ' ', "0.1"};
    SimuEvolArgParse args(cmd);
    cmd.parse(argc, argv);
    Trace trace;

    string output_path{args.output_path.getValue()};
    double mu{args.mu.getValue()};
    double s{args.selection.getValue()};
    double p{args.shuffle_proba.getValue()};
    bool all_sites{args.shuffle_all.getValue()};

    vector<array<double, 20>> fitness_profiles =
        open_preferences(args.preferences_path.getValue(), args.beta.getValue());
    auto nbr_sites = static_cast<unsigned>(fitness_profiles.size());

    NucleotideRateMatrix mutation_rate(args.nuc_matrix_path.getValue(), mu, true);

    Sequence_dna root_sequence(mutation_rate, fitness_profiles, s, p, all_sites);

    std::ifstream tree_stream{args.newick_path.getValue()};
    NHXParser parser{tree_stream};
    std::unique_ptr<const Tree> tree;
    tree = make_from_parser(parser);

    init_alignments(output_path, tree->nb_leaves(), nbr_sites * 3);

    trace.add("Tree nodes", tree->nb_nodes());
    trace.add("Tree leaves", tree->nb_leaves());
    trace.add("Tree length", tree->total_length());
    trace.add("Number codons", nbr_sites);
    trace.add("Number nucleotides", nbr_sites * 3);
    trace.add("Predicted omega",
        root_sequence.predicted_omega(Codon::nucleotides, Codon::nucleotides, false));
    root_sequence.at_equilibrium();
    root_sequence.write_matrices(output_path);

    Process simu_process(*tree, root_sequence);
    simu_process.run(output_path);

    long nbr_non_synonymous, nbr_synonymous;
    tie(nbr_non_synonymous, nbr_synonymous) = simu_process.nbr_substitutions();

    // .txt output
    trace.add("Number substitutions", nbr_synonymous + nbr_non_synonymous);
    trace.add("Number substitutions per site",
        static_cast<double>(nbr_synonymous + nbr_non_synonymous) / nbr_sites);
    trace.add("Synonymous substitutions", nbr_synonymous);
    trace.add("Non-synonymous substitutions", nbr_non_synonymous);
    trace.add(
        "Simulated omega", simu_process.simulated_omega(Codon::nucleotides, Codon::nucleotides));

    trace.write_tsv(output_path);
    cout << "Simulation computed." << endl;
    cout << "Statistics summarized in: " << output_path + ".tsv" << endl;
    cout << "Fasta file in: " << output_path + ".fasta" << endl;
    cout << "Alignment (.ali) file in: " << output_path + ".ali" << endl;
    return 0;
}