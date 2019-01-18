#include <fstream>
#include <iostream>
#include <random>
#include "Eigen/Dense"
#include "argparse.hpp"
#include "codon.hpp"
#include "statistic.hpp"
#include "tree.hpp"

using namespace TCLAP;
using namespace std;

typedef Eigen::Matrix<double, 4, 4> Matrix4x4;
typedef Eigen::Matrix<double, 4, 1> Vector4x1;

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

// Compute the kernel of the mutation-rate matrix.
// This kernel is a vector of the nucleotides frequencies (not normalized to 1) at equilibrium.
Vector4x1 equilibrium_frequencies(Matrix4x4 const &mutation_matrix) {
    Eigen::Matrix<double, 4, Eigen::Dynamic> kernel =
        mutation_matrix.transpose().fullPivLu().kernel();
    Vector4x1 nuc_frequencies;

    if (kernel.cols() > 1) {
        cerr << "The kernel has " << kernel.cols() << " dimensions, this is weird ! " << endl;
        uniform_int_distribution<unsigned> unif_int(0, unsigned(kernel.cols()) - 1);
        unsigned chosen_row = unif_int(Codon::re);
        nuc_frequencies = kernel.col(chosen_row);

    } else {
        nuc_frequencies = kernel.col(0);
    }
    nuc_frequencies /= nuc_frequencies.sum();

    return nuc_frequencies;
}


Matrix4x4 normalize_submatrix(Matrix4x4 mutation_matrix) {
    for (int diag{0}; diag < 4; diag++) { mutation_matrix(diag, diag) = 0.0; }
    mutation_matrix -= mutation_matrix.rowwise().sum().asDiagonal();


    double events{0.};
    Vector4x1 nuc_frequencies = equilibrium_frequencies(mutation_matrix);
    for (int diag{0}; diag < 4; diag++) {
        events -= nuc_frequencies(diag) * mutation_matrix(diag, diag);
    }
    mutation_matrix /= events;
    return mutation_matrix;
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
    Matrix4x4 mutation_rate_matrix;

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
    explicit Sequence_dna(Matrix4x4 const &mutation_rate,
        vector<array<double, 20>> const &fitness_profiles, double const s,
        double const proba_permutation, bool const all_sites)
        : nbr_sites{static_cast<unsigned>(fitness_profiles.size())},
          codon_seq(nbr_sites, 0),
          aa_fitness_profiles{fitness_profiles},
          sel_coef{s},
          proba_permutation{proba_permutation},
          all_sites{all_sites} {
        mutation_rate_matrix = mutation_rate;
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
                    // Selective strength between the mutated and original amino-acids.
                    double s{0.};
                    s = aa_fitness_profiles[site][Codon::codon_to_aa_array[codon_to]];
                    s -= aa_fitness_profiles[site][Codon::codon_to_aa_array[codon_from]];
                    // If the selective strength is 0, the rate of fixation is neutral.
                    // Else, the rate of fixation is computed using population genetic formulas
                    // (Kimura).
                    if (fabs(s) <= Codon::epsilon) {
                        rate_substitution = mutation_rate_matrix(n_from, n_to);
                    } else {
                        // The substitution rate is the mutation rate multiplied by the rate of
                        // fixation.
                        rate_substitution = mutation_rate_matrix(n_from, n_to) * s / (1 - exp(-s));
                    }
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

            unsigned index = substitution_distr(Codon::re);
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

            if (non_syn_opp_matrix.sum() > 10e10 or syn_opp_matrix.sum() > 10e10) {
                printf("Shit !!!!");
            }

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
                double rand_uni = unif_rand_proba(Codon::re);
                if (rand_uni < proba_permutation) {
                    if (all_sites) {
                        for (unsigned site_shuffle{0}; site_shuffle < nbr_sites; site_shuffle++) {
                            shuffle(aa_fitness_profiles[site_shuffle].begin(),
                                aa_fitness_profiles[site_shuffle].end(), Codon::re);
                        }
                    } else {
                        shuffle(aa_fitness_profiles[site].begin(), aa_fitness_profiles[site].end(),
                            Codon::re);
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
        Vector4x1 nuc_frequencies = equilibrium_frequencies(mutation_rate_matrix);

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
            char chosen_codon = freq_codon_distr(Codon::re);
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
                                // Selective strength between the mutated and original amino-acids.
                                double s{0.};
                                // Rate of fixation initialized to 1 (neutral mutation)
                                double rate_fixation{1.};

                                s = aa_fitness_profiles[site][Codon::codon_to_aa_array[codon_to]];
                                s -=
                                    aa_fitness_profiles[site][Codon::codon_to_aa_array[codon_from]];
                                // If the selective strength is 0, the rate of fixation is neutral.
                                // Else, the rate of fixation is computed using population genetic
                                // formulas (Kimura).
                                if (fabs(s) <= Codon::epsilon) {
                                    rate_fixation = 1;
                                } else {
                                    rate_fixation = s / (1 - exp(-s));
                                }

                                dn += codon_freqs[codon_from] * mutation_rate_matrix(n_from, n_to) *
                                      rate_fixation;
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

    void run(string output_filename) { run_recursive(tree.root(), output_filename); }

    // Recursively iterate through the subtree.
    void run_recursive(Tree::NodeIndex node, string output_filename) {
        // Substitutions of the DNA sequence is generated.
        sequences[node]->run_substitutions(tree.node_length(node), substitutions);

        if (tree.is_leaf(node)) {
            // If the node is a leaf, output the DNA sequence and name.
            string dna_str = sequences[node]->get_dna_str();

            // .ali format
            ofstream ali_file;
            ali_file.open(output_filename + ".ali", ios_base::app);
            ali_file << tree.node_name(node) << " " << dna_str << endl;
            ali_file.close();

            // .fasta format
            ofstream fasta_file;
            fasta_file.open(output_filename + ".fasta", ios_base::app);
            fasta_file << ">" << tree.node_name(node) << endl << dna_str << endl;
            fasta_file.close();
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
    TCLAP::ValueArg<double> lambda{
        "l", "lambda", "Strong to weak mutation bias", false, 5.0, "double", cmd};
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

    string preferences_path{args.preferences_path.getValue()};
    string newick_path{args.newick_path.getValue()};
    string output_path{args.output_path.getValue()};

    double mu{args.mu.getValue()};
    double lambda{args.lambda.getValue()};
    double beta{args.beta.getValue()};
    double s{args.selection.getValue()};
    double p{args.shuffle_proba.getValue()};
    bool all_sites{args.shuffle_all.getValue()};

    vector<array<double, 20>> fitness_profiles = open_preferences(preferences_path, beta);
    auto nbr_sites = static_cast<unsigned>(fitness_profiles.size());

    Matrix4x4 mutation_rate;
    mutation_rate << 0, 1, 1, lambda, lambda, 0, 1, lambda, lambda, 1, 0, lambda, lambda, 1, 1, 0;
    mutation_rate = normalize_submatrix(mutation_rate);
    mutation_rate *= mu;

    Sequence_dna root_sequence(mutation_rate, fitness_profiles, s, p, all_sites);

    std::ifstream tree_stream{newick_path};
    NHXParser parser{tree_stream};
    std::unique_ptr<const Tree> tree;
    tree = make_from_parser(parser);

    ofstream ali_file;
    ali_file.open(output_path + ".ali");
    ali_file << tree->nb_leaves() << " " << nbr_sites * 3 << endl;
    ali_file.close();

    ofstream fasta_file;
    fasta_file.open(output_path + ".fasta");
    fasta_file.close();

    // .txt output
    ofstream txt_file;
    txt_file.open(output_path + ".txt");
    txt_file << "The tree contains " << tree->nb_nodes() << " nodes for " << tree->nb_leaves()
             << " species at the tips." << endl;
    txt_file << "The tree has a total branch length of " << tree->total_length() << "." << endl;
    txt_file << "The DNA sequence is " << nbr_sites * 3 << " base pairs long." << endl;
    txt_file << "The mutation transition matrix (" << Codon::nucleotides << ") is: " << endl;
    txt_file << mutation_rate << endl;
    txt_file << "s=" << s << endl;
    txt_file << "all_sites=" << all_sites << endl;
    txt_file << "p=" << p << endl;
    txt_file << "w0="
             << root_sequence.predicted_omega(Codon::nucleotides, Codon::nucleotides, false)
             << endl;
    txt_file << "<w0>="
             << root_sequence.predicted_omega(Codon::nucleotides, Codon::nucleotides, true) << endl;
    txt_file.close();

    root_sequence.at_equilibrium();
    Process simu_process(*tree, root_sequence);
    simu_process.run(output_path);

    long nbr_non_synonymous, nbr_synonymous;
    tie(nbr_non_synonymous, nbr_synonymous) = simu_process.nbr_substitutions();

    // .txt output
    txt_file.open(output_path + ".txt", ios_base::app);
    txt_file << "The simulation mapped " << nbr_synonymous + nbr_non_synonymous
             << " substitutions along the tree." << endl;
    txt_file << "On average, this is "
             << static_cast<double>(nbr_synonymous + nbr_non_synonymous) / nbr_sites
             << " substitutions per site." << endl;
    txt_file << nbr_synonymous << " synonymous and " << nbr_non_synonymous
             << " non-synonymous substitutions." << endl;
    txt_file << "w=" << simu_process.simulated_omega(Codon::nucleotides, Codon::nucleotides)
             << endl;

    /*
        string weak_strong = "WS";
        map<char, string> const WS_map{{'W', "AT"},
                                       {'S', "GC"}};
        for (auto subset_from : weak_strong) {
            for (auto subset_to : weak_strong) {
                txt_file << "w0_" << subset_from << subset_to << "="
                         << root.predicted_omega(WS_map.at(subset_from), WS_map.at(subset_to),
       false) << endl; txt_file << "<w0>_" << subset_from << subset_to << "="
                         << root.predicted_omega(WS_map.at(subset_from), WS_map.at(subset_to), true)
       << endl; txt_file << "w_" << subset_from << subset_to << "="
                         << root.simulated_omega(WS_map.at(subset_from), WS_map.at(subset_to)) <<
       endl;
            }
        }

        for (auto nuc_from : Codon::nucleotides) {
            for (auto nuc_to : Codon::nucleotides) {
                if (nuc_from != nuc_to) {
                    txt_file << "w0_" << nuc_from << nuc_to << "="
                             << root.predicted_omega(char_to_str(nuc_from), char_to_str(nuc_to),
       false) << endl; txt_file << "<w0>_" << nuc_from << nuc_to << "="
                             << root.predicted_omega(char_to_str(nuc_from), char_to_str(nuc_to),
       true) << endl; txt_file << "w_" << nuc_from << nuc_to << "="
                             << root.simulated_omega(char_to_str(nuc_from), char_to_str(nuc_to)) <<
       endl;
                }
            }
        }
    */

    txt_file.close();

    cout << "Simulation computed. Log of the simulation available at: " << endl;
    cout << output_path + ".txt" << endl;
    return 0;
}