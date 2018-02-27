#include <iostream>
#include <random>
#include <fstream>
#include "codon.hpp"

using namespace std;

#include "Eigen/Dense"

typedef Eigen::Matrix<double, 4, 4> Matrix4x4;

#define DOCOPT_HEADER_ONLY

#include "docopt.cpp/docopt.h"

struct Substitution {
    unsigned site{0};
    char codon_from{0};
    char codon_to{0};
    double dn{0.};
    double ds{0.};
};

//Function for average
double sum(vector<double> &v) {
    double return_value = 0.0;
    size_t n = v.size();

    for (int i = 0; i < n; i++) {
        return_value += v[i];
    }

    return return_value;
}

//Function for average
double avg(vector<double> &v) {
    double return_value = 0.0;
    size_t n = v.size();

    for (int i = 0; i < n; i++) {
        return_value += v[i];
    }

    return return_value / n;
}

//Function for variance
double variance(vector<double> &v) {
    double v_squarred = 0.0;
    size_t n = v.size();

    for (int i = 0; i < n; i++) {
        v_squarred += pow(v[i], 2);
    }

    return (v_squarred / n) - pow(avg(v), 2);
}

bool is_synonymous(Substitution const &s) {
    return (Codon::codon_to_aa_array[s.codon_from] == Codon::codon_to_aa_array[s.codon_to]);
}


// Class representing DNA sequence.
class Sequence_dna {
private:
    // The number of sites in the sequence (each position is a codon, thus the DNA sequence is 3 times greater).
    unsigned const nbr_sites;

    // The sequence of codons.
    vector<char> codon_seq;

    // The matrix of mutation rates between nucleotides.
    Matrix4x4 mutation_rate_matrix;

    // The fitness profiles of amino-acids.
    vector<array<double, 20>> aa_fitness_profiles;

    // The probability to randomize the fitness landscape.
    double proba_permutation;

public:
    // Constructor of Sequence_dna.
    // size: the size of the DNA sequence.
    explicit Sequence_dna(unsigned const size) : nbr_sites{size}, codon_seq(size, 0),
                                                 aa_fitness_profiles{0}, proba_permutation{0.} {
        mutation_rate_matrix.Zero();
    }

    // Method translating the codon sequence to amino-acid sequence.
    string translate() const {
        string protein(nbr_sites, ' ');

        // For each site
        for (unsigned site{0}; site < nbr_sites; site++) {
            // Use the codon_to_aa_array to translate site to amino-acid.
            protein[site] = Codon::amino_acids[Codon::codon_to_aa_array[codon_seq[site]]];
        }

        return protein; // return the amino-acid sequence as a string.
    }

    // Method computing the next substitution event to occur, and the time for it to happen.
    // This method takes a time as input. If there is enough time given, the substitution event is computed.
    // This method also returns the time (given as input) decremented by the time needed for the substitution event to occur.
    // If there is not enough time given, no substitution event is computed and the method returns 0.
    // time_left: The time available for a substitution event to occur.
    double next_substitution(double &time_left, vector<Substitution> &substitutions_vec) {

        // Number of possible substitutions is 9 times the number of sites (3 substitutions for each 3 possible positions).
        unsigned nbr_substitutions{9 * nbr_sites};

        // Vector of substitution rates.
        vector<double> substitution_rates(nbr_substitutions, 0);

        // Sum of substitution rates.
        double total_substitution_rates{0.};

        double dn{.0}, ds{.0};

        // For all site of the sequence.
        for (unsigned site{0}; site < nbr_sites; site++) {
            // Codon original before substitution.

            char codon_from = codon_seq[site];

            // Array of neighbors of the original codon (codons differing by only 1 mutation).
            array<tuple<char, char, char>, 9> neighbors = Codon::codon_to_neighbors_array[codon_from];

            // For all possible neighbors.
            for (char neighbor{0}; neighbor < 9; neighbor++) {

                // Codon after mutation, Nucleotide original and Nucleotide after mutation.
                char codon_to{0}, n_from{0}, n_to{0};
                tie(codon_to, n_from, n_to) = neighbors[neighbor];

                // Assign the substitution rate given by the method substitution rate.
                // Rate of substitution initialized to 0 (deleterious mutation)
                double rate_substitution{0.};

                // If the mutated amino-acid is a stop codon, the rate of fixation is 0.
                // Else, if the mutated and original amino-acids are non-synonymous, we compute the rate of fixation.
                // Note that, if the mutated and original amino-acids are synonymous, the rate of fixation is 1.
                if (Codon::codon_to_aa_array[codon_to] == 20) {
                    rate_substitution = 0.;
                } else if (Codon::codon_to_aa_array[codon_from] != Codon::codon_to_aa_array[codon_to]) {
                    // Selective strength between the mutated and original amino-acids.
                    double s{0.};
                    s = aa_fitness_profiles[site][Codon::codon_to_aa_array[codon_to]];
                    s -= aa_fitness_profiles[site][Codon::codon_to_aa_array[codon_from]];
                    // If the selective strength is 0, the rate of fixation is neutral.
                    // Else, the rate of fixation is computed using population genetic formulas (Kimura).
                    if (fabs(s) <= Codon::epsilon) {
                        rate_substitution = mutation_rate_matrix(n_from, n_to);
                    } else {
                        // The substitution rate is the mutation rate multiplied by the rate of fixation.
                        rate_substitution = mutation_rate_matrix(n_from, n_to) * s / (1 - exp(-s));
                    }
                    dn += mutation_rate_matrix(n_from, n_to);
                } else {
                    rate_substitution = mutation_rate_matrix(n_from, n_to);
                    ds += mutation_rate_matrix(n_from, n_to);
                }

                substitution_rates[9 * site + neighbor] = rate_substitution;
                // Increment the sum of substitution rates
                total_substitution_rates += rate_substitution;

            }
        }

        // Decrement the time by inverse of the sum of substitution rates.
        time_left -= 1. / total_substitution_rates;

        // Substitute the sequence if the time is positive, else there is no substitution but the time left is set to 0.
        if (time_left > 0. and total_substitution_rates != 0.) {

            // Uniform random generator between 0 and the total sum of substitution rates.
            uniform_real_distribution<double> unif_rand_total_substitution_rates(0, total_substitution_rates);
            double random_cumulative_substitution_rates = unif_rand_total_substitution_rates(Codon::re);
            double cumulative_substitution_rates{0.};

            unsigned index = 0;
            for (unsigned t{0}; t < substitution_rates.size(); t++) {
                // Iterate through the cumulative substitution rates and break the loop when it is greater than the random cumulative substitution rates
                cumulative_substitution_rates += substitution_rates[t];
                if (random_cumulative_substitution_rates < cumulative_substitution_rates) {
                    index = t;
                    break;
                }
            }
            unsigned site = index / 9;
            char codom_from = codon_seq[site];

            // Array of neighbors of the original codon (codons differing by only 1 mutation).
            array<tuple<char, char, char>, 9> neighbors = Codon::codon_to_neighbors_array[codom_from];

            // Codon after mutation, Nucleotide original and Nucleotide after mutation.
            char codon_to{0}, n_from{0}, n_to{0};

            tie(codon_to, n_from, n_to) = neighbors[index % 9];

            Substitution sub = {site, codom_from, codon_to, 1. / dn, 1. / ds};
            substitutions_vec.push_back(sub);

            // Random shuffle of the fitness landscape
            uniform_real_distribution<double> unif_rand_proba(0, 1);
            double rand_uni = unif_rand_proba(Codon::re);
            if (rand_uni < proba_permutation) {
                shuffle(aa_fitness_profiles[site].begin(), aa_fitness_profiles[site].end(), Codon::re);
            }

            codon_seq[site] = codon_to;

        } else if (time_left < 0.) {
            time_left = 0.;
        }
        return time_left;
    }

    // Method computing all substitution event occuring during a given time-frame.
    // t: time during which substitution events occur (typically branch length).
    // This method is o(n²) where n is the number of sites, but can take into account epistatic effects
    void run_substitutions(double t, vector<Substitution> &subs) {
        while (t > 0) {
            t = next_substitution(t, subs);
        }
    }


    // Method computing the equilibrium frequencies for one site.
    // aa_fitness_profil: The amino-acid fitness profil of the given site.
    array<double, 64> codon_frequencies(array<double, 20> const &aa_fitness_profil) {
        // Compute the kernel of the mutation-rate matrix.
        // This kernel is a vector of the nucleotides frequencies (not normalized to 1) at equilibrium.
        Eigen::Matrix<double, 4, Eigen::Dynamic> kernel = mutation_rate_matrix.transpose().fullPivLu().kernel();
        Eigen::Matrix<double, 4, 1> nuc_frequencies;

        if (kernel.cols() > 1) {
            cerr << "The kernel has " << kernel.cols() << " dimensions, this is weird ! " << endl;
            uniform_int_distribution<unsigned> unif_int(0, unsigned(kernel.cols()) - 1);
            unsigned chosen_row = unif_int(Codon::re);
            nuc_frequencies = kernel.col(chosen_row);

        } else {
            nuc_frequencies = kernel.col(0);
        }

        array<double, 64> codon_frequencies{0};

        // Initialize the total sum of equilibrium frequency at 0.
        double total_frequencies{0.};

        // For each site of the vector of the site frequencies.
        for (char codon{0}; codon < 64; codon++) {
            double codon_freq{1.};

            // For all nucleotides in the codon
            for (auto &nuc: Codon::codon_to_triplet_array[codon]) {
                codon_freq *= nuc_frequencies(nuc);
            }

            if (Codon::codon_to_aa_array[codon] != 20) {
                codon_frequencies[codon] = codon_freq * exp(aa_fitness_profil[Codon::codon_to_aa_array[codon]]);
            } else {
                codon_frequencies[codon] = 0.;
            }


            // Increment the total sum of equilibrium frequencies.
            total_frequencies += codon_frequencies[codon];
        }

        // Normalize the vector of equilibrium frequencies.
        for (char codon{0}; codon < 64; codon++) {
            codon_frequencies[codon] /= total_frequencies;
        }

        return codon_frequencies;
    };

    // Set the the DNA sequence to the mutation-selection equilibrium.
    void at_equilibrium() {

        // For all site of the sequence.
        for (unsigned site{0}; site < nbr_sites; site++) {

            array<double, 64> codon_freqs = codon_frequencies(aa_fitness_profiles[site]);

            // Uniform random generator between 0 and the total sum of equilibrium frequencies.
            uniform_real_distribution<double> unif(0., 1.);
            double random_cumulative_frequencies = unif(Codon::re);
            double cumulative_frequencies{0.};

            char index{0};
            for (char m{0}; m < 64; m++) {
                // Iterate through the cumulative frequencies and break the loop when it is greater than the random cumulative frequencies.
                cumulative_frequencies += codon_freqs[m];
                if (random_cumulative_frequencies < cumulative_frequencies) {
                    index = m;
                    break;
                }
            }

            // Substitute the site with the substitution given by the loop break.
            codon_seq[site] = index;
        }
    }

    // Set attribute method for the codon sequence.
    void set_parameters(Sequence_dna const &sequence_dna) {
        codon_seq = sequence_dna.codon_seq;
        mutation_rate_matrix = sequence_dna.mutation_rate_matrix;
        aa_fitness_profiles = sequence_dna.aa_fitness_profiles;
        proba_permutation = sequence_dna.proba_permutation;
    }

    // Set attribute method for the mutation rate matrix.
    void set_mutation_rate(Matrix4x4 const &mutation_rate) {
        mutation_rate_matrix = mutation_rate;
    }

    // Set attribute method for the amino-acid fitness profil.
    void set_fitness_profiles(vector<array<double, 20>> const &fitness_profiles) {
        assert(nbr_sites == fitness_profiles.size());
        aa_fitness_profiles = fitness_profiles;
    }

    // Set attribute method for the mutation rate matrix.
    void set_proba_permutation(double const proba_randomize) {
        proba_permutation = proba_randomize;
    }

    // Theoretical computation of the predicted omega
    double predicted_omega(bool averaged=true) {
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
                    for (auto &neighbor: Codon::codon_to_neighbors_array[codon_from]) {

                        // Codon after mutation, Nucleotide original and Nucleotide after mutation.
                        char codon_to{0}, n_from{0}, n_to{0};
                        tie(codon_to, n_from, n_to) = neighbor;


                        // If the mutated amino-acid is a stop codon, the rate of fixation is 0.
                        // Else, if the mutated and original amino-acids are non-synonymous, we compute the rate of fixation.
                        // Note that, if the mutated and original amino-acids are synonymous, the rate of fixation is 1.
                        if (Codon::codon_to_aa_array[codon_to] != 20 and
                            Codon::codon_to_aa_array[codon_from] != Codon::codon_to_aa_array[codon_to]) {
                            // Selective strength between the mutated and original amino-acids.
                            double s{0.};
                            // Rate of fixation initialized to 1 (neutral mutation)
                            double rate_fixation{1.};

                            s = aa_fitness_profiles[site][Codon::codon_to_aa_array[codon_to]];
                            s -= aa_fitness_profiles[site][Codon::codon_to_aa_array[codon_from]];
                            // If the selective strength is 0, the rate of fixation is neutral.
                            // Else, the rate of fixation is computed using population genetic formulas (Kimura).
                            if (fabs(s) <= Codon::epsilon) {
                                rate_fixation = 1;
                            } else {
                                rate_fixation = s / (1 - exp(-s));
                            }

                            dn += codon_freqs[codon_from] * mutation_rate_matrix(n_from, n_to) * rate_fixation;
                            d0 += codon_freqs[codon_from] * mutation_rate_matrix(n_from, n_to);
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
        return dna_str; // return the DNA sequence as a string.
    }
};


// Class representing nodes of a tree.
class Node {
private:
    string name;  // The species name of the node.
    double length;  // The length of the branch attached to the node (ascending).
    string newick;  // The newick tree descending from the node.
    vector<Node> children;  // Vector of direct children (first order, no grand-children).
    Sequence_dna sequence_dna;  // The DNA sequence attached to the node.
    vector<Substitution> substitutions;  // The number of substitutions in the branch attached to the node.

public:

    // Constructor
    Node(string name, string const &len, string newick, Sequence_dna seq) :
            name{move(name)}, length{stod(len)}, newick{move(newick)}, sequence_dna{move(seq)} {

        // Parse the newick tree descending of the node.
        parse_newick();
    }

    Node(string newick, unsigned const nbr_sites) :
            name{"Root"}, length{0.}, newick{move(newick)}, sequence_dna(nbr_sites) {

        // Parse the newick tree descending of the node.
        parse_newick();
    }

    // Is true if the node don't have children.
    bool is_leaf() const {
        return newick.length() == 0;
    }

    // Add a node as the vector of children.
    void add_child(Node const &node) {
        children.push_back(node);
    }

    // Recursively iterate through the subtree and count the number of nodes.
    double tot_length() {
        double tot_length = length;

        if (!is_leaf()) {
            // Else, if the node is internal, iterate through the direct children.
            for (auto &child : children) {
                tot_length += child.tot_length();
            }
        }
        return tot_length;
    }

    // Recursively iterate through the subtree and count the number of nodes.
    unsigned nbr_nodes() {
        unsigned nbr_nodes = 1;

        if (!is_leaf()) {
            // Else, if the node is internal, iterate through the direct children.
            for (auto &child : children) {
                nbr_nodes += child.nbr_nodes();
            }
        }
        return nbr_nodes;
    }

    // Recursively iterate through the subtree and count the number of leaves.
    unsigned nbr_leaves() {
        unsigned nbr_leaves = 0;

        if (is_leaf()) {
            // If the node is a leaf, return 1.
            nbr_leaves = 1;
        } else {
            // Else, if the node is internal, iterate through the direct children.
            for (auto &child : children) {
                nbr_leaves += child.nbr_leaves();
            }
        }
        return nbr_leaves;
    }

    // Recursively iterate through the subtree and count the number of substitutions.
    tuple<long, long> nbr_substitutions() {
        long nbr_non_syn{0}, nbr_syn{0};

        nbr_syn = count_if(substitutions.begin(), substitutions.end(), is_synonymous);
        nbr_non_syn = substitutions.size() - nbr_syn;

        if (!is_leaf()) {
            // Else, if the node is internal, iterate through the direct children.
            for (auto &child : children) {
                long child_nbr_syn{0}, child_nbr_non_syn{0};
                tie(child_nbr_non_syn, child_nbr_syn) = child.nbr_substitutions();
                nbr_non_syn += child_nbr_non_syn;
                nbr_syn += child_nbr_syn;
            }
        }
        return make_tuple(nbr_non_syn, nbr_syn);
    }

    // Recursively iterate through the subtree.
    void traverse(string &output_filename) {
        // Substitutions of the DNA sequence is generated.
        sequence_dna.run_substitutions(length, substitutions);

        if (is_leaf()) {
            // If the node is a leaf, output the DNA sequence and name.
            string dna_str = sequence_dna.get_dna_str();

            // .ali format
            ofstream ali_file;
            ali_file.open(output_filename + ".ali", ios_base::app);
            ali_file << name << " " << dna_str << endl;
            ali_file.close();

            // .fasta format
            ofstream fasta_file;
            fasta_file.open(output_filename + ".fasta", ios_base::app);
            fasta_file << ">" << name << endl << dna_str << endl;
            fasta_file.close();
        } else {
            // If the node is internal, iterate through the direct children.
            for (auto &child : children) {
                child.sequence_dna.set_parameters(sequence_dna);
                child.traverse(output_filename);
            }
        }
    }

    // Set method for the parameters of evolution of the sequence
    void set_evolution_parameters(Matrix4x4 const &mutation_rate,
                                  vector<array<double, 20>> const &fitness_profiles,
                                  double const proba_permutation) {
        sequence_dna.set_mutation_rate(mutation_rate);
        sequence_dna.set_fitness_profiles(fitness_profiles);
        sequence_dna.set_proba_permutation(proba_permutation);
    }

    // Theoretical prediction of omega from the parameters of the simulation
    double predicted_omega(bool average) {
        return sequence_dna.predicted_omega(average);
    }

    // Recursively iterate through the subtree and count the number of substitutions.
    tuple<double, double> evolutionary_rates() {
        double dn{0}, ds{0};

        for (auto &substitution: substitutions) {
            if (is_synonymous(substitution)){
                ds += substitution.ds;
            } else {
                dn += substitution.dn;
            }
        }

        if (!is_leaf()) {
            // Else, if the node is internal, iterate through the direct children.
            for (auto &child : children) {
                double dn_child{0.}, ds_child{0.};
                tie(dn_child, ds_child) = child.evolutionary_rates();
                dn += dn_child;
                ds += ds_child;
            }
        }
        return make_tuple(dn, ds);
    }

    // Simulated omega from the substitutions
    double simulated_omega() {
        double dn{0}, ds{0};
        tie(dn, ds) = evolutionary_rates();
        if (ds == .0) {
            cerr << "There is no synonymous substitutions generated by the simulation, dN/dS can't be computed!" << endl;
            return .0;
        } else {
            return dn / ds;
        }
    }

    // Set the the DNA sequence to the mutation-selection equilibrium.
    void at_equilibrium() {
        sequence_dna.at_equilibrium();
    }

    void parse_newick() {

        if (!is_leaf()) {
            // The size of the string of newick tree.
            size_t max_position{newick.size()};
            // The current position in the string of the newick tree.
            size_t position{0};

            // While the current position is lower than the size of the string, their is at least one node to parse.
            while (position < max_position) {
                // 'subtree' is the left hand side of the node name, it can be a subtree or nothing if the node is a leaf.
                string subtree{};
                if (newick[position] == '(') {

                    size_t postpoint{position};
                    unsigned nbr_open{1};

                    for (size_t i{position + 1}; i < max_position; i++) {
                        if (nbr_open == 0) {
                            postpoint = i;
                            break;
                        } else if (newick[i] == '(') {
                            nbr_open++;
                        } else if (newick[i] == ')') {
                            nbr_open--;
                        };
                    }
                    subtree = newick.substr(position + 1, postpoint - position - 2);
                    position = postpoint;
                }

                // 'name_suffix' contains the name of the node and the branch length.
                string name_suffix{};

                size_t next_sep = newick.substr(position).find(',');
                if (next_sep == string::npos) {
                    name_suffix = newick.substr(position);
                    position = max_position;
                } else {
                    name_suffix = newick.substr(position, next_sep);
                    position = position + next_sep + 1;
                }

                // 'length' contains the name of the node.
                string length{};
                // 'name' contains the branch length of the node.
                string name{};
                size_t ddot = name_suffix.rfind(':');
                if (ddot != string::npos) {
                    length = name_suffix.substr(ddot + 1);
                    name = name_suffix.substr(0, ddot);
                } else {
                    name = name_suffix;
                    length = "0";
                }

                // New node from 'subtree', 'name' and 'length' using the DNA sequence of this node.
                add_child(Node(name, length, subtree, sequence_dna));
            }
        }
    }
};


string open_newick(string const &file_name) {
    ifstream input_stream(file_name);
    if (!input_stream) cerr << "Can't open newick file!" << endl;

    string line;
    getline(input_stream, line);

    return line;
}

vector<array<double, 20>> open_preferences(string const &file_name) {
    vector<array<double, 20>> fitness_profiles{0};

    ifstream input_stream(file_name);
    if (!input_stream) cerr << "Can't open preferences file!" << endl;

    string line;

    // skip the header of the file
    getline(input_stream, line);

    while (getline(input_stream, line)) {

        array<double, 20> fitness_profil{0};
        string word;
        istringstream line_stream(line);
        unsigned counter{0};

        while (getline(line_stream, word, ' ')) {
            if (counter > 2) {
                fitness_profil[counter - 3] = log(stod(word));
            }
            counter++;
        }

        fitness_profiles.push_back(fitness_profil);
    }
    return fitness_profiles;
}

static char const USAGE[] =
        R"(
Usage:
      SimuEvol [--preferences=<file_path>] [--newick=<file_path>] [--output=<file_path>] [--mu=<0.5>] [--lambda=<3>] [--p=<0.0>]
      SimuEvol --help
      SimuEvol --version

Options:
-h --help                    show this help message and exit
--version                    show version and exit
--preferences=<file_path>    specify input site-specific preferences file [default: ../data/gal4.txt]
--newick=<file_path>         specify input newick tree [default: ../data/gal4.newick]
--output=<file_path>         specify output protein name [default: ../data/gal4]
--mu=<0.5>                   specify the mutation rate [default: 0.5]
--lambda=<3>                 specify the strong to weak mutation bias [default: 3]
--p=<0.0>                    specify the probability to randomize the fitness landscape [default: 0.0]
)";

int main(int argc, char *argv[]) {

    auto args = docopt::docopt(USAGE,
                               {argv + 1, argv + argc},
                               true,               // show help if requested
                               "SimuEvol 0.1a");  // version string

    string preferences_path{"../data_prefs/np.txt"};
    if (args["--preferences"]) {
        preferences_path = args["--preferences"].asString();
    }

    string newick_path{"../data_trees/np.newick"};
    if (args["--newick"]) {
        newick_path = args["--newick"].asString();
    }

    string output_path{"../data_alignment/np"};
    if (args["--output"]) {
        output_path = args["--output"].asString();
    }

    vector<array<double, 20>> fitness_profiles = open_preferences(preferences_path);
    auto nbr_sites = static_cast<unsigned>(fitness_profiles.size());

    string newick_tree = open_newick(newick_path);
    Node root(newick_tree, nbr_sites);

    Matrix4x4 mutation_rate;
    double mu = 0.5;
    if (args["--mu"]) {
        mu = stod(args["--mu"].asString());
    }
    double lambda = 10.0;
    if (args["--lambda"]) {
        lambda = stod(args["--lambda"].asString());
    }
    mutation_rate << 0, 1, 1, lambda,
            lambda, 0, 1, lambda,
            lambda, 1, 0, lambda,
            lambda, 1, 1, 0;
    mutation_rate *= mu;
    mutation_rate -= mutation_rate.rowwise().sum().asDiagonal();

    double p = 0.0;
    if (args["--p"]) {
        p = stod(args["--p"].asString());
    }

    ofstream ali_file;
    ali_file.open(output_path + ".ali");
    ali_file << root.nbr_leaves() << " " << nbr_sites * 3 << endl;
    ali_file.close();

    ofstream fasta_file;
    fasta_file.open(output_path + ".fasta");
    fasta_file.close();

    // .txt output
    ofstream txt_file;
    txt_file.open(output_path + ".txt");
    txt_file << "The tree contains " << root.nbr_nodes() << " nodes for " << root.nbr_leaves() << " species at the tips."
         << endl;
    txt_file << "The tree has a total branch length of " << root.tot_length() << "." << endl;
    txt_file << "The DNA sequence is " << nbr_sites * 3 << " base pairs long." << endl;
    txt_file << "The mutation transition matrix (" << Codon::nucleotides << ") is: " << endl;
    txt_file << mutation_rate << endl;
    txt_file << "There is a probability of " << p << " to permute the fitness landscape if a substitution occurs." << endl;
    txt_file.close();

    root.set_evolution_parameters(mutation_rate, fitness_profiles, p);
    root.at_equilibrium();
    root.traverse(output_path);

    long nbr_non_synonymous, nbr_synonymous;
    tie(nbr_non_synonymous, nbr_synonymous) = root.nbr_substitutions();

    // .txt output
    txt_file.open(output_path + ".txt", ios_base::app);
    txt_file << "The simulation mapped " << nbr_synonymous + nbr_non_synonymous << " substitutions along the tree." << endl;
    txt_file << "On average, this is " << static_cast<double>(nbr_synonymous + nbr_non_synonymous) / nbr_sites << " substitutions per site." << endl;
    txt_file << nbr_synonymous << " synonymous and " << nbr_non_synonymous << " non-synonymous substitutions." << endl;
    txt_file << "w0_avg=" << root.predicted_omega(true) << endl;
    txt_file << "w0=" << root.predicted_omega(false) << endl;
    txt_file << "w=" << root.simulated_omega() << endl;
    txt_file.close();

    cout << "Simulation computed. Log of the simulation available at: " << endl;
    cout << output_path + ".txt" << endl;
    return 0;
}