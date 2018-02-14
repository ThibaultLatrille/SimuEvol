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
};


bool is_synonymous(Substitution const &s) {
    return (Codon::codon_to_aa_array.at(s.codon_from) == Codon::codon_to_aa_array.at(s.codon_to));
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
            protein.at(site) = Codon::amino_acids.at(Codon::codon_to_aa_array.at(codon_seq.at(site)));
        }

        return protein; // return the amino-acid sequence as a string.
    }

    // Method computing the next substitution event to occur, and the time for it to happen.
    // This method takes a time as input. If there is enough time given, the substitution event is computed.
    // This method also returns the time (given as input) decremented by the time needed for the substitution event to occur.
    // If there is not enough time given, no substitution event is computed and the method returns 0.
    // time_left: The time available for a substitution event to occur.
    double next_substitution(double &time_left, vector<Substitution> &substitutions) {

        // Number of possible substitutions is 9 times the number of sites (3 substitutions for each 3 possible positions).
        unsigned nbr_substitutions{9 * nbr_sites};

        // Vector of substitution rates.
        vector<double> substitution_rates(nbr_substitutions, 0);

        // Sum of substitution rates.
        double total_substitution_rates{0.};

        // For all site of the sequence.
        for (unsigned site{0}; site < nbr_sites; site++) {
            // Codon original before substitution.

            char codon_from = codon_seq.at(site);

            // Array of neighbors of the original codon (codons differing by only 1 mutation).
            array<tuple<char, char, char>, 9> neighbors = Codon::codon_to_neighbors_array.at(codon_from);

            // For all possible neighbors.
            for (char neighbor{0}; neighbor < 9; neighbor++) {

                // Codon after mutation, Nucleotide original and Nucleotide after mutation.
                char codon_to{0}, n_from{0}, n_to{0};
                tie(codon_to, n_from, n_to) = neighbors.at(neighbor);

                // Assign the substitution rate given by the method substitution rate.
                // Rate of fixation initialized to 1 (neutral mutation)
                double rate_fixation{1.};

                // If the mutated amino-acid is a stop codon, the rate of fixation is 0.
                // Else, if the mutated and original amino-acids are non-synonymous, we compute the rate of fixation.
                // Note that, if the mutated and original amino-acids are synonymous, the rate of fixation is 1.
                if (Codon::codon_to_aa_array.at(codon_to) == 20) {
                    rate_fixation = 0.;
                } else if (Codon::codon_to_aa_array.at(codon_from) != Codon::codon_to_aa_array.at(codon_to)) {
                    // Selective strength between the mutated and original amino-acids.
                    double s{0.};
                    s = aa_fitness_profiles.at(site).at(Codon::codon_to_aa_array.at(codon_to));
                    s -= aa_fitness_profiles.at(site).at(Codon::codon_to_aa_array.at(codon_from));
                    // If the selective strength is 0, the rate of fixation is neutral.
                    // Else, the rate of fixation is computed using population genetic formulas (Kimura).
                    if (fabs(s) <= Codon::epsilon) {
                        rate_fixation = 1;
                    } else {
                        rate_fixation = s / (1 - exp(-s));
                    }
                }

                // The substitution rate is the mutation rate multiplied by the rate of fixation.
                substitution_rates.at(9 * site + neighbor) = rate_fixation * mutation_rate_matrix(n_from, n_to);
                // Increment the sum of substitution rates
                total_substitution_rates += substitution_rates.at(9 * site + neighbor);

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
            for (unsigned t{0}; t < nbr_substitutions; t++) {
                // Iterate through the cumulative substitution rates and break the loop when it is greater than the random cumulative substitution rates
                cumulative_substitution_rates += substitution_rates.at(t);
                if (random_cumulative_substitution_rates < cumulative_substitution_rates) {
                    index = t;
                    break;
                }
            }

            // Array of neighbors of the original codon (codons differing by only 1 mutation).
            array<tuple<char, char, char>, 9> neighbors = Codon::codon_to_neighbors_array.at(codon_seq.at(index / 9));
            // Codon after mutation, Nucleotide original and Nucleotide after mutation.
            char codon_to{0}, n_from{0}, n_to{0};

            tie(codon_to, n_from, n_to) = neighbors.at(index % 9);
            unsigned site = index / 9;

            Substitution substitution = {site, codon_seq.at(site), codon_to};
            substitutions.push_back(substitution);


            // Random shuffle of the fitness landscape
            uniform_real_distribution<double> unif_rand_proba(0, 1);
            double rand_uni = unif_rand_proba(Codon::re);
            if (rand_uni < proba_permutation) {
                shuffle(aa_fitness_profiles.at(site).begin(), aa_fitness_profiles.at(site).end(), Codon::re);
            }

            codon_seq.at(site) = codon_to;

        } else if (time_left < 0.) {
            time_left = 0.;
        }
        return time_left;
    }


    // Method computing all substitution event occuring during a given time-frame.
    // t: time during which substitution events occur (typically branch length).
    void run_substitutions(double t, vector<Substitution> &substitutions) {
        while (t > 0) {
            t = next_substitution(t, substitutions);
        }
    }


    // Set the the DNA sequence to the mutation-selection equilibrium.
    void at_equilibrium() {

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

        // For each site of the sequence.
        for (unsigned site{0}; site < nbr_sites; site++) {

            // Initialize the total sum of equilibrium frequency at 0.
            double total_frequencies{0.};

            // For each site of the vector of the site frequencies.
            for (char codon{0}; codon < 64; codon++) {
                double codon_freq{1.};

                // Translate the site to a triplet of DNA nucleotides
                array<char, 3> triplet = Codon::codon_to_triplet_array.at(codon);
                for (char position{0}; position < 3; position++) {
                    codon_freq *= nuc_frequencies(triplet.at(position));
                }

                if (Codon::codon_to_aa_array.at(codon) != 20) {
                    codon_frequencies.at(codon) = exp(aa_fitness_profiles.at(site).at(Codon::codon_to_aa_array.at(codon)));
                } else {
                    codon_frequencies.at(codon) = 0.;
                }


                // Increment the total sum of equilibrium frequencies.
                total_frequencies += codon_frequencies.at(codon);
            }

            // Uniform random generator between 0 and the total sum of equilibrium frequencies.
            uniform_real_distribution<double> unif(0, total_frequencies);
            double random_cumulative_frequencies = unif(Codon::re);
            double cumulative_frequencies{0.};

            char index{0};
            for (char m{0}; m < 64; m++) {
                // Iterate through the cumulative frequencies and break the loop when it is greater than the random cumulative frequencies.
                cumulative_frequencies += codon_frequencies.at(m);
                if (random_cumulative_frequencies < cumulative_frequencies) {
                    index = m;
                    break;
                }
            }

            // Substitute the site with the substitution given by the loop break.
            codon_seq.at(site) = index;
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

    // Method returning the DNA string corresponding to the codon sequence.
    string get_dna_str() const {

        // The DNA string is 3 times larger than the codon sequence.
        string dna_str(nbr_sites * 3, ' ');

        // For each site of the sequence.
        for (unsigned site{0}; site < nbr_sites; site++) {

            // Assert there is no stop in the sequence.
            assert(Codon::codon_to_aa_array.at(codon_seq.at(site)) != 20);

            // Translate the site to a triplet of DNA nucleotides
            array<char, 3> triplet = Codon::codon_to_triplet_array.at(codon_seq.at(site));
            for (char position{0}; position < 3; position++) {
                dna_str.at(3 * site + position) = Codon::nucleotides.at(triplet.at(position));
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
            name{move(name)}, length{stod(len)}, newick{move(newick)}, sequence_dna{move(seq)}, substitutions{0} {

        // Parse the newick tree descending of the node.
        parse_newick();
    }

    Node(string newick, unsigned const nbr_sites) :
            name{"Root"}, length{0.}, newick{move(newick)}, sequence_dna(nbr_sites), substitutions{0} {

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
    long nbr_substitutions(bool synonymous = false) {
        long nbr_substitutions{0};

        if (synonymous) {
            nbr_substitutions = count_if(substitutions.begin(), substitutions.end(), is_synonymous);
        } else {
            nbr_substitutions = substitutions.size();
        }

        if (!is_leaf()) {
            // Else, if the node is internal, iterate through the direct children.
            for (auto &child : children) {
                nbr_substitutions += child.nbr_substitutions(synonymous);
            }
        }
        return nbr_substitutions;
    }

    // Recursively iterate through the subtree.
    void traverse(string &output_filename) {
        // Substitutions of the DNA sequence is generated.
        sequence_dna.run_substitutions(length, substitutions);

        if (is_leaf()) {
            // If the node is a leaf, output the DNA sequence and name.
            ofstream output_file;
            output_file.open(output_filename, ios_base::app);
            output_file << name << " " << sequence_dna.get_dna_str() << endl;
            output_file.close();
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
                if (newick.at(position) == '(') {

                    size_t postpoint{position};
                    unsigned nbr_open{1};

                    for (size_t i{position + 1}; i < max_position; i++) {
                        if (nbr_open == 0) {
                            postpoint = i;
                            break;
                        } else if (newick.at(i) == '(') {
                            nbr_open++;
                        } else if (newick.at(i) == ')') {
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
    if (!input_stream) cerr << "Can't open newick file!";

    string line;
    getline(input_stream, line);

    return line;
}

vector<array<double, 20>> open_preferences(string const &file_name) {
    vector<array<double, 20>> fitness_profiles{0};

    ifstream input_stream(file_name);
    if (!input_stream) cerr << "Can't open preferences file!";

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
                fitness_profil.at(counter - 3) = log(stod(word));
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
      SimuEvol [--preferences=<file_path>] [--newick=<file_path>] [--output=<file_path>] [--mu=<0.5>] [--lambda=<3>] [--p=<0.01>]
      SimuEvol --help
      SimuEvol --version

Options:
-h --help                    show this help message and exit
--version                    show version and exit
--preferences=<file_path>    specify input site-specific preferences file .at(default: ../data/gal4.txt)
--newick=<file_path>         specify input newick tree .at(default: ../data/gal4.newick)
--output=<file_path>         specify output protein name .at(default: ../data/gal4.ali)
--mu=<0.5>                   specify the mutation rate .at(default: 0.5)
--lambda=<3>                 specify the strong to weak mutation bias .at(default: 3)
--p=<0.01>                   specify the probability to randomize the fitness landscape .at(default: 0.01)
)";

int main(int argc, char *argv[]) {

    auto args = docopt::docopt(USAGE,
                               {argv + 1, argv + argc},
                               true,               // show help if requested
                               "SimuEvol 0.1a");  // version string

    string preferences_path{"../data/gal4.txt"};
    if (args.at("--preferences")) {
        preferences_path = args.at("--preferences").asString();
    }

    string newick_path{"../data/gal4.newick"};
    if (args.at("--newick")) {
        newick_path = args.at("--newick").asString();
    }

    string output_path{"../data/gal4.ali"};
    if (args.at("--output")) {
        output_path = args.at("--output").asString();
    }

    vector<array<double, 20>> fitness_profiles = open_preferences(preferences_path);
    auto nbr_sites = static_cast<unsigned>(fitness_profiles.size());

    string newick_tree = open_newick(newick_path);
    Node root(newick_tree, nbr_sites);

    Matrix4x4 mutation_rate;
    double mu = 0.5;
    if (args.at("--mu")) {
        mu = stod(args.at("--mu").asString());
    }
    double lambda = 3.0;
    if (args.at("--lambda")) {
        lambda = stod(args.at("--lambda").asString());
    }
    mutation_rate << 0, 1, 1, lambda,
            lambda, 0, 1, lambda,
            lambda, 1, 0, lambda,
            lambda, 1, 1, 0;
    mutation_rate *= mu;
    mutation_rate -= mutation_rate.rowwise().sum().asDiagonal();

    ofstream output_file;
    output_file.open(output_path);
    output_file << root.nbr_leaves() << " " << nbr_sites * 3 << endl;
    output_file.close();

    cout << "The tree contains " << root.nbr_nodes() << " nodes for " << root.nbr_leaves() << " species at the tips."
         << endl;
    cout << "The tree has a total branch length of " << root.tot_length() << "." << endl;
    cout << "The DNA sequence is " << nbr_sites * 3 << " base pairs long." << endl;

    double p = 0.01;
    if (args.at("--p")) {
        p = stod(args.at("--p").asString());
    }

    root.set_evolution_parameters(mutation_rate, fitness_profiles, p);
    root.at_equilibrium();
    root.traverse(output_path);
    long nbr_substitutions = root.nbr_substitutions();
    long nbr_synonymous = root.nbr_substitutions(true);
    cout << "The simulation mapped " << nbr_substitutions << " substitutions along the tree." << endl;
    cout << nbr_synonymous << " synonymous and " << nbr_substitutions - nbr_synonymous
         << " non-synonymous substitutions." << endl;
    return 0;
}