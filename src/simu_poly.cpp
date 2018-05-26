#include <iostream>
#include <random>
#include <fstream>
#include <set>
#include "codon.hpp"

using namespace std;

#include "Eigen/Dense"

typedef Eigen::Matrix<double, 4, 4> Matrix4x4;
typedef Eigen::Matrix<double, 4, 1> Vector4x1;

#define DOCOPT_HEADER_ONLY

#include "docopt.cpp/docopt.h"

// Compute the kernel of the mutation-rate matrix.
// This kernel is a vector of the nucleotides frequencies (not normalized to 1) at equilibrium.
Vector4x1 equilibrium_frequencies(Matrix4x4 const &mutation_matrix) {
    Eigen::Matrix<double, 4, Eigen::Dynamic> kernel = mutation_matrix.transpose().fullPivLu().kernel();
    Vector4x1 nuc_frequencies;

    if (kernel.cols() > 1) {
        cerr << "The kernel has " << kernel.cols() << " dimensions, this is weird ! " << endl;
        uniform_int_distribution<unsigned> unif_unsigned(0, unsigned(kernel.cols()) - 1);
        unsigned chosen_row = unif_unsigned(Codon::re);
        nuc_frequencies = kernel.col(chosen_row);

    } else {
        nuc_frequencies = kernel.col(0);
    }
    nuc_frequencies /= nuc_frequencies.sum();

    return nuc_frequencies;
}


Matrix4x4 normalize_mutmatrix(Matrix4x4 mutation_matrix) {
    for (unsigned diag{0}; diag < 4; diag++) {
        mutation_matrix(diag, diag) = 0.0;
    }
    mutation_matrix -= mutation_matrix.rowwise().sum().asDiagonal();
    return mutation_matrix;
}

// Class representing DNA sequence.
class Reference {
private:
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

public:
    // The nbr_copies of sites in the sequence (each position is a codon, thus the DNA sequence is 3 times greater).
    unsigned const nbr_sites;
    // The nbr_copies of sites in the sequence (each position is a nucleotide).
    unsigned const nbr_nucleotides;
    // The sequence of codons.
    vector<char> codon_seq;

    // The matrix of mutation rates between nucleotides.
    Matrix4x4 mutation_rate_matrix;

    // The fitness profiles of amino-acids.
    vector<array<double, 20>> aa_fitness_profiles;

    array<unsigned, 4> nucleotides_count;
    array<vector<unsigned>, 4> nuc_rank_to_site;

    // Constructor of Reference_seq.
    // size: the size of the DNA sequence.
    explicit Reference(Matrix4x4 const &mutation_rate,
                       vector<array<double, 20>> const &fitness_profiles) :
            nbr_sites{unsigned(fitness_profiles.size())},
            nbr_nucleotides{unsigned(3 * fitness_profiles.size())},
            codon_seq(fitness_profiles.size(), 0),
            nucleotides_count{},
            nuc_rank_to_site{} {
        aa_fitness_profiles = fitness_profiles;
        mutation_rate_matrix = mutation_rate;
    }

    // Set the the DNA sequence to the mutation-selection equilibrium.
    void at_equilibrium() {
        nucleotides_count.fill(0);
        for (char nuc{0}; nuc < 4; nuc++) {
            nuc_rank_to_site[nuc].clear();
        }
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

            // Substitute the site with the mutation given by the loop break.
            codon_seq[site] = index;
            // For all nucleotides in the codon
        }
        refresh_count();
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

    void refresh_count() {
        nucleotides_count.fill(0);
        for (char nuc{0}; nuc < 4; nuc++) {
            nuc_rank_to_site[nuc].clear();
        }
        // For all site of the sequence.
        for (unsigned site{0}; site < nbr_sites; site++) {
            // For all nucleotides in the codon
            auto triplet = Codon::codon_to_triplet_array[codon_seq[site]];
            for (char position{0}; position < 3; position++) {
                char nuc = triplet[position];
                nucleotides_count[nuc]++;
                nuc_rank_to_site[nuc].push_back(3 * site + position);
            }
        }
        unsigned nbr_nuc{0}, nbr_positions{0};
        for (char nuc{0}; nuc < 4; nuc++) {
            nbr_nuc += nucleotides_count[nuc];
            nbr_positions += nuc_rank_to_site[nuc].size();
        }
        assert(nbr_positions == nbr_sites * 3);
        assert(nbr_nuc == nbr_sites * 3);
    }
};

struct Draw {
    char nuc_from{0};
    char nuc_to{0};
};

struct Mutation {
    char codon_from{0};
    char codon_to{0};
};

bool is_synonymous(Mutation const &s) {
    return (Codon::codon_to_aa_array[s.codon_from] == Codon::codon_to_aa_array[s.codon_to]);
}

class Haplotype {
private:

public:
    // The nbr_copies of sites in the sequence (each position is a codon, thus the DNA sequence is 3 times greater).
    unsigned nbr_copies;
    double fitness;
    map<unsigned, Mutation> mutations;
    array<unsigned, 4> nucleotides_count;
    array<vector<unsigned>, 4> nuc_rank_to_site;

    // Constructor of Reference_seq.
    // size: the size of the DNA sequence.
    explicit Haplotype(unsigned const nbr_copies, double const fitness, Reference const &reference) :
            nbr_copies{nbr_copies}, fitness{fitness}, mutations{},
            nucleotides_count{}, nuc_rank_to_site{} {
        cout << "I should be called only once ! Remember that " << endl;
        nucleotides_count = reference.nucleotides_count;
        nuc_rank_to_site = reference.nuc_rank_to_site;
    }

    void add_mutations(map<unsigned, Draw> &nuc_site_to_draw, Reference const &reference) {
        for (auto &kv_nuc_site_to_draw: nuc_site_to_draw) {
            char nuc_from = kv_nuc_site_to_draw.second.nuc_from;
            char nuc_to = kv_nuc_site_to_draw.second.nuc_to;

            unsigned site = kv_nuc_site_to_draw.first;
            unsigned codon_site = site / 3;
            unsigned nuc_position = site % 3;
            char codon_from{reference.codon_seq[codon_site]};
            array<char, 3> triplet_nuc = Codon::codon_to_triplet_array[codon_from];
            triplet_nuc[nuc_position] = nuc_to;
            char codon_to = Codon::triplet_to_codon(triplet_nuc[0], triplet_nuc[1], triplet_nuc[2]);
            if (mutations.count(codon_site) == 1) {
                // Found
                fitness -= reference.aa_fitness_profiles[codon_site][Codon::codon_to_aa_array[mutations[codon_site].codon_to]];
                mutations[codon_site].codon_to = codon_to;
            } else {
                // Not found
                fitness -= reference.aa_fitness_profiles[codon_site][Codon::codon_to_aa_array[codon_from]];
                mutations[codon_site] = Mutation{codon_from, codon_to};
            }
            nuc_rank_to_site[nuc_from].erase(
                    lower_bound(nuc_rank_to_site[nuc_from].begin(), nuc_rank_to_site[nuc_from].end(), site));
            nuc_rank_to_site[nuc_to].insert(
                    lower_bound(nuc_rank_to_site[nuc_to].begin(), nuc_rank_to_site[nuc_to].end(), site), site);
            nucleotides_count[nuc_from]--;
            nucleotides_count[nuc_to]++;
            fitness += reference.aa_fitness_profiles[codon_site][Codon::codon_to_aa_array[codon_to]];
        }
    }

    void check_consistency(Reference &reference) {
        unsigned nbr_nuc{0};
        set<unsigned> site_set{};
        for (char nuc{0}; nuc < 4; nuc++) {
            nbr_nuc += nucleotides_count[nuc];
            copy(nuc_rank_to_site[nuc].begin(), nuc_rank_to_site[nuc].end(), inserter(site_set, site_set.end()));
            assert(nuc_rank_to_site[nuc].size() == nucleotides_count[nuc]);
        }
        assert(site_set.size() == nbr_nuc);
        assert(nbr_nuc == reference.nbr_sites * 3);
    }
};

// Class representing nodes of a tree.
class Population {
private:
    unsigned population_size;
    vector<Haplotype> haplotype_vector;
    Reference reference_seq;

public:
    // Constructor
    Population(unsigned population_size, Reference &ref_seq) :
            population_size{population_size}, haplotype_vector{}, reference_seq{ref_seq} {
        haplotype_vector.emplace_back(Haplotype(population_size, 0.0, reference_seq));
    }

    void check_consistency() {
        assert(haplotype_vector.size() <= population_size);
        assert(!haplotype_vector.empty());
        cout << haplotype_vector.size() << " haplotypes in the population" << endl;

        unsigned nbr_copies{0};
        for (auto &haplotype: haplotype_vector) {
            haplotype.check_consistency(reference_seq);
            nbr_copies += haplotype.nbr_copies;
        }
        assert(nbr_copies == population_size);
    }

    void mutation() {
        map<unsigned, map<unsigned, map<unsigned, Draw>>> hap_to_copy_to_nuc_site_to_draw{};

        for (char nuc_from{0}; nuc_from < 4; nuc_from++) {
            vector<unsigned> nucleotide_copies_array(haplotype_vector.size(), 0);

            for (unsigned hap_id{0}; hap_id < haplotype_vector.size(); hap_id++) {
                nucleotide_copies_array[hap_id] =
                        haplotype_vector[hap_id].nucleotides_count[nuc_from] * haplotype_vector[hap_id].nbr_copies;
            }
            unsigned nucleotide_copies = accumulate(nucleotide_copies_array.begin(), nucleotide_copies_array.end(),
                                                    static_cast<unsigned>(0));
            assert(nucleotide_copies <= (population_size * reference_seq.nbr_nucleotides));

            array<double, 4> mutation_rate_array{0};
            for (char nuc_to{0}; nuc_to < 4; nuc_to++) {
                if (nuc_from == nuc_to) {
                    mutation_rate_array[nuc_to] = 0;
                } else {
                    mutation_rate_array[nuc_to] = reference_seq.mutation_rate_matrix(nuc_from, nuc_to);
                }
            }
            double mutation_rate = accumulate(mutation_rate_array.begin(), mutation_rate_array.end(), 0.0);
            assert(mutation_rate < 1.0);

            binomial_distribution<unsigned> binomial_distr(nucleotide_copies, mutation_rate);
            unsigned binomial_draw = binomial_distr(Codon::re);

            if (binomial_draw > 0) {
                discrete_distribution<char> nuc_distr(mutation_rate_array.begin(), mutation_rate_array.end());
                discrete_distribution<unsigned> haplotype_freq_distr(nucleotide_copies_array.begin(),
                                                                     nucleotide_copies_array.end());
                unsigned nbr_draws{0};

                while (nbr_draws < binomial_draw) {
                    char nuc_draw = nuc_distr(Codon::re);
                    unsigned haplotype_draw = haplotype_freq_distr(Codon::re);
                    unsigned haplotype_nucleotide_copies{nucleotide_copies_array[haplotype_draw]};

                    uniform_int_distribution<unsigned> copy_and_site_distr(0, haplotype_nucleotide_copies - 1);
                    unsigned copy_and_site_draw = copy_and_site_distr(Codon::re);

                    unsigned haplotype_nucleotide_count{haplotype_vector[haplotype_draw].nucleotides_count[nuc_from]};
                    unsigned nuc_rank = copy_and_site_draw % (haplotype_nucleotide_count);
                    unsigned nuc_site = haplotype_vector[haplotype_draw].nuc_rank_to_site[nuc_from][nuc_rank];
                    unsigned copy = copy_and_site_draw / haplotype_nucleotide_count;

                    bool flag{true};
                    if (hap_to_copy_to_nuc_site_to_draw.count(haplotype_draw) == 1) {
                        if (hap_to_copy_to_nuc_site_to_draw[haplotype_draw].count(copy) == 1) {
                            if (hap_to_copy_to_nuc_site_to_draw[haplotype_draw][copy].count(nuc_site) == 1) {
                                flag = false;
                            }
                        }
                    }
                    if (flag) {
                        hap_to_copy_to_nuc_site_to_draw[haplotype_draw][copy][nuc_site] = Draw{nuc_from, nuc_draw};
                        nbr_draws++;
                    }
                }
            }
        }

        for (auto &kv_hap_to_copy_to_nuc_site_to_draw: hap_to_copy_to_nuc_site_to_draw) {
            for (auto &kv_copy_to_nuc_site_to_draw: kv_hap_to_copy_to_nuc_site_to_draw.second) {
                Haplotype new_haplotype = haplotype_vector[kv_hap_to_copy_to_nuc_site_to_draw.first];
                haplotype_vector[kv_hap_to_copy_to_nuc_site_to_draw.first].nbr_copies--;
                new_haplotype.nbr_copies = 1;
                new_haplotype.add_mutations(kv_copy_to_nuc_site_to_draw.second, reference_seq);
                haplotype_vector.emplace_back(new_haplotype);
            }
        }
    }

    void selection_and_drift() {
        vector<double> haplotypes_fitness_array(haplotype_vector.size(), 0);
        for (unsigned hap_id{0}; hap_id < haplotype_vector.size(); hap_id++) {
            haplotypes_fitness_array[hap_id] =
                    haplotype_vector[hap_id].nbr_copies * abs(1 + haplotype_vector[hap_id].fitness);
            haplotype_vector[hap_id].nbr_copies = 0;
        }
        discrete_distribution<unsigned> haplotype_freq_distr(haplotypes_fitness_array.begin(),
                                                             haplotypes_fitness_array.end());
        for (unsigned indiv{0}; indiv < population_size; indiv++) {
            unsigned hap_id = haplotype_freq_distr(Codon::re);
            haplotype_vector[hap_id].nbr_copies++;
        }
    }

    void fixation() {
        // remove haplotypes with 0 copies
        auto iter = haplotype_vector.begin();
        while (iter != haplotype_vector.end()) {
            if (iter->nbr_copies == 0) {
                haplotype_vector.erase(iter);
            } else {
                iter++;
            }
        }

        // Fix the reference
        bool flag{false};
        for (unsigned site{0}; site < reference_seq.nbr_sites; site++) {
            unsigned copy{0};
            char codon_to{-1};
            for (auto haplotype: haplotype_vector) {
                if (haplotype.mutations.count(site) == 1) {
                    if (codon_to == -1) {
                        codon_to = haplotype.mutations[site].codon_to;
                        copy += haplotype.nbr_copies;
                    } else if (codon_to == haplotype.mutations[site].codon_to) {
                        copy += haplotype.nbr_copies;
                    } else {
                        break;
                    }
                } else {
                    break;
                }
            }
            if (copy == population_size) {
                flag = true;
                for (auto haplotype: haplotype_vector) {
                    haplotype.mutations.erase(site);
                }
                reference_seq.codon_seq[site] = codon_to;
            }
        }
        if (flag) {
            reference_seq.refresh_count();
        }

    }

    void forward() {
        mutation();
        selection_and_drift();
    }

    void run_forward(unsigned t_max) {
        check_consistency();
        for (unsigned time{0}; time < t_max; time++) {
            forward();
        }
        fixation();
        check_consistency();
    }
};

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
      SimuEvol [--preferences=<file_path>] [--mu=<1e-4>] [--lambda=<5.0>] [--pop_size=<2000>]
      SimuEvol --help
      SimuEvol --version

Options:
-h --help                    show this help message and exit
--version                    show version and exit
--preferences=<file_path>    specify input site-specific preferences file [default: ../data_prefs/np.txt]
--mu=<1e-4>                   specify the mutation rate [default: 1e-4]
--lambda=<5.0>                 specify the strong to weak mutation bias [default: 5.0]
--pop_size=<2000>                 specify the strong to weak mutation bias [default: 2000]
)";

int main(int argc, char *argv[]) {

    auto args = docopt::docopt(USAGE,
                               {argv + 1, argv + argc},
                               true,              // show help if requested
                               "SimuPoly 0.1a");  // version string

    string preferences_path{"../data_prefs/np.txt"};
    if (args["--preferences"]) {
        preferences_path = args["--preferences"].asString();
    }


    vector<array<double, 20>> fitness_profiles = open_preferences(preferences_path);

    double mu = 1e-6;
    if (args["--mu"]) {
        mu = stod(args["--mu"].asString());
    }

    double lambda = 5.0;
    if (args["--lambda"]) {
        lambda = stod(args["--lambda"].asString());
    }

    unsigned population_size = 1000;
    if (args["--pop_size"]) {
        population_size = static_cast<unsigned>(stoi(args["--pop_size"].asString()));
    }

    Matrix4x4 mutation_rate;
    mutation_rate << 0, 1, 1, lambda,
            lambda, 0, 1, lambda,
            lambda, 1, 0, lambda,
            lambda, 1, 1, 0;
    mutation_rate *= mu;
    mutation_rate = normalize_mutmatrix(mutation_rate);


    // output
    cout << "The DNA sequence is " << fitness_profiles.size() * 3 << " base pairs unsigned." << endl;
    cout << "The mutation transition matrix (" << Codon::nucleotides << ") is: " << endl;
    cout << mutation_rate << endl;

    Reference reference = Reference(mutation_rate, fitness_profiles);
    reference.at_equilibrium();
    cout << reference.get_dna_str() << endl;

    Population pop = Population(population_size, reference);
    pop.run_forward(1000);
    cout << reference.get_dna_str() << endl;
    return 0;
}