#include <iostream>
#include <random>
#include <fstream>
#include <set>
#include <iomanip>
#include "codon.hpp"

using namespace std;

#define DOCOPT_HEADER_ONLY

#include "docopt.cpp/docopt.h"

#include <chrono>

typedef chrono::high_resolution_clock::time_point TimeVar;

#define duration(a) chrono::duration_cast<std::chrono::nanoseconds>(a).count()
#define timeNow() chrono::high_resolution_clock::now()
double time_mutation{0.0}, time_selection{0.0}, time_extinction{0.0}, time_fixation{0.0}, time_export{0.0};
string version_string = "SimuPoly 0.1";
unsigned total_generation = 0;
typedef map<tuple<unsigned, unsigned, unsigned>, tuple<char, char>>::iterator IterMap;

//Function for sum
double sum(vector<unsigned> &v) {
    return accumulate(v.begin(), v.end(), static_cast<unsigned>(0));
}

//Function for average
double mean(vector<unsigned> &v) {
    double return_value = 0.0;

    for (int i = 0; i < v.size(); i++) {
        return_value += i * v[i];
    }

    return return_value / sum(v);
};

string join(vector<unsigned> &v, char sep) {
    return accumulate(v.begin() + 1, v.end(), to_string(v[0]),
                      [sep](const string &acc, int b) {
                          return acc + sep + to_string(b);
                      });
};

vector<tuple<char, char>> fixations_vector;
vector<tuple<char, char>> mutations_vector;

bool is_synonymous(tuple<char, char> &pair) {
    return Codon::codon_to_aa_array[get<0>(pair)] == Codon::codon_to_aa_array[get<1>(pair)];
};

double non_syn_fixation_bias() {
    unsigned min_last = 100000;
    unsigned last = min(static_cast<unsigned>(fixations_vector.size()), min_last);
    long nbr_syn_fixations = count_if(fixations_vector.end() - last, fixations_vector.end(), is_synonymous);
    long nbr_non_syn_fixation = last - nbr_syn_fixations;
    long nbr_syn_mutations = count_if(mutations_vector.end() - last, mutations_vector.end(), is_synonymous);
    long nbr_non_syn_mutations = last - nbr_syn_mutations;

    if (nbr_non_syn_mutations > 0 and nbr_syn_fixations > 0) {
        double bias = nbr_non_syn_fixation * nbr_syn_mutations;
        bias /= (nbr_syn_fixations * nbr_non_syn_mutations);
        return bias;
    } else {
        return 0.0;
    }
};

// Class representing DNA sequence.
class Reference {
private:
    // Method computing the equilibrium frequencies for one site.
    // aa_fitness_profil: The amino-acid fitness profil of the given site.

    array<double, 64> codon_frequencies(array<double, 20> const &aa_fitness_profil) {

        array<double, 64> codon_frequencies{0};

        // Initialize the total sum of equilibrium frequency at 0.
        double total_frequencies{0.};

        // For each site of the vector of the site frequencies.
        for (char codon{0}; codon < 64; codon++) {
            double codon_freq{1.};

            // For all nucleotides in the codon
            for (auto const &nuc: Codon::codon_to_triplet_array[codon]) {
                codon_freq *= nuc_frequencies[nuc];
            }

            if (Codon::codon_to_aa_array[codon] != 20) {
                codon_frequencies[codon] = codon_freq * exp(aa_fitness_profil[Codon::codon_to_aa_array[codon]]);
            } else {
                codon_frequencies[codon] = 0.;
            }


            // Increment the total sum of equilibrium frequencies.
            total_frequencies += codon_frequencies[codon];
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
    // The matrix of mutation rates between nucleotides.
    double lambda;
    array<array<double, 4>, 4> mutation_rate_matrix;
    array<double, 4> nuc_frequencies;

    // The fitness profiles of amino-acids.
    vector<array<double, 20>> aa_fitness_profiles;

    array<unsigned, 4> nucleotides_count;
    array<vector<unsigned>, 4> nuc_rank_to_site;

    // Constructor of Reference_seq.
    // size: the size of the DNA sequence.
    explicit Reference(double lambda, double mu, vector<array<double, 20>> const &fitness_profiles) :
            lambda{lambda}, nbr_sites{unsigned(fitness_profiles.size())},
            nbr_nucleotides{unsigned(3 * fitness_profiles.size())},
            codon_seq(fitness_profiles.size(), 0),
            nuc_frequencies{{lambda, 1, 1, lambda}},
            nucleotides_count{}, nuc_rank_to_site{}, mutation_rate_matrix{} {
        aa_fitness_profiles = fitness_profiles;

        mutation_rate_matrix[0] = {0, 1, 1, lambda};
        mutation_rate_matrix[1] = {lambda, 0, 1, lambda};
        mutation_rate_matrix[2] = {lambda, 1, 0, lambda};
        mutation_rate_matrix[3] = {lambda, 1, 1, 0};

        double sum_freq = accumulate(nuc_frequencies.begin(), nuc_frequencies.end(), 0.0);
        for (char nuc_from{0}; nuc_from < 4; nuc_from++) {
            nuc_frequencies[nuc_from] /= sum_freq;
            for (char nuc_to{0}; nuc_to < 4; nuc_to++) {
                mutation_rate_matrix[nuc_from][nuc_to] *= mu;
            }
        }
    }

    double mutation_rate() {
        double total = 0;
        for (char nuc_from{0}; nuc_from < 4; nuc_from++) {
            for (char nuc_to{0}; nuc_to < 4; nuc_to++) {
                total += nuc_frequencies[nuc_from] * mutation_rate_matrix[nuc_from][nuc_to];
            }
        }
    };

    // Set the the DNA sequence to the mutation-selection equilibrium.
    void at_equilibrium() {
        nucleotides_count.fill(0);
        for (char nuc{0}; nuc < 4; nuc++) {
            nuc_rank_to_site[nuc].clear();
        }
        // For all site of the sequence.
        for (unsigned site{0}; site < nbr_sites; site++) {

            array<double, 64> codon_freqs = codon_frequencies(aa_fitness_profiles[site]);
            discrete_distribution<char> freq_codon_distr(codon_freqs.begin(), codon_freqs.end());
            codon_seq[site] = freq_codon_distr(Codon::re);

        }
        refresh_count();
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

    void const set_parameters(Reference const &ref) {
        codon_seq = ref.codon_seq;
        nucleotides_count = ref.nucleotides_count;
        nuc_rank_to_site = ref.nuc_rank_to_site;
    };
};

class Haplotype {
public:
    // The nbr_copies of sites in the sequence (each position is a codon, thus the DNA sequence is 3 times greater).
    unsigned nbr_copies;
    double fitness;
    vector<char> codon_from_vector;
    vector<char> codon_to_vector;
    array<unsigned, 4> nucleotides_count;
    array<vector<unsigned>, 4> nuc_rank_to_site;

    // Constructor of Reference_seq.
    // size: the size of the DNA sequence.
    explicit Haplotype(unsigned const nbr_copies, double const fitness, Reference const &reference) :
            nbr_copies{nbr_copies}, fitness{fitness},
            codon_from_vector{reference.codon_seq},
            codon_to_vector{reference.codon_seq},
            nucleotides_count{}, nuc_rank_to_site{} {
        nucleotides_count = reference.nucleotides_count;
        nuc_rank_to_site = reference.nuc_rank_to_site;
    }

    char mutate_codon(char codon_from, char nuc_position, char nuc_to) {
        array<char, 3> triplet_nuc = Codon::codon_to_triplet_array[codon_from];
        triplet_nuc[nuc_position] = nuc_to;
        return Codon::triplet_to_codon(triplet_nuc[0], triplet_nuc[1], triplet_nuc[2]);
    }

    void add_mutations(IterMap iter, IterMap end, Reference const &reference) {
        while (iter != end) {
            char nuc_from = get<0>(iter->second);
            char nuc_to = get<1>(iter->second);
            assert(nuc_from != nuc_to);
            unsigned site = get<2>(iter->first);
            unsigned codon_site = site / 3;
            auto nuc_position = static_cast<char>(site % 3);

            char codon_from = codon_to_vector[codon_site];
            char codon_to = mutate_codon(codon_from, nuc_position, nuc_to);
            if (Codon::codon_to_aa_array[codon_to] != 20) {
                mutations_vector.emplace_back(make_tuple(codon_from, codon_to));
                codon_from_vector[codon_site] = codon_from;
                codon_to_vector[codon_site] = codon_to;

                nuc_rank_to_site[nuc_from].erase(
                        lower_bound(nuc_rank_to_site[nuc_from].begin(), nuc_rank_to_site[nuc_from].end(), site));
                nuc_rank_to_site[nuc_to].insert(
                        lower_bound(nuc_rank_to_site[nuc_to].begin(), nuc_rank_to_site[nuc_to].end(), site), site);
                nucleotides_count[nuc_from]--;
                nucleotides_count[nuc_to]++;
                fitness -= reference.aa_fitness_profiles[codon_site][Codon::codon_to_aa_array[codon_from]];
                fitness += reference.aa_fitness_profiles[codon_site][Codon::codon_to_aa_array[codon_to]];
            }
            iter++;
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
        for (unsigned site{0}; site < reference.nbr_sites; site++) {
            assert(Codon::codon_to_aa_array[codon_to_vector[site]] != 20);
            assert(Codon::codon_to_aa_array[codon_from_vector[site]] != 20);
        }
        assert(site_set.size() == nbr_nuc);
        assert(nbr_nuc == reference.nbr_sites * 3);
    }

    static struct {
        bool operator()(const Haplotype &left, const Haplotype &right) {
            return left.nbr_copies > right.nbr_copies;
        }

        bool operator()(const Haplotype &left, float right) {
            return left.nbr_copies > right;
        }

        bool operator()(float left, const Haplotype &right) {
            return left > right.nbr_copies;
        }
    } GreaterThan;

};

// Class representing a population
class Population {
public:
    unsigned population_size;
    Reference reference_seq;
    unsigned elapsed;
    double beta;
    vector<Haplotype> haplotype_vector;

    // Constructor
    Population(double mutation_rate, double mu, vector<array<double, 20>> const &fitness_profiles,
               unsigned population_size, double beta, string &output_path) :
            population_size{population_size}, beta{beta}, haplotype_vector{},
            reference_seq(mutation_rate, mu, fitness_profiles), elapsed{0} {
        reference_seq.at_equilibrium();
        haplotype_vector.emplace_back(Haplotype(population_size, 0.0, reference_seq));
        check_consistency();
        for (unsigned time{0}; time < population_size * 2; time++) {
            forward();
        }
        mutations_vector.clear();
        fixations_vector.clear();
        cout << "Population created " << endl;
        assert(!haplotype_vector.empty());
    }

    void set_parameters(Population &pop) {
        reference_seq.set_parameters(pop.reference_seq);
        haplotype_vector = pop.haplotype_vector;
        elapsed = pop.elapsed;
    };

    void check_consistency() {
        assert(haplotype_vector.size() <= population_size);
        assert(!haplotype_vector.empty());

        unsigned nbr_copies{0};
        for (auto &haplotype: haplotype_vector) {
            haplotype.check_consistency(reference_seq);
            nbr_copies += haplotype.nbr_copies;
        }
        assert(nbr_copies == population_size);
    }

    void forward() {
        mutation();
        selection_and_drift();
        extinction();
    }

    void run_forward(unsigned t_max) {
        assert(!haplotype_vector.empty());
        for (unsigned time{0}; time <= t_max; time++) {
            forward();
            if (time % population_size == 0) {
                fixation();
            }
        }
        if (t_max > 0) {
            fixation();
        }
        check_consistency();
        elapsed += t_max;
    }

    void mutation() {
        TimeVar t_start = timeNow();
        map<tuple<unsigned, unsigned, unsigned>, tuple<char, char>> cood_to_mutation{};

        for (char nuc_from{0}; nuc_from < 4; nuc_from++) {
            vector<unsigned> nucleotide_copies_array(haplotype_vector.size(), 0);

            for (unsigned hap_id{0}; hap_id < haplotype_vector.size(); hap_id++) {
                nucleotide_copies_array[hap_id] =
                        haplotype_vector[hap_id].nucleotides_count[nuc_from] * haplotype_vector[hap_id].nbr_copies;
            }
            unsigned nucleotide_copies = accumulate(nucleotide_copies_array.begin(), nucleotide_copies_array.end(),
                                                    static_cast<unsigned>(0));
            assert(nucleotide_copies <= (population_size * reference_seq.nbr_nucleotides));

            double mutation_rate = accumulate(reference_seq.mutation_rate_matrix[nuc_from].begin(),
                                              reference_seq.mutation_rate_matrix[nuc_from].end(), 0.0);
            assert(mutation_rate < 1.0);
            binomial_distribution<unsigned> binomial_distr(nucleotide_copies, mutation_rate);
            unsigned binomial_draw = binomial_distr(Codon::re);

            if (binomial_draw > 0) {
                discrete_distribution<char> nuc_distr(reference_seq.mutation_rate_matrix[nuc_from].begin(),
                                                      reference_seq.mutation_rate_matrix[nuc_from].end());
                discrete_distribution<unsigned> haplotype_freq_distr(nucleotide_copies_array.begin(),
                                                                     nucleotide_copies_array.end());
                unsigned nbr_draws{0};

                while (nbr_draws < binomial_draw) {
                    unsigned haplotype_draw = haplotype_freq_distr(Codon::re);
                    unsigned haplotype_nucleotide_copies = nucleotide_copies_array[haplotype_draw];

                    uniform_int_distribution<unsigned> copy_and_site_distr(0, haplotype_nucleotide_copies - 1);
                    unsigned copy_and_site_draw = copy_and_site_distr(Codon::re);

                    unsigned haplotype_nucleotide_count = haplotype_vector[haplotype_draw].nucleotides_count[nuc_from];
                    unsigned nuc_rank = copy_and_site_draw % (haplotype_nucleotide_count);
                    unsigned nuc_site = haplotype_vector[haplotype_draw].nuc_rank_to_site[nuc_from][nuc_rank];
                    unsigned copy = copy_and_site_draw / haplotype_nucleotide_count;

                    auto cood = make_tuple(haplotype_draw, copy, nuc_site);
                    if (cood_to_mutation.count(cood) == 0) {
                        cood_to_mutation[cood] = make_tuple(nuc_from, nuc_distr(Codon::re));
                        nbr_draws++;
                    }
                }
            }
        }

        auto begin = cood_to_mutation.begin();
        auto next = cood_to_mutation.begin();
        while (begin != cood_to_mutation.end()) {
            unsigned hap_id = get<0>(begin->first);
            while (true) {
                next++;
                if (next == cood_to_mutation.end()) {
                    break;
                } else {
                    if ((get<0>(next->first) != hap_id) or (get<1>(next->first) != get<1>(begin->first))) {
                        break;
                    }
                }
            }
            Haplotype new_haplotype = haplotype_vector[hap_id];
            haplotype_vector[hap_id].nbr_copies--;
            new_haplotype.nbr_copies = 1;
            new_haplotype.add_mutations(begin, next, reference_seq);
            haplotype_vector.emplace_back(new_haplotype);
            begin = next;
        }
        time_mutation += duration(timeNow() - t_start);
    }

    void selection_and_drift() {
        TimeVar t_start = timeNow();

        vector<double> hap_fit_array(haplotype_vector.size(), 0);
        for (unsigned hap_id{0}; hap_id < haplotype_vector.size(); hap_id++) {
            double fitness = max(0.0, 1.0 + beta * haplotype_vector[hap_id].fitness / population_size);
            hap_fit_array[hap_id] = haplotype_vector[hap_id].nbr_copies * fitness;
            haplotype_vector[hap_id].nbr_copies = 0;
        }
        discrete_distribution<unsigned> haplotype_freq_distr(hap_fit_array.begin(), hap_fit_array.end());
        for (unsigned indiv{0}; indiv < population_size; indiv++) {
            haplotype_vector[haplotype_freq_distr(Codon::re)].nbr_copies++;
        }
        total_generation++;
        time_selection += duration(timeNow() - t_start);
    }

    void extinction() {
        TimeVar t_start = timeNow();

        // remove haplotypes with 0 copies
        sort(haplotype_vector.begin(), haplotype_vector.end(), Haplotype::GreaterThan);
        auto low_bound = lower_bound(haplotype_vector.begin(), haplotype_vector.end(), 0, Haplotype::GreaterThan);
        haplotype_vector.erase(low_bound, haplotype_vector.end());

        time_extinction += duration(timeNow() - t_start);
    }

    void fixation() {
        TimeVar t_start = timeNow();
        map<tuple<char, char>, unsigned> codon_from_to_copy{};
        for (unsigned site{0}; site < reference_seq.nbr_sites; site++) {
            codon_from_to_copy.clear();
            for (auto const &haplotype: haplotype_vector) {
                char codon_from = haplotype.codon_from_vector[site];
                char codon_to = haplotype.codon_to_vector[site];
                if (codon_to != codon_from) {
                    codon_from_to_copy[make_tuple(codon_from, codon_to)] += haplotype.nbr_copies;
                } else {
                    break;
                }
            }

            if (codon_from_to_copy.size() == 1) {
                for (auto const &kv_tuple_copy: codon_from_to_copy) {
                    if (kv_tuple_copy.second == population_size) {
                        char codon_from = get<0>(kv_tuple_copy.first);
                        char codon_to = get<1>(kv_tuple_copy.first);
                        assert(codon_to != codon_from);
                        for (auto &haplotype: haplotype_vector) {
                            haplotype.codon_from_vector[site] = codon_to;
                        }
                        fixations_vector.emplace_back(make_tuple(codon_from, codon_to));
                        reference_seq.codon_seq[site] = codon_to;
                    }
                }
            }
        }

        time_fixation += duration(timeNow() - t_start);
    }

    void output_vcf(string &output_filename, string &node_name) {
        TimeVar t_start = timeNow();

        vector<unsigned> sfs_non_syn(population_size, 0), sfs_syn(population_size, 0);
        unsigned non_syn_nbr = 0, syn_nbr = 0;
        unsigned complex_sites{0};
        assert(sfs_non_syn.size() == population_size);
        assert(sum(sfs_non_syn) == 0);
        assert(sfs_syn.size() == population_size);
        assert(sum(sfs_syn) == 0);

        ofstream vcf_file;
        vcf_file.open(output_filename + node_name + ".vcf");
        vcf_file << "##fileformat=VCFv4.0" << endl;
        vcf_file << "##source=" << version_string << endl;
        vcf_file << "##nodeName=" << node_name << endl;
        vcf_file << "##sequenceSize=" << reference_seq.nbr_sites << endl;
        vcf_file << "##populationSize=" << population_size << endl;
        vcf_file << "##numberHaplotypes=" << haplotype_vector.size() << endl;
        vcf_file << "##reference=" << get_dna_str() << endl;
        vcf_file << "##FILTER=<ID=stop,Description=\"The mutation is a stop codon\">" << endl;
        vcf_file << "##INFO=<ID=REFCODON,Number=1,Type=String,Description=\"Codon of the reference\">" << endl;
        vcf_file << "##INFO=<ID=ALTCODON,Number=1,Type=String,Description=\"Codon of the alternative\">" << endl;
        vcf_file << "##INFO=<ID=REFAA,Number=1,Type=Character,Description=\"Amino-acid of the reference\">" << endl;
        vcf_file << "##INFO=<ID=ALTAA,Number=1,Type=Character,Description=\"Amino-acid of the alternative\">" << endl;
        vcf_file << "##INFO=<ID=POSITION,Number=1,Type=Integer,Description=\"Mutated codon position\">" << endl;
        vcf_file << "##INFO=<ID=ALTCOUNT,Number=1,Type=Integer,Description=\"Number of alternative copy\">" << endl;
        vcf_file << "##INFO=<ID=SYN,Number=1,Type=boolean,Description=\"Is a synonymous mutation\">" << endl;
        vcf_file << "##FORMAT=<ID=GT,Number=1,Type=String,Description=\"Genotype\">" << endl;
        string header{"#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT"};
        unsigned indiv{1};
        for (unsigned hap_id{0}; hap_id < haplotype_vector.size(); hap_id++) {
            for (unsigned copy{0}; copy < haplotype_vector[hap_id].nbr_copies; copy++) {
                header += "\tId";
                header += to_string(indiv);
                header += "Hap";
                header += to_string(hap_id);
                indiv++;
            }
        }
        vcf_file << header << endl;

        for (unsigned site{0}; site < reference_seq.nbr_sites; site++) {
            map<tuple<char, char>, unsigned> codon_from_to_copy{};
            for (auto const &haplotype: haplotype_vector) {
                char codon_from = haplotype.codon_from_vector[site];
                char codon_to = haplotype.codon_to_vector[site];
                if (codon_to != codon_from) {
                    codon_from_to_copy[make_tuple(codon_from, codon_to)] += haplotype.nbr_copies;
                }
            }

            if (codon_from_to_copy.size() == 1) {
                for (auto const &kv_tuple_copy: codon_from_to_copy) {
                    char codon_from = get<0>(kv_tuple_copy.first);
                    char codon_to = get<1>(kv_tuple_copy.first);
                    if (codon_to != codon_from) {
                        char position{0};
                        while (position < 3) {
                            if (Codon::codon_to_nuc(codon_from, position) != Codon::codon_to_nuc(codon_to, position)) {
                                break;
                            } else {
                                position++;
                            }
                        }
                        assert(position != 3);
                        string line{};
                        line += ".\t";
                        line += to_string(3 * site + position);
                        line += "\t.\t";
                        line += Codon::codon_to_nuc(codon_from, position);
                        line += "\t";
                        line += Codon::codon_to_nuc(codon_to, position);
                        line += "\t100\t";
                        if (Codon::codon_to_aa_array[codon_to] == 20) {
                            line += "stop";
                        } else {
                            line += "PASS";
                        }
                        line += "\tREFCODON=";
                        line += Codon::codon_string(codon_from);
                        line += ";ALTCODON=";
                        line += Codon::codon_string(codon_to);
                        line += ";REFAA=";
                        line += Codon::codon_aa_string(codon_from);
                        line += ";ALTAA=";
                        line += Codon::codon_aa_string(codon_to);
                        line += ";POSITION=";
                        line += to_string(position);
                        line += ";ALTCOUNT=";
                        line += to_string(kv_tuple_copy.second);
                        line += ";SYN=";

                        assert(kv_tuple_copy.second < population_size);
                        if (Codon::codon_to_aa_array[codon_from] == Codon::codon_to_aa_array[codon_to]) {
                            sfs_syn[kv_tuple_copy.second]++;
                            syn_nbr++;
                            line += "TRUE";
                        } else {
                            sfs_non_syn[kv_tuple_copy.second]++;
                            non_syn_nbr++;
                            line += "FALSE";
                        }

                        line += "\tGT";
                        for (auto const &haplotype: haplotype_vector) {
                            for (unsigned copy{0}; copy < haplotype.nbr_copies; copy++) {
                                line += "\t";
                                line += Codon::codon_to_nuc(haplotype.codon_to_vector[site], position);
                            }
                        }
                        vcf_file << line << endl;
                    }
                }
            } else if (codon_from_to_copy.size() > 1) {
                complex_sites++;
            }
        }

        vcf_file.close();

        assert(sfs_non_syn.size() == population_size);
        assert(sfs_syn.size() == population_size);
        assert(sum(sfs_non_syn) == non_syn_nbr);
        assert(sum(sfs_syn) == syn_nbr);
        double fixation_bias = non_syn_fixation_bias();
        double mean_fitness = accumulate(haplotype_vector.begin(), haplotype_vector.end(), 0.0,
                                         [](double acc, Haplotype &h) { return acc + h.fitness; }) /
                              haplotype_vector.size();
        // .ali format
        ofstream tsv_file;
        tsv_file.open(output_filename + ".tsv", ios_base::app);
        tsv_file << node_name << "\t" << elapsed << "\tFixedNbr\t" << fixations_vector.size() << endl;
        tsv_file << node_name << "\t" << elapsed << "\tNonSynFixBias\t" << fixation_bias << endl;
        tsv_file << node_name << "\t" << elapsed << "\tMeanFitness\t" << mean_fitness << endl;
        tsv_file << node_name << "\t" << elapsed << "\tHapNbr\t" << haplotype_vector.size() << endl;
        tsv_file << node_name << "\t" << elapsed << "\tCpxSites\t" << complex_sites << endl;
        tsv_file << node_name << "\t" << elapsed << "\tPn\t" << non_syn_nbr << endl;
        tsv_file << node_name << "\t" << elapsed << "\tSFSn\t" << join(sfs_non_syn, ' ') << endl;
        tsv_file << node_name << "\t" << elapsed << "\tE[SFSn]\t" << mean(sfs_non_syn) << endl;
        tsv_file << node_name << "\t" << elapsed << "\tPs\t" << syn_nbr << endl;
        tsv_file << node_name << "\t" << elapsed << "\tSFSs\t" << join(sfs_syn, ' ') << endl;
        tsv_file << node_name << "\t" << elapsed << "\tE[SFSs]\t" << mean(sfs_syn) << endl;
        tsv_file.close();

        time_export += duration(timeNow() - t_start);
    }


    // Method returning the DNA string corresponding to the codon sequence.
    string get_dna_str() const {

        // The DNA string is 3 times larger than the codon sequence.
        string dna_str(reference_seq.nbr_sites * 3, ' ');

        // For each site of the sequence.
        for (unsigned site{0}; site < reference_seq.nbr_sites; site++) {

            // Assert there is no stop in the sequence.
            assert(Codon::codon_to_aa_array[reference_seq.codon_seq[site]] != 20);

            // Translate the site to a triplet of DNA nucleotides
            array<char, 3> triplet = Codon::codon_to_triplet_array[reference_seq.codon_seq[site]];
            for (char position{0}; position < 3; position++) {
                dna_str[3 * site + position] = Codon::nucleotides[triplet[position]];
            }
        }
        return dna_str; // return the DNA sequence as a string.
    }
};

double length_computed = 0;

// Class representing nodes of a tree.
class Node {
private:
    string name;  // The species name of the node.
    double length;  // The length of the branch attached to the node (ascending).
    string newick;  // The newick tree descending from the node.
    vector<Node> children;  // Vector of direct children (first order, no grand-children).
    Population population;  // The DNA sequence attached to the node.

public:

    // Constructor
    Node(string name, string const &len, string newick, Population &pop) :
            name{move(name)}, length{stod(len)}, newick{move(newick)}, population{pop} {

        // Parse the newick tree descending of the node.
        parse_newick();
    }

    Node(string newick, Population &pop) :
            name{"Root"}, length{0.}, newick{move(newick)}, population{pop} {

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

    // Recursively iterate through the subtree.
    void traverse(string &output_filename) {
        // Substitutions of the DNA sequence is generated.
        auto time = static_cast<unsigned>(length / (population.reference_seq.mutation_rate()));
        population.run_forward(time);

        length_computed += length;
        cout << length_computed << " length computed" << endl;

        if (is_leaf()) {
            // If the node is a leaf, output the DNA sequence and name.
            string dna_str = population.get_dna_str();

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

            population.output_vcf(output_filename, name);

        } else {
            // If the node is internal, iterate through the direct children.
            for (auto &child : children) {
                child.population.set_parameters(population);
                child.traverse(output_filename);
            }
        }
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
                add_child(Node(name, length, subtree, population));
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
      SimuPoly [--preferences=<file_path>] [--newick=<file_path>] [--mu=<1e-7>] [--lambda=<5.0>] [--pop_size=<1000>] [--beta=<1.0>]
      SimuPoly --help
      SimuPoly --version

Options:
-h --help                    show this help message and exit
--version                    show version and exit
--preferences=<file_path>    specify input site-specific preferences file [default: ../data_prefs/np.txt]
--newick=<file_path>         specify input newick tree [default: ../data/np.newick]
--mu=<1e-7>                  specify the mutation rate [default: 1e-7]
--lambda=<5.0>               specify the strong to weak mutation bias [default: 5.0]
--pop_size=<1000>            specify the population size [default: 1000]
--beta=<1.0>                 specify the strength of selection [default: 1.0]
)";

int main(int argc, char *argv[]) {

    auto args = docopt::docopt(USAGE,
                               {argv + 1, argv + argc},
                               true,              // show help if requested
                               version_string);  // version string

    string preferences_path{"../data_prefs/np.txt"};
    if (args["--preferences"]) {
        preferences_path = args["--preferences"].asString();
    }
    vector<array<double, 20>> fitness_profiles = move(open_preferences(preferences_path));

    string newick_path{"../data_trees/gal4.newick"};
    if (args["--newick"]) {
        newick_path = args["--newick"].asString();
    }
    string newick_tree = open_newick(newick_path);

    double mu = 1e-6;
    if (args["--mu"]) {
        mu = stod(args["--mu"].asString());
    }

    double lambda = 5.0;
    if (args["--lambda"]) {
        lambda = stod(args["--lambda"].asString());
    }

    unsigned pop_size = 1000;
    if (args["--pop_size"]) {
        pop_size = static_cast<unsigned>(stoi(args["--pop_size"].asString()));
    }

    double beta = 1.0;
    if (args["--beta"]) {
        mu = stod(args["--beta"].asString());
    }

    // output
    cout << "The DNA sequence is " << fitness_profiles.size() * 3 << " base pairs." << endl;

    string output_path{"../data_sfs/np"};
    if (args["--output"]) {
        output_path = args["--output"].asString();
    }
    output_path += to_string(mu) + "_" + to_string(pop_size) + "_" + to_string(lambda) + "_" + to_string(beta);
    ofstream txt_file;
    txt_file.open(output_path + ".tsv");
    txt_file << "mu\t" << mu << endl;
    txt_file << "ne\t" << pop_size << endl;
    txt_file << "lambda\t" << lambda << endl;
    txt_file << "n\t" << fitness_profiles.size() * 3 << endl;
    txt_file.close();

    Population population(lambda, mu, fitness_profiles, pop_size, beta, output_path);
    Node root(newick_tree, population);
    cout << "The tree has " << root.nbr_leaves() << " leaves"<< endl;
    cout << "The tree is of length " << root.tot_length() << endl;
    ofstream ali_file;
    ali_file.open(output_path + ".ali");
    ali_file << root.nbr_leaves() << " " << fitness_profiles.size() * 3 << endl;
    ali_file.close();

    ofstream fasta_file;
    fasta_file.open(output_path + ".fasta");
    fasta_file.close();

    ofstream tsv_file;
    tsv_file.open(output_path + ".fasta");
    tsv_file.close();

    root.traverse(output_path);

    cout << total_generation << " generations computed" << endl;
    cout << fixations_vector.size() << " fixations" << endl;
    cout << mutations_vector.size() << " mutations" << endl;
    cout << "Non-synonymous fixation bias " << non_syn_fixation_bias() << endl;

    double total_time = time_mutation + time_selection + time_extinction + time_fixation + time_export;
    cout << setprecision(3) << total_time / 1e9 << "s total time" << endl;
    cout << 100 * time_mutation / total_time << "% of time spent in calculating mutation ("
         << time_mutation / 1e9 << "s)" << endl;
    cout << 100 * time_selection / total_time << "% of time spent in calculating selection ("
         << time_selection / 1e9 << "s)" << endl;
    cout << 100 * time_extinction / total_time << "% of time spent in calculating extinction ("
         << time_extinction / 1e9 << "s)" << endl;
    cout << 100 * time_fixation / total_time << "% of time spent in calculating fixation ("
         << time_fixation / 1e9 << "s)" << endl;
    cout << 100 * time_export / total_time << "% of time spent in exporting vcf ("
         << time_export / 1e9 << "s)" << endl;
    return 0;
}