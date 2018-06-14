#include <iostream>
#include <random>
#include <fstream>
#include <set>
#include <iomanip>
#include "codon.hpp"

using namespace std;

string version_string = "SimuPoly 0.1";

#define DOCOPT_HEADER_ONLY

#include "docopt.cpp/docopt.h"

#include <chrono>

#define duration(a) chrono::duration_cast<std::chrono::nanoseconds>(a).count()
#define timeNow() chrono::high_resolution_clock::now()

typedef chrono::high_resolution_clock::time_point TimeVar;

//Function for sum
double sum(vector<unsigned> const &v) {
    return accumulate(v.begin(), v.end(), static_cast<unsigned>(0));
}

//Function for mean
double mean(vector<unsigned> const &v) {
    double return_value = 0.0;

    for (unsigned i = 0; i < v.size(); i++) {
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

struct Change {
    unsigned long id;
    char codon_from;
    char codon_to;
};

struct Change_table {
    long non_syn_mut;
    long syn_mut;
    long non_syn_fix;
    long syn_fix;
};

double dn_from_table(Change_table const &table) {
    auto dn = static_cast<double>(table.non_syn_fix);
    dn /= table.non_syn_mut;
    return dn;
}

double ds_from_table(Change_table const &table) {
    auto ds = static_cast<double>(table.syn_fix);
    ds /= table.syn_mut;
    return ds;
};

struct TimeElapsed {
    double mutation;
    double selection;
    double extinction;
    double fixation;
    double exportation;
};

struct SummaryStatistics {
    double dn{0}, ds{0}, omega{0}, omega_predicted;
    double dn_sample{0}, ds_sample{0}, omega_sample{0};
    double pn{0}, ps{0}, pnps{0};
    double pn_sample{0}, ps_sample{0}, pnps_sample{0};
};

bool isfinite(SummaryStatistics const &s) {
    return (isfinite(s.omega) and isfinite(s.omega_sample));
};

//Function for avg of a vector ignoring NaN
double average(vector<SummaryStatistics> const &vector_stats, function<double(SummaryStatistics)> func) {
    double total{0};
    unsigned elements_non_nan{0};
    for (auto &stat: vector_stats) {
        double x = func(stat);
        if (isfinite(x)) {
            total += x;
            elements_non_nan++;
        }
    }
    return total / elements_non_nan;
}

SummaryStatistics average(vector<SummaryStatistics> const &vector_stats) {
    SummaryStatistics avg_stats;
    avg_stats.dn = average(vector_stats, [](SummaryStatistics stat) { return stat.dn; });
    avg_stats.ds = average(vector_stats, [](SummaryStatistics stat) { return stat.ds; });
    avg_stats.omega = average(vector_stats, [](SummaryStatistics stat) { return stat.omega; });
    avg_stats.dn_sample = average(vector_stats, [](SummaryStatistics stat) { return stat.dn_sample; });
    avg_stats.ds_sample = average(vector_stats, [](SummaryStatistics stat) { return stat.ds_sample; });
    avg_stats.omega_sample = average(vector_stats, [](SummaryStatistics stat) { return stat.omega_sample; });
    avg_stats.omega_predicted = average(vector_stats, [](SummaryStatistics stat) { return stat.omega_predicted; });
    return avg_stats;
}

bool operator<(const Change &lhs, const Change &rhs) {
    return lhs.id < rhs.id;
}

bool operator==(const Change &lhs, const Change &rhs) {
    return lhs.id == rhs.id;
}

bool is_synonymous(const Change &change) {
    return Codon::codon_to_aa_array[change.codon_from] == Codon::codon_to_aa_array[change.codon_to];
};

class Haplotype {
public:
    // The nbr_copies of sites in the sequence (each position is a codon, thus the DNA sequence is 3 times greater).
    unsigned nbr_copies{0};
    double fitness{0.0};
    map<unsigned, Change> change_map;

    // Constructor of Reference_seq.
    // size: the size of the DNA sequence.
    Haplotype(unsigned const nbr_copies,
              double const fitness) :
            nbr_copies{nbr_copies}, fitness{fitness} {

    };

    explicit Haplotype() = default;

    void check_consistency(unsigned nbr_sites) {
        for (auto &kv_change: change_map) {
            assert(Codon::codon_to_aa_array[kv_change.second.codon_to] != 20);
            assert(Codon::codon_to_aa_array[kv_change.second.codon_from] != 20);
        }
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
    // TimeElapsed
    unsigned elapsed{0};

    // Population
    unsigned population_size;

    // Mutation
    double lambda;
    array<array<double, 4>, 4> mutation_rate_matrix;
    array<double, 4> nuc_frequencies;

    // Selection
    vector<array<double, 20>> aa_fitness_profiles;

    // Reference sequence
    unsigned const nbr_sites;
    unsigned const nbr_nucleotides;
    vector<char> codon_seq;

    // Haplotypes
    vector<Haplotype> haplotype_vector;

    // Keep track of all changes
    vector<Change> fixations_vector;
    vector<Change> mutations_vector;

    // Statics variables (shared by all instances)
    static TimeElapsed time;
    static unsigned total_generation;
    static unsigned long nbr_mutations, nbr_fixations;
    static vector<SummaryStatistics> stats_vector;

    // Constructor
    Population(double mutation_bias, double mu, vector<array<double, 20>> const &fitness_profiles,
               unsigned population_size, string &output_path) :
            population_size{population_size},
            lambda{mutation_bias}, mutation_rate_matrix{}, nuc_frequencies{{mutation_bias, 1, 1, mutation_bias}},
            aa_fitness_profiles{fitness_profiles},
            nbr_sites{unsigned(fitness_profiles.size())}, nbr_nucleotides{unsigned(3 * fitness_profiles.size())},
            codon_seq(fitness_profiles.size(), 0),
            haplotype_vector{},
            fixations_vector{}, mutations_vector{} {

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

        // Draw codon from codon frequencies
        for (unsigned site{0}; site < nbr_sites; site++) {
            array<double, 64> codon_freqs = codon_frequencies(aa_fitness_profiles[site]);
            discrete_distribution<char> freq_codon_distr(codon_freqs.begin(), codon_freqs.end());
            codon_seq[site] = freq_codon_distr(Codon::re);

        }

        haplotype_vector.emplace_back(Haplotype(population_size, 0.0));
        check_consistency();
        for (unsigned gen{0}; gen < population_size * 2; gen++) {
            forward();
        }
        cout << "Population created " << endl;
        assert(!haplotype_vector.empty());
    }

    void set_parameters(Population const &pop) {
        codon_seq = pop.codon_seq;
        haplotype_vector = pop.haplotype_vector;
        elapsed = pop.elapsed;
        mutations_vector.clear();
        fixations_vector.clear();
    };

    // Method computing the equilibrium frequencies for one site.
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
        }

        return codon_frequencies;
    };

    double mutation_rate() {
        double total = 0;
        for (char nuc_from{0}; nuc_from < 4; nuc_from++) {
            for (char nuc_to{0}; nuc_to < 4; nuc_to++) {
                total += nuc_frequencies[nuc_from] * mutation_rate_matrix[nuc_from][nuc_to];
            }
        }
        return total;
    };

    array<double, 4> sum_mutation_rates() {
        array<double, 4> sum_rates{0};
        for (char nuc_from{0}; nuc_from < 4; nuc_from++) {
            for (char nuc_to{0}; nuc_to < 4; nuc_to++) {
                sum_rates[nuc_from] += mutation_rate_matrix[nuc_from][nuc_to];
            }
        }
        return sum_rates;
    };

    double max_mutation_rate() {
        double maximum = 0.0;
        array<double, 4> sum_rates = sum_mutation_rates();
        for (char nuc{0}; nuc < 4; nuc++) {
            maximum = max(maximum, sum_rates[nuc]);
        }
        return maximum;
    };

    char draw_nuc_to(char nuc_from) {
        double maximum = max_mutation_rate();
        double sum_from = sum_mutation_rates()[nuc_from];
        if (maximum != sum_from) {
            uniform_real_distribution<double> uni_distr(maximum);
            if (uni_distr(Codon::re) > sum_from) {
                return nuc_from;
            }
        }

        discrete_distribution<char> mutation_distr(mutation_rate_matrix[nuc_from].begin(),
                                                   mutation_rate_matrix[nuc_from].end());
        return mutation_distr(Codon::re);
    };

    void check_consistency() {
        assert(haplotype_vector.size() <= population_size);
        assert(!haplotype_vector.empty());

        unsigned nbr_copies{0};
        for (auto &haplotype: haplotype_vector) {
            haplotype.check_consistency(nbr_sites);
            nbr_copies += haplotype.nbr_copies;
        }
        assert(nbr_copies == population_size);
    }

    void forward() {
        mutation();
        selection_and_drift();
        extinction();
        elapsed++;
    }

    void run_forward(unsigned t_max) {
        assert(!haplotype_vector.empty());
        for (unsigned gen{0}; gen <= t_max; gen++) {
            forward();
            if (gen % population_size == 0) {
                fixation();
            }
        }
        if (t_max > 0) {
            fixation();
        }
        check_consistency();
    }

    void mutation() {
        TimeVar t_start = timeNow();
        set<tuple<unsigned, unsigned, unsigned>> coordinates_set{};

        binomial_distribution<unsigned> binomial_distr(population_size * nbr_nucleotides, max_mutation_rate());
        unsigned binomial_draw = binomial_distr(Codon::re);

        if (binomial_draw > 0) {
            vector<unsigned> nucleotide_copies_array(haplotype_vector.size(), 0);
            for (unsigned hap_id{0}; hap_id < haplotype_vector.size(); hap_id++) {
                nucleotide_copies_array[hap_id] = haplotype_vector[hap_id].nbr_copies;
            }
            discrete_distribution<unsigned> haplotype_freq_distr(nucleotide_copies_array.begin(),
                                                                 nucleotide_copies_array.end());
            unsigned nbr_draws{0};
            while (nbr_draws < binomial_draw) {
                unsigned haplotype_draw = haplotype_freq_distr(Codon::re);

                uniform_int_distribution<unsigned> copy_and_site_distr(0, haplotype_vector[haplotype_draw].nbr_copies *
                                                                          nbr_nucleotides - 1);
                unsigned copy_and_site_draw = copy_and_site_distr(Codon::re);
                unsigned nuc_site = copy_and_site_draw % nbr_nucleotides;
                unsigned copy = copy_and_site_draw / nbr_nucleotides;

                auto coordinate = make_tuple(haplotype_draw, copy, nuc_site);
                if (coordinates_set.count(coordinate) == 0) {
                    coordinates_set.insert(coordinate);
                    nbr_draws++;
                }
            }
        }

        auto iter = coordinates_set.begin();
        while (iter != coordinates_set.end()) {
            unsigned hap_id = get<0>(*iter);
            unsigned copy_id = get<1>(*iter);
            bool at_least_one_mutation{false};
            Haplotype haplotype{};
            while (true) {
                unsigned site = get<2>(*iter);
                unsigned codon_site = site / 3;
                auto nuc_position = static_cast<char>(site % 3);

                char codon_from{-1};
                bool has_change = haplotype_vector[hap_id].change_map.count(codon_site) == 1;
                if (has_change) {
                    codon_from = haplotype_vector[hap_id].change_map[codon_site].codon_to;
                } else {
                    codon_from = codon_seq[codon_site];
                }

                array<char, 3> triplet_nuc = Codon::codon_to_triplet_array[codon_from];
                char nuc_from = triplet_nuc[nuc_position];
                char nuc_to = draw_nuc_to(nuc_from);

                if (nuc_from != nuc_to) {
                    triplet_nuc[nuc_position] = nuc_to;
                    char codon_to = Codon::triplet_to_codon(triplet_nuc[0], triplet_nuc[1], triplet_nuc[2]);
                    if (Codon::codon_to_aa_array[codon_to] != 20) {
                        Change change{mutations_vector.size() + 1, codon_from, codon_to};
                        if (not at_least_one_mutation) {
                            haplotype = haplotype_vector[hap_id];
                        }
                        haplotype.change_map[codon_site] = change;
                        haplotype.fitness -= aa_fitness_profiles[codon_site][Codon::codon_to_aa_array[codon_from]];
                        haplotype.fitness += aa_fitness_profiles[codon_site][Codon::codon_to_aa_array[codon_to]];
                        mutations_vector.emplace_back(change);
                        nbr_mutations++;
                        at_least_one_mutation = true;
                    }
                }

                iter++;
                if (iter == coordinates_set.end() or (get<0>(*iter) != hap_id) or (get<1>(*iter) != copy_id)) {
                    if (at_least_one_mutation) {
                        haplotype_vector[hap_id].nbr_copies--;
                        haplotype.nbr_copies = 1;
                        haplotype_vector.emplace_back(haplotype);
                    }
                    break;
                }
            }

        }
        time.mutation += duration(timeNow() - t_start);
    }

    void selection_and_drift() {
        TimeVar t_start = timeNow();

        vector<double> hap_fit_array(haplotype_vector.size(), 0);
        for (unsigned hap_id{0}; hap_id < haplotype_vector.size(); hap_id++) {
            double fitness = max(0.0, 1.0 + haplotype_vector[hap_id].fitness / population_size);
            hap_fit_array[hap_id] = haplotype_vector[hap_id].nbr_copies * fitness;
            haplotype_vector[hap_id].nbr_copies = 0;
        }
        discrete_distribution<unsigned> haplotype_freq_distr(hap_fit_array.begin(), hap_fit_array.end());
        for (unsigned indiv{0}; indiv < population_size; indiv++) {
            haplotype_vector[haplotype_freq_distr(Codon::re)].nbr_copies++;
        }
        total_generation++;

        time.selection += duration(timeNow() - t_start);
    }

    void extinction() {
        TimeVar t_start = timeNow();

        // remove haplotypes with 0 copies
        sort(haplotype_vector.begin(), haplotype_vector.end(), Haplotype::GreaterThan);
        auto low_bound = lower_bound(haplotype_vector.begin(), haplotype_vector.end(), 0, Haplotype::GreaterThan);
        haplotype_vector.erase(low_bound, haplotype_vector.end());

        time.extinction += duration(timeNow() - t_start);
    }

    void fixation() {
        TimeVar t_start = timeNow();
        set<unsigned> site_set;
        transform(haplotype_vector[0].change_map.begin(), haplotype_vector[0].change_map.end(),
                  inserter(site_set, site_set.end()),
                  [](auto pair) { return pair.first; });
        unsigned hap_id{1};
        while (!site_set.empty() and hap_id < haplotype_vector.size()) {
            set<unsigned> hap_site_set;
            transform(haplotype_vector[hap_id].change_map.begin(), haplotype_vector[hap_id].change_map.end(),
                      inserter(hap_site_set, hap_site_set.end()),
                      [](auto pair) { return pair.first; });
            set<unsigned> intersect;
            set_intersection(site_set.begin(), site_set.end(), hap_site_set.begin(), hap_site_set.end(),
                             inserter(intersect, intersect.begin()));
            site_set = intersect;
            hap_id++;
        }
        for (auto site: site_set) {
            set<Change> set_changes;
            for (auto &haplotype: haplotype_vector) {
                set_changes.insert(haplotype.change_map[site]);
            }

            if (set_changes.size() == 1) {
                for (auto &haplotype: haplotype_vector) {
                    haplotype.change_map.erase(site);
                }
                codon_seq[site] = (*set_changes.begin()).codon_to;
                fixations_vector.emplace_back(*set_changes.begin());
                nbr_fixations++;
            }
        }

        time.fixation += duration(timeNow() - t_start);
    }

    void output_vcf(string &output_filename, string &node_name, unsigned gen) {
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
        vcf_file << "##sequenceSize=" << nbr_nucleotides << endl;
        vcf_file << "##populationSize=" << population_size << endl;
        vcf_file << "##numberHaplotypes=" << haplotype_vector.size() << endl;
        vcf_file << "##reference=" << get_dna_str() << endl;
        vcf_file << "##FILTER=<ID=s50,Description=\"The alternative is the major allele\">" << endl;
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

        for (unsigned site{0}; site < nbr_sites; site++) {
            map<tuple<char, char>, unsigned> codon_from_to_copy{};
            for (auto &haplotype: haplotype_vector) {
                if (haplotype.change_map.count(site) == 1) {
                    char codon_from = haplotype.change_map[site].codon_from;
                    char codon_to = haplotype.change_map[site].codon_to;
                    if (codon_to != codon_from) {
                        codon_from_to_copy[make_tuple(codon_from, codon_to)] += haplotype.nbr_copies;
                    }
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
                        if (kv_tuple_copy.second > population_size / 2) {
                            line += "s50";
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
                        for (auto &haplotype: haplotype_vector) {
                            for (unsigned copy{0}; copy < haplotype.nbr_copies; copy++) {
                                line += "\t";
                                if (haplotype.change_map.count(site) == 1) {
                                    line += Codon::codon_to_nuc(haplotype.change_map[site].codon_to, position);
                                } else {
                                    line += Codon::codon_to_nuc(codon_seq[site], position);
                                }
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
        double mean_fitness = accumulate(haplotype_vector.begin(), haplotype_vector.end(), 0.0,
                                         [](double acc, Haplotype &h) { return acc + h.fitness; }) /
                              haplotype_vector.size();
        // .ali format
        ofstream tsv_file;
        tsv_file.open(output_filename + ".tsv", ios_base::app);
        SummaryStatistics stats = compute_summary_stats();
        stats_vector.push_back(stats);
        tsv_file << node_name << "\t" << elapsed << "\tdN\t" << stats.dn << endl;
        tsv_file << node_name << "\t" << elapsed << "\tdS\t" << stats.ds << endl;
        tsv_file << node_name << "\t" << elapsed << "\tomega\t" << stats.omega << endl;
        tsv_file << node_name << "\t" << elapsed << "\tFixedNbr\t" << fixations_vector.size() << endl;

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

        time.exportation += duration(timeNow() - t_start);
    }

    SummaryStatistics compute_summary_stats() {
        SummaryStatistics stats;
        Change_table table = compute_change_table();
        stats.dn = dn_from_table(table);
        stats.ds = ds_from_table(table);
        stats.omega = stats.dn / stats.ds;

        double copies = 0.0;
        for (auto const &haplotype: haplotype_vector) {
            auto hap_table = table;
            for (auto const &kv_change: haplotype.change_map) {
                if (is_synonymous(kv_change.second)) {
                    hap_table.syn_fix++;
                } else {
                    hap_table.non_syn_fix++;
                }
            }
            auto hap_dn = dn_from_table(hap_table);
            auto hap_ds = ds_from_table(hap_table);
            double hap_omega = hap_dn / hap_ds;
            if (isfinite(hap_omega)) {
                copies += haplotype.nbr_copies;
                stats.dn_sample += hap_dn * haplotype.nbr_copies;
                stats.ds_sample += hap_ds * haplotype.nbr_copies;
                stats.omega_sample += haplotype.nbr_copies * hap_omega;
            }
        }
        stats.dn_sample /= copies;
        stats.ds_sample /= copies;
        stats.omega_sample /= copies;

        stats.omega_predicted = predicted_omega();
        return stats;
    };

    Change_table compute_change_table(bool rmv_mutations = false) {
        long syn_mut{0}, non_syn_mut{0};
        if (rmv_mutations) {
            vector<Change> mutation_absorbed(mutations_vector.size());

            set<Change> polymorphic_set;
            for (auto const &haplotype: haplotype_vector) {
                for (auto const &kv_change: haplotype.change_map) {
                    polymorphic_set.insert(kv_change.second);
                }
            }

            assert(is_sorted(polymorphic_set.begin(), polymorphic_set.end()));
            if (!mutations_vector.empty()) {
                assert(is_sorted(mutations_vector.begin(), mutations_vector.end()));
                for (auto &change: polymorphic_set) {
                    assert(lower_bound(mutations_vector.begin(), mutations_vector.end(), change) !=
                           mutations_vector.end());
                }
            }

            auto diff_it = set_difference(mutations_vector.begin(), mutations_vector.end(),
                                          polymorphic_set.begin(), polymorphic_set.end(),
                                          mutation_absorbed.begin());
            mutation_absorbed.resize(diff_it - mutation_absorbed.begin());
            assert(mutation_absorbed.size() <= mutations_vector.size());
            assert(mutation_absorbed.size() >= fixations_vector.size());
            syn_mut = count_if(mutation_absorbed.begin(), mutation_absorbed.end(), is_synonymous);
            non_syn_mut = mutation_absorbed.size() - syn_mut;
        } else {
            syn_mut = count_if(mutations_vector.begin(), mutations_vector.end(), is_synonymous);
            non_syn_mut = mutations_vector.size() - syn_mut;
        }
        long syn_fix = count_if(fixations_vector.begin(), fixations_vector.end(), is_synonymous);
        long non_syn_fix = fixations_vector.size() - syn_fix;

        return Change_table{non_syn_mut, syn_mut, non_syn_fix, syn_fix};
    }

    // Theoretical computation of the predicted omega
    double predicted_omega() {
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

                            dn += codon_freqs[codon_from] * mutation_rate_matrix[n_from][n_to] * rate_fixation;
                            d0 += codon_freqs[codon_from] * mutation_rate_matrix[n_from][n_to];
                        }
                    }
                }

            }
            dn_per_site[site] = dn;
            ds_per_site[site] = d0;
        }
        double flow_non_synonymous = accumulate(dn_per_site.begin(), dn_per_site.end(), 0.0);
        double flow_synonymous = accumulate(ds_per_site.begin(), ds_per_site.end(), 0.0);
        return flow_non_synonymous / flow_synonymous;
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

    static void avg_summary_statistics() {
        cout << total_generation << " generations computed" << endl;
        cout << nbr_fixations << " fixations" << endl;
        cout << nbr_mutations << " mutations" << endl;
        SummaryStatistics avg_stats = average(stats_vector);
        cout << "dN=" << avg_stats.dn << endl;
        cout << "dS=" << avg_stats.ds << endl;
        cout << "omega=" << avg_stats.omega << endl;
        cout << "Random dN=" << avg_stats.dn_sample << endl;
        cout << "Random dS=" << avg_stats.ds_sample << endl;
        cout << "Random omega=" << avg_stats.omega_sample << endl;
        cout << "Predicted omega=" << avg_stats.omega_predicted << endl;
        double total_time =
                time.mutation + time.selection + time.extinction +
                time.fixation + time.exportation;
        cout << setprecision(3) << total_time / 1e9 << "s total time" << endl;
        cout << 100 * time.mutation / total_time << "% of time spent in calculating mutation ("
             << time.mutation / 1e9 << "s)" << endl;
        cout << 100 * time.selection / total_time << "% of time spent in calculating selection ("
             << time.selection / 1e9 << "s)" << endl;
        cout << 100 * time.extinction / total_time << "% of time spent in calculating extinction ("
             << time.extinction / 1e9 << "s)" << endl;
        cout << 100 * time.fixation / total_time << "% of time spent in calculating fixation ("
             << time.fixation / 1e9 << "s)" << endl;
        cout << 100 * time.exportation / total_time << "% of time spent in exportationing vcf ("
             << time.exportation / 1e9 << "s)" << endl;
    }
};

// Initialize static variables
unsigned Population::total_generation = 0;
unsigned long Population::nbr_mutations = 0, Population::nbr_fixations = 0;
TimeElapsed Population::time = {0.0};
vector<SummaryStatistics> Population::stats_vector{};

// Class representing nodes of a tree.
class Node {
public:
    static double length_computed;

private:
    string name;  // The species name of the node.
    double length;  // The length of the branch attached to the node (ascending).
    string newick;  // The newick tree descending from the node.
    vector<Node> children;  // Vector of direct children (first order, no grand-children).

public:
    Population population;  // The population attached to the node.
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
        children.emplace_back(node);
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
    void traverse(string &output_filename, double scale) {
        // Substitutions of the DNA sequence is generated.
        auto gen = static_cast<unsigned>(scale * length / (population.mutation_rate()));
        population.run_forward(gen);

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

            population.output_vcf(output_filename, name, gen);

        } else {
            // If the node is internal, iterate through the direct children.
            for (auto &child : children) {
                child.population.set_parameters(population);
                child.traverse(output_filename, scale);
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

// Initialize static variables
double Node::length_computed = 0.0;

string open_newick(string const &file_name) {
    ifstream input_stream(file_name);
    if (!input_stream) cerr << "Can't open newick file!" << endl;

    string line;
    getline(input_stream, line);

    return line;
}


vector<array<double, 20>> open_preferences(string const &file_name, double const &beta) {
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
                fitness_profil[counter - 3] = beta * log(stod(word));
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

    double beta = 1.0;
    if (args["--beta"]) {
        beta = stod(args["--beta"].asString());
    }
    vector<array<double, 20>> fitness_profiles = open_preferences(preferences_path, beta);

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

    Population population(lambda, mu, fitness_profiles, pop_size, output_path);

    string newick_path{"../data_trees/mammals.newick"};
    if (args["--newick"]) {
        newick_path = args["--newick"].asString();
    }
    bool newick_is_int{true};
    try {
        stoi(newick_path);
    } catch (const invalid_argument &ia) {
        newick_is_int = false;
    }
    if (!newick_is_int) {
        string newick_tree = open_newick(newick_path);
        Node root(newick_tree, population);
        unsigned nbr_leaves = root.nbr_leaves();
        cout << "The tree has " << nbr_leaves << " leaves" << endl;
        cout << "The tree is of length " << root.tot_length() << endl;
        ofstream ali_file;
        ali_file.open(output_path + ".ali");
        ali_file << nbr_leaves << " " << fitness_profiles.size() * 3 << endl;
        ali_file.close();

        ofstream fasta_file;
        fasta_file.open(output_path + ".fasta");
        fasta_file.close();

        ofstream tsv_file;
        tsv_file.open(output_path + ".fasta");
        tsv_file.close();

        root.traverse(output_path, 0.01);
        Population::avg_summary_statistics();
    } else {
        auto generations = static_cast<unsigned>(stoi(newick_path));
        cout << "the tree given is a number and not a newick file (a tree), the simulator will run for "
             << generations << " and output the summary statistics" << endl;
        population.run_forward(generations);
        Population::avg_summary_statistics();
    }
    return 0;
}