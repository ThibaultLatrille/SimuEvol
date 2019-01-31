#include <chrono>
#include "argparse.hpp"
#include "codon.hpp"
#include "io.hpp"
#include "matrices.hpp"
#include "random.hpp"
#include "statistic.hpp"
#include "tree.hpp"

using namespace TCLAP;
using namespace std;

#define duration(a) chrono::duration_cast<std::chrono::nanoseconds>(a).count()
#define timeNow() chrono::high_resolution_clock::now()

typedef chrono::high_resolution_clock::time_point TimeVar;

unsigned long binomial_coefficient(unsigned long n, unsigned long k) {
    assert(0 <= k);
    assert(k <= n);
    unsigned long i;
    unsigned long b;
    if (0 == k || n == k) { return 1; }
    if (k > (n - k)) { k = n - k; }
    if (1 == k) { return n; }
    b = 1;
    for (i = 1; i <= k; ++i) {
        b *= (n - (k - i));
        assert(b > 0); /* Overflow */
        b /= i;
    }
    return b;
}

unsigned pow_int(unsigned x, unsigned p) {
    if (p == 0) return 1;
    if (p == 1) return x;

    unsigned tmp = pow_int(x, p / 2);
    if (p % 2 == 0)
        return tmp * tmp;
    else
        return x * tmp * tmp;
}

string join(vector<unsigned> &v, char sep) {
    return accumulate(v.begin() + 1, v.end(), to_string(v[0]),
        [sep](const string &acc, int b) { return acc + sep + to_string(b); });
};

struct Change {
    unsigned site;
    char codon_from{-1};
    char codon_to{-1};
};

bool operator<(const Change &lhs, const Change &rhs) { return lhs.site < rhs.site; }

bool operator==(const Change &lhs, const Change &rhs) { return lhs.site == rhs.site; }

bool is_synonymous(const Change &change) {
    return Codon::codon_to_aa_array[change.codon_from] == Codon::codon_to_aa_array[change.codon_to];
};

struct Change_table {
    unsigned non_syn_mut;
    unsigned syn_mut;
    unsigned non_syn_fix;
    unsigned syn_fix;
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
    double mutation{0.0};
    double selection{0.0};
    double extinction{0.0};
    double fixation{0.0};
    double exportation{0.0};
    double correlation{0.0};
    double normal{0.0};
    double matrix{0.0};
    double exp{0.0};
};

TimeElapsed time_elapsed;
Trace trace;

struct SummaryStatistics {
    double dn{0}, ds{0}, omega{0}, omega_predicted{0};
    double dn_sample{0}, ds_sample{0}, omega_sample{0};
    double pn_sample{0}, ps_sample{0}, pnps_sample{0};
    double pin_sample{0}, pis_sample{0}, pinpis_sample{0};
    double at_pct{0}, mean_fitness{0};
};

// Function for avg of a vector ignoring NaN
double average(vector<SummaryStatistics> const &vector_stats,
    function<double(SummaryStatistics)> const &func) {
    double total{0};
    unsigned elements_non_nan{0};
    for (auto &stat : vector_stats) {
        double x = func(stat);
        if (isfinite(x)) {
            total += x;
            elements_non_nan++;
        }
    }
    return total / elements_non_nan;
}

void average(vector<SummaryStatistics> const &vector_stats) {
    cout << "dN=" << average(vector_stats, [](SummaryStatistics stat) { return stat.dn; }) << endl;
    cout << "dS=" << average(vector_stats, [](SummaryStatistics stat) { return stat.ds; }) << endl;
    cout << "omega=" << average(vector_stats, [](SummaryStatistics stat) { return stat.omega; })
         << endl;
    cout << "Sample omega="
         << average(vector_stats, [](SummaryStatistics stat) { return stat.omega_sample; }) << endl;
    cout << "Pred omega=" << average(vector_stats, [](SummaryStatistics stat) {
        return stat.omega_predicted;
    }) << endl;
    cout << "Sample piN/piS="
         << average(vector_stats, [](SummaryStatistics stat) { return stat.pnps_sample; }) << endl;
    cout << "Sample pN/pS=" << average(vector_stats, [](SummaryStatistics stat) {
        return stat.pinpis_sample;
    }) << endl;
    cout << "%AT=" << average(vector_stats, [](SummaryStatistics stat) { return stat.at_pct; })
         << endl;
}

class Haplotype {
  public:
    // The nbr_copies of sites in the sequence (each position is a codon, thus the DNA sequence is 3
    // times greater).
    unsigned nbr_copies{0};
    double fitness{0.0};
    set<Change> set_change;

    // Constructor of Reference_seq.
    // size: the size of the DNA sequence.
    Haplotype(unsigned const nbr_copies, double const fitness)
        : nbr_copies{nbr_copies}, fitness{fitness} {};

    explicit Haplotype() = default;

    void check_consistency(unsigned nbr_sites) const {
        for (auto &change : set_change) {
            assert(Codon::codon_to_aa_array[change.codon_to] != 20);
            assert(Codon::codon_to_aa_array[change.codon_from] != 20);
        }
    }

    static struct {
        bool operator()(const Haplotype &left, const Haplotype &right) {
            return left.nbr_copies > right.nbr_copies;
        }

        bool operator()(const Haplotype &left, float right) { return left.nbr_copies > right; }

        bool operator()(float left, const Haplotype &right) { return left > right.nbr_copies; }
    } GreaterThan;
};

// Class representing a genetically linked sequences (blocks are unlinked between them)
class Block {
  public:
    // Block
    unsigned &population_size;
    unsigned position;

    // Selection
    vector<array<double, 20>> aa_fitness_profiles;

    // Reference sequence
    unsigned nbr_sites;
    unsigned nbr_nucleotides;
    vector<char> codon_seq;

    // Haplotypes
    vector<Haplotype> haplotype_vector;

    // Keep track of mutations and fixations
    Change_table change_table{0, 0, 0, 0};

    // Statics variables (shared by all instances)
    static unsigned nbr_mutations, nbr_fixations;

    // Constructor
    explicit Block(vector<array<double, 20>> const &fitness_profiles, const unsigned &position,
        unsigned &population_size, NucleotideRateMatrix const &nuc_matrix)
        : population_size{population_size},
          position{position},
          aa_fitness_profiles{fitness_profiles},
          nbr_sites{unsigned(fitness_profiles.size())},
          nbr_nucleotides{unsigned(3 * fitness_profiles.size())},
          codon_seq(fitness_profiles.size(), 0),
          haplotype_vector{} {
        // Draw codon from codon frequencies
        for (unsigned site{0}; site < nbr_sites; site++) {
            array<double, 64> codon_freqs =
                codon_frequencies(aa_fitness_profiles[site], nuc_matrix, 4 * population_size);
            discrete_distribution<char> freq_codon_distr(codon_freqs.begin(), codon_freqs.end());
            codon_seq[site] = freq_codon_distr(generator);
        }

        haplotype_vector.emplace_back(Haplotype(2 * population_size, 0.0));
        check_consistency();
        assert(!haplotype_vector.empty());
    }

    void check_consistency() const {
        assert(haplotype_vector.size() <= 2 * population_size);
        assert(!haplotype_vector.empty());

        unsigned nbr_copies{0};
        for (auto &haplotype : haplotype_vector) {
            haplotype.check_consistency(nbr_sites);
            nbr_copies += haplotype.nbr_copies;
        }
        assert(nbr_copies == 2 * population_size);
    }

    void forward(NucleotideRateMatrix const &nuc_matrix) {
        mutation(nuc_matrix);
        selection_and_drift();
        extinction();
        fixation();
    }

    discrete_distribution<unsigned> haplotype_freq_distr() const {
        vector<unsigned> nucleotide_copies_array(haplotype_vector.size(), 0);
        for (unsigned i_hap{0}; i_hap < haplotype_vector.size(); i_hap++) {
            nucleotide_copies_array[i_hap] = haplotype_vector[i_hap].nbr_copies;
        }
        discrete_distribution<unsigned> haplotype_distr(
            nucleotide_copies_array.begin(), nucleotide_copies_array.end());
        return haplotype_distr;
    }

    void mutation(NucleotideRateMatrix const &p) {
        TimeVar t_start = timeNow();

        binomial_distribution<unsigned> binomial_distr(
            2 * population_size * nbr_nucleotides, p.max_sum_mutation_rates);
        unsigned binomial_draw = binomial_distr(generator);

        if (binomial_draw > 0) {
            set<tuple<unsigned, unsigned, unsigned>> coordinates_set{};
            unsigned nbr_draws{0};

            discrete_distribution<unsigned> haplotype_distr{};
            if (haplotype_vector.size() > 1) { haplotype_distr = haplotype_freq_distr(); }

            while (nbr_draws < binomial_draw) {
                unsigned haplotype_draw = 0;

                if (haplotype_vector.size() > 1) { haplotype_draw = haplotype_distr(generator); }

                uniform_int_distribution<unsigned> copy_and_site_distr(
                    0, haplotype_vector[haplotype_draw].nbr_copies * nbr_nucleotides - 1);
                unsigned copy_and_site_draw = copy_and_site_distr(generator);
                unsigned nuc_site = copy_and_site_draw % nbr_nucleotides;
                unsigned copy = copy_and_site_draw / nbr_nucleotides;

                auto coordinate = make_tuple(haplotype_draw, copy, nuc_site);
                if (coordinates_set.count(coordinate) == 0) {
                    coordinates_set.insert(coordinate);
                    nbr_draws++;
                }
            }

            auto iter = coordinates_set.begin();
            while (iter != coordinates_set.end()) {
                unsigned i_hap = get<0>(*iter);
                unsigned copy_id = get<1>(*iter);
                bool at_least_one_mutation{false};
                Haplotype haplotype{};
                while (true) {
                    unsigned site = get<2>(*iter);
                    unsigned codon_site = site / 3;
                    auto nuc_position = static_cast<char>(site % 3);

                    char codon_from{0};
                    auto it = haplotype_vector[i_hap].set_change.find(Change{codon_site});
                    if (it != haplotype_vector[i_hap].set_change.end()) {
                        assert(it->site == codon_site);
                        codon_from = it->codon_to;
                    } else {
                        codon_from = codon_seq[codon_site];
                    }

                    array<char, 3> triplet_nuc = Codon::codon_to_triplet_array[codon_from];
                    char nuc_from = triplet_nuc[nuc_position];

                    bool draw_mutation = true;

                    if (p.max_sum_mutation_rates != p.sum_mutation_rates(nuc_from)) {
                        double sum_unif = p.max_real_distr(generator);
                        if (sum_unif > p.sum_mutation_rates(nuc_from)) { draw_mutation = false; }
                    }

                    if (draw_mutation) {
                        triplet_nuc[nuc_position] = p.mutation_distr[nuc_from](generator);
                        char codon_to =
                            Codon::triplet_to_codon(triplet_nuc[0], triplet_nuc[1], triplet_nuc[2]);
                        if (Codon::codon_to_aa_array[codon_to] != 20) {
                            Change change{codon_site + position, codon_from, codon_to};
                            if (not at_least_one_mutation) { haplotype = haplotype_vector[i_hap]; }
                            haplotype.set_change.insert(change);
                            haplotype.fitness -=
                                aa_fitness_profiles[codon_site]
                                                   [Codon::codon_to_aa_array[codon_from]];
                            haplotype.fitness +=
                                aa_fitness_profiles[codon_site][Codon::codon_to_aa_array[codon_to]];
                            if (is_synonymous(change)) {
                                change_table.syn_mut++;
                            } else {
                                change_table.non_syn_mut++;
                            }
                            nbr_mutations++;
                            at_least_one_mutation = true;
                        }
                    }

                    iter++;
                    if (iter == coordinates_set.end() or (get<0>(*iter) != i_hap) or
                        (get<1>(*iter) != copy_id)) {
                        if (at_least_one_mutation) {
                            haplotype_vector[i_hap].nbr_copies--;
                            haplotype.nbr_copies = 1;
                            haplotype_vector.push_back(haplotype);
                        }
                        break;
                    }
                }
            }
        }
        time_elapsed.mutation += duration(timeNow() - t_start);
    }

    void selection_and_drift() {
        TimeVar t_start = timeNow();
        if (haplotype_vector.size() > 1) {
            vector<double> fitness_vector(haplotype_vector.size(), 0);

            for (unsigned i_hap{0}; i_hap < haplotype_vector.size(); i_hap++) {
                fitness_vector[i_hap] =
                    (1 + haplotype_vector[i_hap].fitness) * haplotype_vector[i_hap].nbr_copies;
                haplotype_vector[i_hap].nbr_copies = 0;
            }

            double fit_tot = accumulate(fitness_vector.begin(), fitness_vector.end(), 0.0);
            unsigned children_tot = 2 * population_size;

            for (unsigned i_hap{0}; i_hap < haplotype_vector.size() - 1; i_hap++) {
                binomial_distribution<unsigned> bin(children_tot, fitness_vector[i_hap] / fit_tot);
                unsigned children = bin(generator);
                haplotype_vector[i_hap].nbr_copies = children;
                children_tot -= children;
                fit_tot -= fitness_vector[i_hap];
            }
            assert(children_tot <= 2 * population_size);
            assert(children_tot >= 0);
            haplotype_vector[haplotype_vector.size() - 1].nbr_copies = children_tot;
        }
        time_elapsed.selection += duration(timeNow() - t_start);
    }

    void extinction() {
        TimeVar t_start = timeNow();

        // remove haplotypes with 0 copies
        if (haplotype_vector.size() > 1) {
            if (!is_sorted(
                    haplotype_vector.begin(), haplotype_vector.end(), Haplotype::GreaterThan)) {
                sort(haplotype_vector.begin(), haplotype_vector.end(), Haplotype::GreaterThan);
            }
            if (haplotype_vector.rbegin()->nbr_copies == 0) {
                auto low_bound = lower_bound(
                    haplotype_vector.begin(), haplotype_vector.end(), 0, Haplotype::GreaterThan);
                if (low_bound != haplotype_vector.end()) {
                    haplotype_vector.erase(low_bound, haplotype_vector.end());
                }
            }
        }

        time_elapsed.extinction += duration(timeNow() - t_start);
    }

    void fixation() {
        TimeVar t_start = timeNow();
        if (nbr_sites == 1 and haplotype_vector.size() == 1 and
            haplotype_vector.begin()->set_change.size() == 1) {
            Change change = *(haplotype_vector.begin()->set_change.begin());
            haplotype_vector.begin()->set_change.clear();
            codon_seq[0] = change.codon_to;
            if (is_synonymous(change)) {
                change_table.syn_fix++;
            } else {
                change_table.non_syn_fix++;
            }
            nbr_fixations++;
        } else if (nbr_sites > 1 or (nbr_sites == 1 and haplotype_vector.size() > 1)) {
            set<Change> current_change_set = haplotype_vector[0].set_change;
            unsigned i_hap{1};
            while (!current_change_set.empty() and i_hap < haplotype_vector.size()) {
                set<Change> intersect;
                set_intersection(current_change_set.begin(), current_change_set.end(),
                    haplotype_vector[i_hap].set_change.begin(),
                    haplotype_vector[i_hap].set_change.end(),
                    inserter(intersect, intersect.begin()));
                current_change_set = intersect;
                i_hap++;
            }
            for (auto change : current_change_set) {
                set<Change> set_changes;
                for (auto &haplotype : haplotype_vector) {
                    set_changes.insert(*haplotype.set_change.find(change));
                }

                if (set_changes.size() == 1) {
                    for (auto &haplotype : haplotype_vector) {
                        haplotype.set_change.erase(*set_changes.begin());
                    }
                    assert(change.site - position >= 0);
                    assert(change.site - position < codon_seq.size());
                    codon_seq[change.site - position] = set_changes.begin()->codon_to;
                    if (is_synonymous(*set_changes.begin())) {
                        change_table.syn_fix++;
                    } else {
                        change_table.non_syn_fix++;
                    }
                    nbr_fixations++;
                }
            }
        }
        time_elapsed.fixation += duration(timeNow() - t_start);
    }
};

class Population {
  public:
    // TimeElapsed
    double time_from_root{0};

    // Blocks
    vector<Block> blocks;
    unsigned sample_size{0};

    LogMultivariate log_multivariate;
    unsigned population_size;
    double generation_time;
    NucleotideRateMatrix nuc_matrix;

    CorrelationMatrix const &cor_matrix;

    // Statics variables (shared by all instances)
    static vector<SummaryStatistics> stats_vector;

    explicit Population() = default;

    Population(vector<array<double, 20>> const &fitness_profiles, unsigned sample_size,
        LogMultivariate &log_multi, bool linked_sites,
        NucleotideRateMatrix const &nucleotide_matrix, CorrelationMatrix const &cor_matrix)
        : blocks{},
          sample_size{sample_size},
          log_multivariate{log_multi},
          population_size{log_multivariate.population_size()},
          generation_time{log_multivariate.generation_time()},
          nuc_matrix{nucleotide_matrix},
          cor_matrix{cor_matrix} {
        if (linked_sites) {
            blocks.emplace_back(Block(fitness_profiles, 0, population_size, nuc_matrix));
        } else {
            blocks.reserve(fitness_profiles.size());
            for (unsigned site{0}; site < fitness_profiles.size(); site++) {
                vector<array<double, 20>> site_fitness_profile = {fitness_profiles[site]};
                blocks.emplace_back(Block(site_fitness_profile, site, population_size, nuc_matrix));
            }
        }
        cout << blocks.size() << " blocks created" << endl;
    }

    void check_consistency() {
        for (auto const &block : blocks) { block.check_consistency(); }
    };

    void run_forward(double t_max) {
        double time = 0.0;
        TimeVar t_start = timeNow();
        Eigen::SelfAdjointEigenSolver<EMatrix> eigen_solver(cor_matrix);
        EMatrix transform =
            eigen_solver.eigenvectors() * eigen_solver.eigenvalues().cwiseSqrt().asDiagonal();
        EVector sampled_vector = EVector::Zero(cor_matrix.dimensions);
        time_elapsed.correlation += duration(timeNow() - t_start);

        while (time < t_max) {
            t_start = timeNow();
            for (int dim = 0; dim < cor_matrix.dimensions; dim++) {
                sampled_vector(dim) = generation_time * normal_distrib(generator);
            }
            time_elapsed.normal += duration(timeNow() - t_start);

            t_start = timeNow();
            log_multivariate += transform * sampled_vector;
            time_elapsed.matrix += duration(timeNow() - t_start);

            t_start = timeNow();
            population_size = log_multivariate.population_size();
            generation_time = log_multivariate.generation_time();
            nuc_matrix.set_mutation_rate(log_multivariate.mu());
            time_elapsed.exp += duration(timeNow() - t_start);

            if (population_size < sample_size) {
                cerr << "The population size became lower than sample size (reflecting back)."
                     << endl;
                population_size = 2 * sample_size - population_size;
                log_multivariate.set_population_size(population_size);
            };

            for (auto &block : blocks) {
                assert(!block.haplotype_vector.empty());
                block.forward(nuc_matrix);
            }
            check_consistency();
            time += generation_time;
            time_from_root += generation_time;
        }
        check_consistency();
    }

    void burn_in(unsigned burn_in_length) {
        cout << "Burn-in for " << burn_in_length << " generations." << endl;
        for (unsigned gen{1}; gen <= burn_in_length; gen++) {
            for (auto &block : blocks) {
                assert(!block.haplotype_vector.empty());
                block.forward(nuc_matrix);
            }
            check_consistency();
        }
        cout << "Burn-in completed." << endl;
    }

    void theoretical_dnds(SummaryStatistics &stats) const {
        double sub_flow{0.}, mut_flow{0.};

        for (auto const &block : blocks) {
            double dn{0}, d0{0};
            tie(dn, d0) =
                predicted_dn_d0(block.aa_fitness_profiles, nuc_matrix, 4 * population_size);
            sub_flow += dn;
            mut_flow += d0;
        }
        stats.omega_predicted = sub_flow / mut_flow;
    };

    void output_vcf(string &output_filename, string &node_name) const {
        TimeVar t_start = timeNow();

        SummaryStatistics stats;
        // Predicted dN/dS and %AT
        theoretical_dnds(stats);

        double at_sites_obs{0};
        unsigned nbr_nucleotides{0};

        Change_table table{0, 0, 0, 0};
        unsigned sum_population_size{0};

        for (auto &block : blocks) {
            table.non_syn_mut += block.change_table.non_syn_mut;
            table.syn_mut += block.change_table.syn_mut;
            table.non_syn_fix += block.change_table.non_syn_fix;
            table.syn_fix += block.change_table.syn_fix;
            sum_population_size += block.population_size;
            nbr_nucleotides += block.nbr_nucleotides;
        }

        stats.dn += 2 * sum_population_size * dn_from_table(table) / blocks.size();
        stats.ds += 2 * sum_population_size * ds_from_table(table) / blocks.size();
        stats.omega = stats.dn / stats.ds;

        // VCF file on the sample
        vector<unsigned> sfs_non_syn(2 * sample_size + 1, 0), sfs_syn(2 * sample_size + 1, 0);
        unsigned non_syn_nbr = 0, syn_nbr = 0;
        unsigned complex_sites{0};
        assert(sfs_non_syn.size() == 2 * sample_size + 1);
        assert(sum(sfs_non_syn) == 0);
        assert(sfs_syn.size() == 2 * sample_size + 1);
        assert(sum(sfs_syn) == 0);

        string out;
        out += "##fileformat=VCFv4.0";
        out += "\n##source=SimuPoly";
        out += "\n##nodeName=" + node_name;
        out += "\n##sequenceSize=" + to_string(nbr_nucleotides);
        out += "\n##ploidyLevel=diploid";
        out += "\n##numberIndividuals=" + to_string(sample_size);
        out += "\n##numberGenotypes=" + to_string(2 * sample_size);
        out += "\n##reference=" + get_dna_str();
        out += "\n##FILTER=<ID=s50,Description=\"The alternative is the major allele\">";
        out += "\n##FILTER=<ID=s100,Description=\"The reference has not been sampled\">";
        out += "\n##INFO=<ID=REFCODON,Number=1,Type=String,Description=\"Codon of the reference\">";
        out +=
            "\n##INFO=<ID=ALTCODON,Number=1,Type=String,Description=\"Codon of the alternative\">";
        out +=
            "\n##INFO=<ID=REFAA,Number=1,Type=Character,Description=\"Amino-acid of the "
            "reference\">";
        out +=
            "\n##INFO=<ID=ALTAA,Number=1,Type=Character,Description=\"Amino-acid of the "
            "alternative\">";
        out +=
            "\n##INFO=<ID=POSITION,Number=1,Type=Integer,Description=\"Mutated codon position\">";
        out +=
            "\n##INFO=<ID=ALTCOUNT,Number=1,Type=Integer,Description=\"Number of alternative "
            "copy\">";
        out += "\n##INFO=<ID=SYN,Number=1,Type=boolean,Description=\"Is a synonymous mutation\">";
        out += "\n##FORMAT=<ID=GT,Number=1,Type=String,Description=\"Genotype\">";
        out += "\n#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT";

        for (unsigned indiv{1}; indiv <= sample_size; indiv++) {
            out += "\tId";
            out += to_string(indiv);
        }

        for (auto const &block : blocks) {
            // Draw the sample of individuals
            vector<unsigned> haplotypes_sample(2 * block.population_size, 0);
            unsigned sum_copies = 0;
            for (unsigned i_hap{0}; i_hap < block.haplotype_vector.size(); i_hap++) {
                for (unsigned copy_id{0}; copy_id < block.haplotype_vector[i_hap].nbr_copies;
                     copy_id++) {
                    assert(sum_copies + copy_id < 2 * block.population_size);
                    haplotypes_sample[sum_copies + copy_id] = i_hap;
                }
                sum_copies += block.haplotype_vector[i_hap].nbr_copies;
            }
            shuffle(haplotypes_sample.begin(), haplotypes_sample.end(), generator);
            haplotypes_sample.resize(2 * sample_size);

            // dN/dS computed using one individual of the sample
            uniform_int_distribution<unsigned> chosen_distr(0, 2 * sample_size - 1);
            unsigned chosen = chosen_distr(generator);
            for (auto const &change :
                block.haplotype_vector[haplotypes_sample[chosen]].set_change) {
                if (is_synonymous(change)) {
                    table.syn_fix++;
                } else {
                    table.non_syn_fix++;
                }
            }

            for (unsigned site{0}; site < block.nbr_sites; site++) {
                map<tuple<char, char>, unsigned> codon_from_to_copy{};
                for (auto const &i_hap : haplotypes_sample) {
                    char codon_to = block.codon_seq[site];

                    auto it = block.haplotype_vector[i_hap].set_change.find(
                        Change{block.position + site});
                    if (it != block.haplotype_vector[i_hap].set_change.end()) {
                        assert(it->site == site + block.position);
                        char codon_from = it->codon_from;
                        codon_to = it->codon_to;
                        if (codon_to != codon_from) {
                            codon_from_to_copy[make_tuple(codon_from, codon_to)]++;
                        }
                    }

                    array<char, 3> triplet = Codon::codon_to_triplet_array[codon_to];
                    for (char position{0}; position < 3; position++) {
                        if (triplet[position] == 0 or triplet[position] == 3) { at_sites_obs++; }
                    }
                }

                if (codon_from_to_copy.size() == 1) {
                    char codon_from = get<0>(codon_from_to_copy.begin()->first);
                    char codon_to = get<1>(codon_from_to_copy.begin()->first);
                    if (codon_to != codon_from) {
                        unsigned alt_freq = codon_from_to_copy.begin()->second;
                        char position{0};
                        while (position < 3) {
                            if (Codon::codon_to_nuc(codon_from, position) !=
                                Codon::codon_to_nuc(codon_to, position)) {
                                break;
                            } else {
                                position++;
                            }
                        }
                        assert(position != 3);
                        string line{"\n"};
                        line += ".\t";
                        line += to_string(3 * (block.position + site) + position);
                        line += "\t.\t";
                        line += Codon::codon_to_nuc(codon_from, position);
                        line += "\t";
                        line += Codon::codon_to_nuc(codon_to, position);
                        line += "\t100\t";
                        if (alt_freq == 2 * sample_size) {
                            line += "s100";
                        } else if (alt_freq > sample_size) {
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
                        line += to_string(alt_freq);
                        line += ";SYN=";

                        if (Codon::codon_to_aa_array[codon_from] ==
                            Codon::codon_to_aa_array[codon_to]) {
                            sfs_syn[alt_freq]++;
                            syn_nbr++;
                            line += "TRUE";
                        } else {
                            sfs_non_syn[alt_freq]++;
                            non_syn_nbr++;
                            line += "FALSE";
                        }

                        line += "\tGT";
                        for (unsigned indiv{0}; indiv < sample_size; indiv++) {
                            line += "\t";
                            for (unsigned ploidy{0}; ploidy < 2; ploidy++) {
                                unsigned i_hap = haplotypes_sample[indiv * 2 + ploidy];
                                auto it = block.haplotype_vector[i_hap].set_change.find(
                                    Change{block.position + site});
                                char nuc{0};
                                if (it != block.haplotype_vector[i_hap].set_change.end()) {
                                    assert(it->site == block.position + site);
                                    nuc = Codon::codon_to_nuc(it->codon_to, position);
                                } else {
                                    nuc = Codon::codon_to_nuc(block.codon_seq[site], position);
                                }
                                line += nuc;
                                if (ploidy == 0) { line += "|"; }
                            }
                        }
                        out += line;
                    }
                } else if (codon_from_to_copy.size() > 1) {
                    complex_sites++;
                }
            }

            // piN/piS computed on the sample
            for (unsigned i{0}; i < 2 * sample_size; i++) {
                for (unsigned j{i + 1}; j < 2 * sample_size; j++) {
                    if (i != j) {
                        unsigned hap_i = haplotypes_sample[i];
                        unsigned hap_j = haplotypes_sample[j];
                        auto it_first = block.haplotype_vector[hap_i].set_change.begin();
                        auto end_first = block.haplotype_vector[hap_i].set_change.end();
                        auto it_second = block.haplotype_vector[hap_j].set_change.begin();
                        auto end_second = block.haplotype_vector[hap_j].set_change.end();
                        while (it_first != end_first and it_second != end_second) {
                            Change diff{};
                            if (it_second == end_second or (*it_first) < (*it_second)) {
                                diff = *it_first;
                                it_first++;
                            } else if (it_first == end_first or (*it_second) < (*it_first)) {
                                diff = *it_second;
                                it_second++;
                            } else if ((*it_first) == (*it_second)) {
                                diff.codon_to = it_first->codon_to;
                                diff.codon_from = it_second->codon_from;
                                it_first++;
                                it_second++;
                            }
                            if (is_synonymous(diff)) {
                                stats.pis_sample++;
                            } else {
                                stats.pin_sample++;
                            }
                        }
                    }
                }
            }

            stats.mean_fitness +=
                accumulate(block.haplotype_vector.begin(), block.haplotype_vector.end(), 0.0,
                    [](double acc, Haplotype const &h) { return acc + h.nbr_copies * h.fitness; }) /
                block.population_size;
        }

        ofstream vcf;
        vcf.open(output_filename + "_" + node_name + ".vcf");
        vcf << out << endl;
        vcf.close();

        // pN/pS computed on the sample

        stats.pn_sample = non_syn_nbr;
        stats.ps_sample = syn_nbr;
        stats.pnps_sample = stats.pn_sample / stats.ps_sample;

        stats.pin_sample *= 2.0 / (2 * sample_size * (2 * sample_size - 1));
        stats.pis_sample *= 2.0 / (2 * sample_size * (2 * sample_size - 1));
        stats.pinpis_sample = stats.pin_sample / stats.pis_sample;

        stats.dn_sample = 2 * sum_population_size * dn_from_table(table) / blocks.size();
        stats.ds_sample = 2 * sum_population_size * ds_from_table(table) / blocks.size();
        stats.omega_sample = stats.dn_sample / stats.ds_sample;


        stats.mean_fitness /= blocks.size();
        stats.at_pct = at_sites_obs / (nbr_nucleotides * 2 * sample_size);
        assert(sfs_non_syn.size() == 2 * sample_size + 1);
        assert(sfs_syn.size() == 2 * sample_size + 1);
        assert(sum(sfs_non_syn) == non_syn_nbr);
        assert(sum(sfs_syn) == syn_nbr);

        stats_vector.push_back(stats);

        unsigned nbr_mutations = table.syn_mut + table.non_syn_mut;
        unsigned nbr_fixations = table.syn_fix + table.non_syn_fix;

        trace.add("taxon_name", node_name);
        trace.add("#generations_from_root", time_from_root);
        trace.add("population_size", population_size);
        trace.add("generation_time_in_year", generation_time);
        trace.add("mutation_rate_per_generation", nuc_matrix.mutation_rate);
        log_multivariate.add_to_trace(trace);
        trace.add("dN", stats.dn);
        trace.add("dS", stats.ds);
        trace.add("dN_sample", stats.dn_sample);
        trace.add("dS_sample", stats.ds_sample);
        trace.add("dN/dS", stats.omega);
        trace.add("dN/dS_sample", stats.omega_sample);
        trace.add("dN/dS_predicted", stats.omega_predicted);
        trace.add("pN_sample", stats.pn_sample);
        trace.add("pS_sample", stats.ps_sample);
        trace.add("pN/pS_sample", stats.pnps_sample);
        trace.add("piN_sample", stats.pin_sample);
        trace.add("piS_sample", stats.pis_sample);
        trace.add("piN/piS_sample", stats.pinpis_sample);
        trace.add("#mutations", nbr_mutations);
        trace.add("#fixations", nbr_fixations);
        trace.add("#non_synonymous_mutations", table.non_syn_mut);
        trace.add("#synonymous_mutations", table.syn_mut);
        trace.add("#non_synonymous_fixations", table.non_syn_fix);
        trace.add("#synonymous_fixations", table.syn_fix);
        trace.add("MeanFitness", stats.mean_fitness);
        trace.add("#complex_sites", complex_sites);
        trace.add("%AT", stats.at_pct);
        trace.add("SFSn", join(sfs_non_syn, ' '));
        trace.add("E[SFSn]", mean(sfs_non_syn));
        trace.add("SFSs", join(sfs_syn, ' '));
        trace.add("E[SFSs]", mean(sfs_syn));

        time_elapsed.exportation += duration(timeNow() - t_start);
    }

    vector<double> theoretial_sfs(
        NucleotideRateMatrix const &p, unsigned nbr_sample, bool synonymous) {
        unsigned precision = 8;
        unsigned nbr_points = pow_int(2, precision) + 1;

        vector<unsigned> sample_range(nbr_sample - 1, 0);
        iota(sample_range.begin(), sample_range.end(), 1);

        vector<unsigned long> binom_coeff(sample_range.size(), 0);
        vector<double> sample_sfs(sample_range.size(), 0);

        for (size_t index{0}; index < sample_range.size(); index++) {
            binom_coeff[index] = binomial_coefficient(nbr_sample, sample_range[index]);
        }

        for (auto const &block : blocks) {
            double theta = 4 * block.population_size * p.mutation_rate;

            for (auto const &aa_fitness_profil : block.aa_fitness_profiles) {
                array<double, 64> codon_freqs =
                    codon_frequencies(aa_fitness_profil, p, 4 * population_size);

                double x = 0.0;
                double x_max = 1.0;
                double h = (x_max - x) / (nbr_points - 1);

                vector<double> x_array(nbr_points, 0);
                vector<double> y_array(nbr_points, 0);

                for (unsigned point{0}; point < nbr_points; point++) {
                    x += h;
                    x_array[point] = x;

                    double res = 0;
                    for (char codon_from{0}; codon_from < 64; codon_from++) {
                        if (Codon::codon_to_aa_array[codon_from] != 20) {
                            double tmp_res = 0;

                            for (auto &neighbor : Codon::codon_to_neighbors_array[codon_from]) {
                                char codon_to{0}, n_from{0}, n_to{0};
                                tie(codon_to, n_from, n_to) = neighbor;
                                assert(n_from != n_to);

                                char aa_from = Codon::codon_to_aa_array[codon_from];
                                char aa_to = Codon::codon_to_aa_array[codon_to];

                                if (((synonymous) and (aa_from == aa_to)) or
                                    ((!synonymous) and (aa_from != aa_to))) {
                                    double pij = theta;

                                    if (synonymous) {
                                        pij *= 1 - x;
                                    } else {
                                        double s = aa_fitness_profil[aa_to];
                                        s -= aa_fitness_profil[aa_from];
                                        s *= 4 * population_size;
                                        if (fabs(s) <= Codon::epsilon) {
                                            pij *= 1 - x;
                                        } else {
                                            pij *= (1 - exp(-s * (1 - x))) / (1 - exp(-s));
                                        }
                                    }
                                    pij *= p(n_from, n_to);
                                    tmp_res += pij;
                                }
                            }
                            res += codon_freqs[codon_from] * tmp_res;
                        }
                    }

                    y_array[point] = res;

                    for (size_t index{0}; index < sample_range.size(); index++) {
                        unsigned a = sample_range[index];
                        double y = y_array[point] * pow(x, a - 1) * pow(1 - x, nbr_sample - a - 1);
                        sample_sfs[index] += binom_coeff[index] * h * y;
                    }
                }
            }
        }
        return sample_sfs;
    }

    // Method returning the DNA string corresponding to the codon sequence.
    string get_dna_str() const {
        unsigned nbr_sites = 0;
        for (auto &block : blocks) { nbr_sites += block.nbr_sites; }

        // The DNA string is 3 times larger than the codon sequence.
        string dna_str{};
        dna_str.reserve(nbr_sites * 3);
        // For each site of the sequence.
        for (auto const &block : blocks) {
            for (unsigned site{0}; site < block.nbr_sites; site++) {
                // Assert there is no stop in the sequence.
                assert(Codon::codon_to_aa_array[block.codon_seq[site]] != 20);

                // Translate the site to a triplet of DNA nucleotides
                array<char, 3> triplet = Codon::codon_to_triplet_array[block.codon_seq[site]];
                for (char position{0}; position < 3; position++) {
                    dna_str += Codon::nucleotides[triplet[position]];
                }
            }
        }
        return dna_str;  // return the DNA sequence as a string.
    }

    static void avg_summary_statistics() {
        cout << Block::nbr_fixations << " fixations" << endl;
        cout << Block::nbr_mutations << " mutations" << endl;
        average(stats_vector);
        double total_time = time_elapsed.mutation + time_elapsed.selection +
                            time_elapsed.extinction + time_elapsed.fixation +
                            time_elapsed.exportation + time_elapsed.correlation +
                            time_elapsed.matrix + time_elapsed.exp + time_elapsed.normal;
        cout << setprecision(3) << total_time / 1e9 << "s total time" << endl;
        cout << 100 * time_elapsed.mutation / total_time
             << "% of time spent in calculating mutation (" << time_elapsed.mutation / 1e9 << "s)"
             << endl;
        cout << 100 * time_elapsed.selection / total_time
             << "% of time spent in calculating selection (" << time_elapsed.selection / 1e9 << "s)"
             << endl;
        cout << 100 * time_elapsed.extinction / total_time
             << "% of time spent in calculating extinction (" << time_elapsed.extinction / 1e9
             << "s)" << endl;
        cout << 100 * time_elapsed.fixation / total_time
             << "% of time spent in calculating fixation (" << time_elapsed.fixation / 1e9 << "s)"
             << endl;
        cout << 100 * time_elapsed.exportation / total_time
             << "% of time spent in exportationing vcf (" << time_elapsed.exportation / 1e9 << "s)"
             << endl;
        cout << 100 * time_elapsed.correlation / total_time
             << "% of time spent in the correlation (" << time_elapsed.correlation / 1e9 << "s)"
             << endl;
        cout << 100 * time_elapsed.normal / total_time << "% of time spent in the normal draw ("
             << time_elapsed.normal / 1e9 << "s)" << endl;
        cout << 100 * time_elapsed.matrix / total_time << "% of time spent in the matrix * ("
             << time_elapsed.matrix / 1e9 << "s)" << endl;
        cout << 100 * time_elapsed.exp / total_time << "% of time spent in the exp ("
             << time_elapsed.exp / 1e9 << "s)" << endl;
    }
};

// Initialize static variables
unsigned Block::nbr_mutations = 0, Block::nbr_fixations = 0;
vector<SummaryStatistics> Population::stats_vector{};

class Process {
  private:
    static double years_computed;
    const Tree &tree;
    vector<Population *> populations;  // Vector of sequence DNA.

  public:
    // Constructor
    Process(const Tree &intree, Population &root_pop) : tree{intree}, populations() {
        populations.resize(tree.nb_nodes());
        populations[tree.root()] = &root_pop;
    }

    void run(string &output_filename) { run_recursive(tree.root(), output_filename); }


    // Recursively iterate through the subtree.
    void run_recursive(Tree::NodeIndex node, string &output_filename) {
        // Substitutions of the DNA sequence is generated.

        populations[node]->run_forward(tree.node_length(node));
        years_computed += tree.node_length(node);
        cout << years_computed << " years computed ("
             << static_cast<int>(100 * years_computed / tree.total_length()) << "%)." << endl;

        string name = tree.node_name(node);
        if (tree.is_leaf(node)) {
            // If the node is a leaf, output the DNA sequence and name.
            write_sequence(output_filename, name, populations[node]->get_dna_str());
            populations[node]->output_vcf(output_filename, name);
        } else {
            // If the node is internal, iterate through the direct children.
            for (auto &child : tree.children(node)) {
                populations[child] = new Population(*populations[node]);
                run_recursive(child, output_filename);
            }
        }
    }
};

// Initialize static variables
double Process::years_computed = 0.0;

class SimuPolyArgParse : public SimuArgParse {
  public:
    explicit SimuPolyArgParse(CmdLine &cmd) : SimuArgParse(cmd) {}

    TCLAP::ValueArg<std::string> correlation_path{
        "c", "correlation_matrix", "input correlation matrix path", false, "", "string", cmd};
    TCLAP::ValueArg<double> mu{"m", "mu", "Mutation rate", false, 1e-8, "double", cmd};
    TCLAP::ValueArg<double> root_age{
        "a", "root_age", "Age of the root", false, 50e6, "double", cmd};
    TCLAP::ValueArg<double> generation_time{
        "g", "generation_time", "The number of year between generations", false, 40, "double", cmd};
    TCLAP::ValueArg<unsigned> pop_size{
        "n", "pop_size", "Population size", false, 500, "unsigned", cmd};
    TCLAP::ValueArg<unsigned> sample_size{
        "p", "sample_size", "Sample size", false, 20, "unsigned", cmd};
    TCLAP::ValueArg<double> beta{
        "b", "beta", "Effective population size (relative)", false, 1.0, "double", cmd};
    SwitchArg linked{"l", "linked", "Sites are genetically linked", cmd, false};
};

int main(int argc, char *argv[]) {
    CmdLine cmd{"SimuPoly", ' ', "0.1"};
    SimuPolyArgParse args(cmd);
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
    unsigned pop_size{args.pop_size.getValue()};
    unsigned sample_size{args.sample_size.getValue()};
    assert(sample_size <= pop_size);
    assert(sample_size > 0);
    double beta{args.beta.getValue()};
    assert(beta > 0.0);
    bool linked_sites{args.linked.getValue()};


    vector<array<double, 20>> fitness_profiles =
        open_preferences(preferences_path, beta / (4 * pop_size));

    Tree tree(newick_path);
    tree.set_root_age(root_age);

    unsigned burn_in = 100 * pop_size;
    NucleotideRateMatrix nuc_matrix(nuc_matrix_path, mu, true);

    LogMultivariate log_multivariate(pop_size, generation_time, mu);
    CorrelationMatrix correlation_matrix(correlation_path);
    correlation_matrix /= root_age;

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
    parameters.add("#generations_burn_in", burn_in);
    parameters.add("population_size", pop_size);
    parameters.add("sample_size", sample_size);
    parameters.add("linked_sites", linked_sites);
    correlation_matrix.add_to_trace(parameters);
    parameters.write_tsv(output_path + ".parameters");

    init_alignments(output_path, tree.nb_leaves(), fitness_profiles.size() * 3);
    Population root_population(fitness_profiles, sample_size, log_multivariate, linked_sites,
        nuc_matrix, correlation_matrix);
    root_population.burn_in(burn_in);

    Process simu_process(tree, root_population);
    simu_process.run(output_path);

    trace.write_tsv(output_path);

    Population::avg_summary_statistics();

    cout << "Simulation computed." << endl;
    cout << "Statistics summarized in: " << output_path + ".tsv" << endl;
    cout << "Fasta file in: " << output_path + ".fasta" << endl;
    cout << "Alignment (.ali) file in: " << output_path + ".ali" << endl;
    return 0;
}