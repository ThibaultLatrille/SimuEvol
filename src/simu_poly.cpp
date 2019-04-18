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

static string info =
    "##FILTER=<ID=s50,Description=\"The alternative is the major allele\">\n"
    "##FILTER=<ID=s100,Description=\"The reference has not been sampled\">\n"
    "##INFO=<ID=REFCODON,Number=1,Type=String,Description=\"Codon of the reference\">\n"
    "##INFO=<ID=ALTCODON,Number=1,Type=String,Description=\"Codon of the alternative\">\n"
    "##INFO=<ID=REFAA,Number=1,Type=Character,Description=\"Amino-acid of the reference\">\n"
    "##INFO=<ID=ALTAA,Number=1,Type=Character,Description=\"Amino-acid of the alternative\">\n"
    "##INFO=<ID=POSITION,Number=1,Type=Integer,Description=\"Mutated codon position\">\n"
    "##INFO=<ID=ALTCOUNT,Number=1,Type=Integer,Description=\"Number of alternative copy\">\n"
    "##INFO=<ID=SYN,Number=1,Type=boolean,Description=\"Is a synonymous mutation\">\n"
    "##FORMAT=<ID=GT,Number=1,Type=String,Description=\"Genotype\">\n"
    "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT";

bool is_synonymous(char codon_from, char codon_to) {
    return Codon::codon_to_aa_array.at(codon_from) == Codon::codon_to_aa_array.at(codon_to);
}

typedef binomial_distribution<u_long> BinomialDistr;
typedef unordered_map<u_long, BinomialDistr> MapBinomialDistr;

class Events {
  public:
    u_long non_syn_mut{0};
    u_long syn_mut{0};
    u_long non_syn_fix{0};
    u_long syn_fix{0};

    explicit Events() = default;

    void add(Events const &e) {
        non_syn_mut += e.non_syn_mut;
        syn_mut += e.syn_mut;
        non_syn_fix += e.non_syn_fix;
        syn_fix += e.syn_fix;
    }

    void clear() {
        non_syn_mut = 0;
        syn_mut = 0;
        non_syn_fix = 0;
        syn_fix = 0;
    }

    void add_to_trace(Trace &trace) const {
        trace.add("Branch_dS", ds());
        trace.add("Branch_dN", dn());
        trace.add("Branch_dNdS", dnds());
        trace.add("Branch_#mutations", mutations());
        trace.add("Branch_#fixations", fixations());
        trace.add("Branch_#non_synonymous_mutations", non_syn_mut);
        trace.add("Branch_#synonymous_mutations", syn_mut);
        trace.add("Branch_#non_synonymous_fixations", non_syn_fix);
        trace.add("Branch_#synonymous_fixations", syn_fix);
    }

    void add_to_tree(Tree &tree, Tree::NodeIndex node) const {
        tree.set_tag(node, "Branch_dS", d_to_string(ds()));
        tree.set_tag(node, "Branch_dN", d_to_string(dn()));
        tree.set_tag(node, "Branch_dNdS", d_to_string(dnds()));
    }

  private:
    double dn() const { return static_cast<double>(non_syn_fix) / non_syn_mut; }

    double ds() const { return static_cast<double>(syn_fix) / syn_mut; };

    double dnds() const {
        if (syn_fix == 0 or non_syn_mut == 0) {
            return 0.0;
        } else {
            return dn() / ds();
        }
    }

    u_long fixations() const { return non_syn_fix + syn_fix; }

    u_long mutations() const { return non_syn_mut + syn_mut; }
};

class Substitution {
  public:
    Substitution(char codon_from, char codon_to, double time_event, double time_between,
        u_long non_syn_mut_flow, u_long syn_mut_flow, double delta_logfitness,
        double mean_logfitness)
        : codon_from{codon_from},
          codon_to{codon_to},
          time_event{time_event},
          time_between{time_between},
          non_syn_mut_flow{non_syn_mut_flow},
          syn_mut_flow{syn_mut_flow},
          delta_logfitness{delta_logfitness},
          mean_logfitness{mean_logfitness} {}

    bool is_synonymous() const { return ::is_synonymous(codon_from, codon_to); }

    bool is_non_synonymous() const { return !::is_synonymous(codon_from, codon_to); }

    void add_to_trace(Trace &trace) const {
        trace.add("codon_from", Codon::codon_string(codon_from));
        trace.add("codon_to", Codon::codon_string(codon_to));
        trace.add("aa_from", Codon::codon_aa_string(codon_from));
        trace.add("aa_to", Codon::codon_aa_string(codon_to));
        trace.add("synonymous", is_synonymous());
        trace.add("selection_coefficient", delta_logfitness);
        trace.add("mean_fitness", mean_logfitness);
        trace.add("time_event", time_event);
        trace.add("time_between", time_between);
    }
    char codon_from;
    char codon_to;
    double time_event;
    double time_between;
    u_long non_syn_mut_flow;
    u_long syn_mut_flow;
    double delta_logfitness;
    double mean_logfitness;
};

class Polymorphism {
  private:
    u_long sample_size;
    double inv_harmonic_sum;
    double inv_choice;

  public:
    u_long pairwise_syn{0};
    u_long pairwise_non_syn{0};
    u_long non_syn_nbr{0};
    u_long syn_nbr{0};
    u_long complex_sites{0};
    double syn_flow{0.0};
    double non_syn_flow{0.0};

    explicit Polymorphism(u_long insample_size)
        : sample_size{insample_size} {
        double harmonic_sum = 0;
        for (u_long i = 1; i < 2 * sample_size; i++) { harmonic_sum += 1.0 / i; }
        inv_harmonic_sum = 1.0 / harmonic_sum;
        inv_choice = 2.0 / (2 * sample_size * (2 * sample_size - 1));
    };

    void add(Polymorphism const &poly) {
        pairwise_syn += poly.pairwise_syn;
        pairwise_non_syn += poly.pairwise_non_syn;
        non_syn_nbr += poly.non_syn_nbr;
        syn_nbr += poly.syn_nbr;
        complex_sites += poly.complex_sites;
        syn_flow += poly.syn_flow;
        non_syn_flow += poly.non_syn_flow;
    }

    void add_to_trace(Trace &trace) const {
        trace.add("theta_watterson_non_syn", theta_watterson_non_syn());
        trace.add("theta_watterson_syn", theta_watterson_syn());
        trace.add("theta_watterson", theta_watterson());
        trace.add("theta_pairwise_non_syn", theta_pairwise_non_syn());
        trace.add("theta_pairwise_syn", theta_pairwise_syn());
        trace.add("theta_pairwise", theta_pairwise());
    }

    void add_to_tree(Tree &tree, Tree::NodeIndex node) const {
        tree.set_tag(node, "theta_watterson_non_syn", d_to_string(theta_watterson_non_syn()));
        tree.set_tag(node, "theta_watterson_syn", d_to_string(theta_watterson_syn()));
        tree.set_tag(node, "theta_watterson", d_to_string(theta_watterson()));
        tree.set_tag(node, "theta_pairwise_non_syn", d_to_string(theta_pairwise_non_syn()));
        tree.set_tag(node, "theta_pairwise_syn", d_to_string(theta_pairwise_syn()));
        tree.set_tag(node, "theta_pairwise", d_to_string(theta_pairwise()));
    }

  private:
    double theta_watterson_non_syn() const { return non_syn_nbr * inv_harmonic_sum / non_syn_flow; }
    double theta_watterson_syn() const { return syn_nbr * inv_harmonic_sum / syn_flow; }
    double theta_watterson() const { return theta_watterson_non_syn() + theta_watterson_syn(); }
    double theta_pairwise_non_syn() const { return pairwise_non_syn * inv_choice / non_syn_flow; }
    double theta_pairwise_syn() const { return pairwise_syn * inv_choice / syn_flow; }
    double theta_pairwise() const { return theta_pairwise_non_syn() + theta_pairwise_syn(); }
};

struct TimeElapsed {
    double mutation{0.0};
    double selection{0.0};
    double extinction{0.0};
    double fixation{0.0};
    double exportation{0.0};
    double correlation{0.0};
};

TimeElapsed timer;
Trace tracer_nodes;
Trace tracer_leaves;
Trace tracer_substitutions;
Trace tracer_traits;

class Haplotype {
  public:
    // The nbr_copies of sites in the sequence (each position is a codon, thus the DNA sequence is 3
    // times greater).
    u_long nbr_copies{0};
    double fitness{0.0};
    unordered_map<u_long, char> diff_sites;

    explicit Haplotype() = default;

    // Constructor of Reference_seq.
    // size: the size of the DNA sequence.
    Haplotype(u_long const nbr_copies, double const fitness)
        : nbr_copies{nbr_copies}, fitness{fitness} {};

    bool check_consistency(u_long nbr_sites) const {
        for (auto &diff : diff_sites) {
            if (Codon::codon_to_aa_array.at(diff.second) == 20) {
                cerr << "The haplotype contains a stop codon" << endl;
                return false;
            }
        }
        return true;
    }

    static struct {
        bool operator()(const Haplotype &left, const Haplotype &right) {
            return left.nbr_copies > right.nbr_copies;
        }

        bool operator()(const Haplotype &left, float right) { return left.nbr_copies > right; }

        bool operator()(float left, const Haplotype &right) { return left > right.nbr_copies; }
    } GreaterThan;
};

// Class representing a genetically linked sequences (exons are unlinked between them)
class Exon {
  public:
    // Exon
    u_long position;

    // Selection
    vector<array<double, 20>> aa_fitness_profiles;

    // Reference sequence
    u_long nbr_sites;
    u_long nbr_nucleotides;
    vector<char> codon_seq;

    // Haplotypes
    vector<Haplotype> haplotype_vector;

    // Keep track of mutations and fixations
    Events events{};

    // Statics variables (shared by all instances)
    static u_long nbr_mutations, nbr_fixations;

    // Constructor
    explicit Exon(vector<array<double, 20>> const &fitness_profiles, u_long const &position,
        u_long const &population_size, NucleotideRateMatrix const &nuc_matrix)
        : position{position},
          aa_fitness_profiles{fitness_profiles},
          nbr_sites{u_long(fitness_profiles.size())},
          nbr_nucleotides{u_long(3 * fitness_profiles.size())},
          codon_seq(fitness_profiles.size(), 0),
          haplotype_vector{} {
        // Draw codon from codon frequencies
        for (u_long site{0}; site < nbr_sites; site++) {
            array<double, 64> codon_freqs =
                codon_frequencies(aa_fitness_profiles.at(site), nuc_matrix, 4 * population_size);
            discrete_distribution<char> freq_codon_distr(codon_freqs.begin(), codon_freqs.end());
            codon_seq[site] = freq_codon_distr(generator);
        }

        haplotype_vector.emplace_back(2 * population_size, 0.0);
        assert(check_consistency(population_size));
    }

    bool check_consistency(u_long const &population_size) const {
        if (haplotype_vector.size() > 2 * population_size) {
            cerr << "Too many haplotypes" << endl;
            return false;
        }
        if (haplotype_vector.empty()) {
            cerr << "No haplotype" << endl;
            return false;
        }

        u_long nbr_copies{0};
        for (auto &haplotype : haplotype_vector) {
            if (!haplotype.check_consistency(nbr_sites)) {
                cerr << "The haplotype is not consistent" << endl;
                return false;
            }
            nbr_copies += haplotype.nbr_copies;
        }
        if (nbr_copies != 2 * population_size) {
            cerr << "The number of copies is not equal to the population size." << endl;
            return false;
        }
        return true;
    }

    void forward(NucleotideRateMatrix const &nuc_matrix, u_long const &population_size,
        BinomialDistr &binomial_distrib, double time_current, vector<Substitution> &substitutions,
        u_long &non_syn_mutations, u_long &syn_mutations, bool trace) {
        mutation(nuc_matrix, binomial_distrib, non_syn_mutations, syn_mutations, trace);
        selection_and_drift(population_size);
        extinction();
        fixation(time_current, substitutions, non_syn_mutations, syn_mutations, trace);
    }

    u_long random_nuc_site() const {
        return uniform_int_distribution<u_long>(0, nbr_nucleotides - 1)(generator);
    }

    void mutation(NucleotideRateMatrix const &p, BinomialDistr &binomial_distrib,
        u_long &non_syn_mutations, u_long &syn_mutations, bool trace) {
        TimeVar t_start = timeNow();

        // Randomly choose the number of mutations at this generation (for this exon)
        u_long binomial_draw = binomial_distrib(generator);

        // Early break if 0 mutation are drawn
        if (binomial_draw == 0) { return; }

        // Compute the distribution of haplotype frequency in the population
        discrete_distribution<u_long> haplotype_distri;
        if (haplotype_vector.size() != 1) {
            vector<u_long> nbr_copies(haplotype_vector.size(), 0);
            std::transform(haplotype_vector.begin(), haplotype_vector.end(), nbr_copies.begin(),
                [](const Haplotype &h) { return h.nbr_copies; });
            haplotype_distri = discrete_distribution<u_long>(nbr_copies.begin(), nbr_copies.end());
        }

        for (u_long nbr_draws{0}; nbr_draws < binomial_draw; nbr_draws++) {
            // Randomly choose the haplotype
            u_long hap_id{0};
            if (haplotype_vector.size() != 1) { hap_id = haplotype_distri(generator); }

            // Early continue (next iteration) if the selected haplotype is not held by any
            // individual
            if (haplotype_vector.at(hap_id).nbr_copies == 0) { continue; }

            // Randomly choose the site
            u_long nuc_site = random_nuc_site();
            u_long codon_site = nuc_site / 3;
            auto nuc_pos = static_cast<char>(nuc_site % 3);

            // The corresponding codon given the selected haplotype and site
            char codon_from{0};
            if (haplotype_vector.at(hap_id).diff_sites.count(codon_site) > 0) {
                codon_from = haplotype_vector.at(hap_id).diff_sites.at(codon_site);
            } else {
                codon_from = codon_seq.at(codon_site);
            }

            // The selected nucleotide
            array<char, 3> triplet_nuc = Codon::codon_to_triplet_array.at(codon_from);
            char nuc_from = triplet_nuc.at(nuc_pos);

            // Early continue (next iteration) if the rate away from this nucleotide is too low
            // (random uniform) compared to the maximum rate away
            bool no_mutation = false;
            if (p.max_sum_mutation_rates != p.sum_mutation_rates(nuc_from)) {
                double sum_unif = p.max_real_distr(generator);
                if (sum_unif > p.sum_mutation_rates(nuc_from)) { no_mutation = true; }
            }
            if (no_mutation) { continue; }

            // Randomly choose the target nucleotide
            triplet_nuc[nuc_pos] = p.mutation_distr.at(nuc_from)(generator);
            char codon_to =
                Codon::triplet_to_codon(triplet_nuc.at(0), triplet_nuc.at(1), triplet_nuc.at(2));
            if (Codon::codon_to_aa_array.at(codon_to) == 20) { continue; }

            // Depending on whether the mutation is back to the reference sequence
            Haplotype haplotype = haplotype_vector.at(hap_id);
            if (codon_to == codon_seq.at(codon_site)) {
                assert(haplotype.diff_sites.count(codon_site) > 0);
                haplotype.diff_sites.erase(codon_site);
            } else {
                haplotype.diff_sites[codon_site] = codon_to;
            }

            // Update the fitness of the new haplotype
            double df =
                aa_fitness_profiles.at(codon_site).at(Codon::codon_to_aa_array.at(codon_to)) -
                aa_fitness_profiles.at(codon_site).at(Codon::codon_to_aa_array.at(codon_from));
            haplotype.fitness += df;

            // Update the vector of haplotypes
            haplotype_vector.at(hap_id).nbr_copies--;
            haplotype.nbr_copies = 1;
            haplotype_vector.push_back(haplotype);

            // Track the created mutation
            if (!trace) { continue; }
            if (is_synonymous(codon_from, codon_to)) {
                syn_mutations++;
                events.syn_mut++;
            } else {
                non_syn_mutations++;
                events.non_syn_mut++;
            }
            nbr_mutations++;
        }
        timer.mutation += duration(timeNow() - t_start);
    }

    void selection_and_drift(u_long const &population_size) {
        TimeVar t_start = timeNow();
        u_long children_tot = 2 * population_size;

        // Early break if only 1 haplotype
        if (haplotype_vector.size() == 1) {
            haplotype_vector.front().nbr_copies = children_tot;
            return;
        }

        // The fitness associated to each haplotype (weigthed by the number of copies)
        vector<double> fitnesses(haplotype_vector.size(), 0);
        std::transform(haplotype_vector.begin(), haplotype_vector.end(), fitnesses.begin(),
            [](const Haplotype &h) { return (1.0 + h.fitness) * h.nbr_copies; });
        double fit_tot = accumulate(fitnesses.begin(), fitnesses.end(), 0.0);

        // Random draws from the multinomial distribution
        for (size_t i_hap{0}; i_hap < haplotype_vector.size() - 1; i_hap++) {
            haplotype_vector.at(i_hap).nbr_copies = binomial_distribution<u_long>(
                children_tot, fitnesses.at(i_hap) / fit_tot)(generator);
            children_tot -= haplotype_vector.at(i_hap).nbr_copies;
            fit_tot -= fitnesses.at(i_hap);
        }
        assert(children_tot <= 2 * population_size);
        assert(children_tot >= 0);
        haplotype_vector.back().nbr_copies = children_tot;
        timer.selection += duration(timeNow() - t_start);
    }

    void extinction() {
        TimeVar t_start = timeNow();
        // Early break if only 1 haplotype
        if (haplotype_vector.size() == 1) { return; }

        // Sort the vector of haplotypes by number of copies
        if (!is_sorted(haplotype_vector.begin(), haplotype_vector.end(), Haplotype::GreaterThan)) {
            sort(haplotype_vector.begin(), haplotype_vector.end(), Haplotype::GreaterThan);
        }

        // Remove haplotypes with 0 copies
        if (haplotype_vector.back().nbr_copies == 0) {
            auto low_bound = lower_bound(
                haplotype_vector.begin(), haplotype_vector.end(), 0, Haplotype::GreaterThan);
            if (low_bound != haplotype_vector.end()) {
                haplotype_vector.erase(low_bound, haplotype_vector.end());
            }
        }
        timer.extinction += duration(timeNow() - t_start);
    }

    void fixation(double time_current, vector<Substitution> &substitutions,
        u_long &non_syn_mutations, u_long &syn_mutations, bool trace) {
        TimeVar t_start = timeNow();

        // ! Copy before the loop is necessary, since it can be updated during the loop
        auto diffs = haplotype_vector.begin()->diff_sites;

        // For each site (different to the reference) of the most common haplotype
        for (auto const &diff : diffs) {
            u_long site = diff.first;
            char codon_to = diff.second;
            bool polymorphic_site = false;
            set<char> poly_codons{codon_to};
            for (size_t hap_id{1}; hap_id < haplotype_vector.size(); hap_id++) {
                if (haplotype_vector.at(hap_id).diff_sites.count(site) == 0) {
                    // In the case the site of this haplotype is the same than the reference
                    polymorphic_site = true;
                    break;
                } else {
                    // In the case the site of this haplotype is different from the reference
                    poly_codons.insert(haplotype_vector.at(hap_id).diff_sites.at(site));
                }
            }

            // Early continue (next iteration) if the site is polymorphic
            if (polymorphic_site) { continue; }

            char codon_from = codon_seq.at(site);
            if (poly_codons.size() > 1) {
                // In the case the reference codon is not present in the alternative haplotypes,
                // but they are several alternative codons, we have to find the most common.
                unordered_map<char, u_long> codon_copies{};
                for (char codon : poly_codons) { codon_copies[codon] = 0; }
                for (auto const &haplotype : haplotype_vector) {
                    codon_copies.at(haplotype.diff_sites.at(site)) += haplotype.nbr_copies;
                }
                codon_to = std::max_element(codon_copies.begin(), codon_copies.end(),
                    [](auto const &p1, auto const &p2) { return p1.second < p2.second; })
                               ->first;
            }
            codon_seq[site] = codon_to;

            // Update the vector of haplotypes
            double df = aa_fitness_profiles.at(site).at(Codon::codon_to_aa_array.at(codon_to)) -
                        aa_fitness_profiles.at(site).at(Codon::codon_to_aa_array.at(codon_from));

            double mean_fitness = 0;
            for (auto &haplotype : haplotype_vector) {
                mean_fitness += haplotype.fitness;
                haplotype.fitness -= df;
                if (haplotype.diff_sites.at(site) == codon_to) {
                    haplotype.diff_sites.erase(site);
                } else {
                    assert(poly_codons.size() > 1);
                    assert(poly_codons.find(haplotype.diff_sites.at(site)) != poly_codons.end());
                }
            }
            mean_fitness /= haplotype_vector.size();

            // Track the created substitution
            if (!trace) { continue; }
            if (is_synonymous(codon_from, codon_to)) {
                events.syn_fix++;
            } else {
                events.non_syn_fix++;
            }
            nbr_fixations++;
            double t_between = time_current;
            if (!substitutions.empty()) {
                t_between -= substitutions.back().time_event;
                if (substitutions.back().time_event == time_current) {
                    non_syn_mutations = substitutions.back().non_syn_mut_flow;
                    syn_mutations = substitutions.back().syn_mut_flow;
                }
            }
            auto sub = Substitution(codon_from, codon_to, time_current, t_between,
                non_syn_mutations, syn_mutations, df, mean_fitness);
            non_syn_mutations = 0;
            syn_mutations = 0;
            sub.add_to_trace(tracer_substitutions);
            substitutions.push_back(sub);
        }
        timer.fixation += duration(timeNow() - t_start);
    }

    void sample_one_individual() {
        // Compute the distribution of haplotype frequency in the population
        vector<u_long> nbr_copies(haplotype_vector.size(), 0);
        std::transform(haplotype_vector.begin(), haplotype_vector.end(), nbr_copies.begin(),
            [](const Haplotype &h) { return h.nbr_copies; });
        u_long rand_hap =
            discrete_distribution<u_long>(nbr_copies.begin(), nbr_copies.end())(generator);

        auto diffs = haplotype_vector[rand_hap].diff_sites;
        // For each site (different to the reference) of the most common haplotype
        for (auto const &diff : diffs) {
            // Update the vector of haplotypes
            u_long site = diff.first;
            char codon_from = codon_seq[diff.first];
            char codon_to = diff.second;

            double df = aa_fitness_profiles.at(site).at(Codon::codon_to_aa_array.at(codon_to)) -
                        aa_fitness_profiles.at(site).at(Codon::codon_to_aa_array.at(codon_from));

            codon_seq[site] = codon_to;
            for (auto &haplotype : haplotype_vector) {
                haplotype.fitness -= df;
                if (haplotype.diff_sites.count(site) > 0 and
                    haplotype.diff_sites.at(site) == codon_to) {
                    haplotype.diff_sites.erase(site);
                } else if (haplotype.diff_sites.count(site) == 0) {
                    haplotype.diff_sites[site] = codon_from;
                }
            }
        }
        assert(haplotype_vector[rand_hap].diff_sites.empty());
    }

    tuple<double, double> flow(NucleotideRateMatrix const &nuc_matrix) const {
        double non_syn_mut_flow{0.0}, syn_mut_flow{0.0};
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

                // If the mutated amino-acid is a stop codon, the rate of fixation is 0.
                // Else, if the mutated and original amino-acids are non-synonymous, we compute the
                // rate of fixation. Note that, if the mutated and original amino-acids are
                // synonymous, the rate of fixation is 1.
                if (Codon::codon_to_aa_array[codon_to] != 20) {
                    if (Codon::codon_to_aa_array[codon_from] !=
                        Codon::codon_to_aa_array[codon_to]) {
                        non_syn_mut_flow += nuc_matrix.normalized_rate(n_from, n_to);
                    } else {
                        syn_mut_flow += nuc_matrix.normalized_rate(n_from, n_to);
                    }
                }
            }
        }
        return make_tuple(non_syn_mut_flow, syn_mut_flow);
    }
};

class Population {
  public:
    // TimeElapsed
    double time_from_root{0};
    // Exons
    vector<Exon> exons;

    u_long sample_size{0};
    LogMultivariate log_multivariate;

    u_long population_size{0};
    double generation_time{0};
    NucleotideRateMatrix nuc_matrix;
    EMatrix const &transform_matrix;
    bool branch_wise{};

    vector<Substitution> substitutions{};
    u_long non_syn_mutations{0}, syn_mutations{0};
    PieceWiseMultivariate piecewise_multivariate{};

    mutable MapBinomialDistr binomial_distribs;

    Population(vector<array<double, 20>> const &fitness_profiles, u_long sample_size,
        LogMultivariate &log_multi, u_long exon_size, NucleotideRateMatrix const &nucleotide_matrix,
        EMatrix const &transform_matrix, bool branch_wise)
        : exons{},
          sample_size{sample_size},
          log_multivariate{log_multi},
          population_size{log_multivariate.population_size()},
          generation_time{log_multivariate.generation_time()},
          nuc_matrix{nucleotide_matrix},
          transform_matrix{transform_matrix},
          branch_wise{branch_wise} {
        auto dv = std::div(static_cast<int>(fitness_profiles.size()), exon_size);
        if (dv.rem != 0) { dv.quot++; }
        exons.reserve(dv.quot);
        for (int exon{0}; exon < dv.quot; exon++) {
            size_t begin_exon = exon * exon_size;
            size_t end_exon = min(begin_exon + exon_size, fitness_profiles.size());

            std::vector<array<double, 20>> exon_profiles(
                fitness_profiles.begin() + begin_exon, fitness_profiles.begin() + end_exon);

            exons.emplace_back(exon_profiles, begin_exon, population_size, nuc_matrix);
        }
        assert(nbr_sites() == fitness_profiles.size());
        cout << exons.size() << " exons created." << endl;
        if (dv.rem != 0) { cout << "Last exon is " << dv.rem << " sites long." << endl; }
        update_binomial_distribs();
    }

    void update_binomial_distribs() const {
        binomial_distribs.clear();
        for (auto const &exon : exons) {
            if (binomial_distribs.count(exon.nbr_sites) == 0) {
                binomial_distribs[exon.nbr_sites] =
                    binomial_distribution<u_long>(2 * population_size * exon.nbr_nucleotides,
                        nuc_matrix.max_sum_mutation_rates * nuc_matrix.mutation_rate);
            }
        }
    }

    u_long nbr_sites() const {
        u_long sites = 0;
        for (auto const &exon : exons) { sites += exon.nbr_sites; }
        return sites;
    }

    u_long nbr_nucleotides() const { return 3 * nbr_sites(); }

    bool check_consistency() const {
        for (auto const &exon : exons) {
            if (!exon.check_consistency(population_size)) {
                cerr << "The exon is not consistent" << endl;
                return false;
            }
        }
        for (size_t i = 1; i < substitutions.size(); ++i) {
            if (abs(substitutions.at(i).time_between -
                    (substitutions.at(i).time_event - substitutions.at(i - 1).time_event)) > 1e-6) {
                cerr << "The time between substitutions is not consistent" << endl;
                return false;
            }
            if (substitutions.at(i).time_between < 0.0) {
                cerr << "The substitutions are not ordered (in increasing time)" << endl;
                return false;
            }
        }
        return true;
    };

    EVector delta_log_multivariate(double distance) {
        TimeVar t_start = timeNow();
        EVector normal_vector = EVector::Zero(dimensions);
        for (int dim = 0; dim < dimensions; dim++) {
            normal_vector(dim) = normal_distrib(generator);
        }
        timer.correlation += duration(timeNow() - t_start);
        return sqrt(distance) * (transform_matrix * normal_vector);
    }

    void update_brownian(EVector const &delta) {
        TimeVar t_start = timeNow();
        log_multivariate += delta;

        population_size = log_multivariate.population_size();
        generation_time = log_multivariate.generation_time();
        nuc_matrix.set_mutation_rate(log_multivariate.mutation_rate_per_generation());

        if (population_size < sample_size) {
            cerr << "The population size (" << population_size
                 << ") became lower than sample size (" << sample_size << ")" << endl;
            population_size = 2 * sample_size - population_size;
            log_multivariate.set_population_size(population_size);
        }
        update_binomial_distribs();
        timer.correlation += duration(timeNow() - t_start);
    }

    void run_and_trace(
        string const &output_filename, Tree::NodeIndex node, Tree &tree, Population const *parent) {
        double t_max = tree.node_length(node);
        double time_current = 0.0;
        EVector delta;

        if (branch_wise) {
            delta = delta_log_multivariate(t_max / tree.max_distance_to_root());
            update_brownian(delta / 2);
            piecewise_multivariate.AddMultivariate(t_max, log_multivariate);
        }
        while (time_current < t_max) {
            if (!branch_wise) {
                delta = delta_log_multivariate(generation_time / tree.max_distance_to_root());
                update_brownian(delta);
                piecewise_multivariate.AddMultivariate(generation_time, log_multivariate);
            }
            for (auto &exon : exons) {
                exon.forward(nuc_matrix, population_size, binomial_distribs.at(exon.nbr_sites),
                    time_current, substitutions, non_syn_mutations, syn_mutations, true);
            }
            assert(check_consistency());
            time_current += generation_time;
            time_from_root += generation_time;
        }
        node_trace(output_filename, node, tree, parent);
        if (branch_wise) { update_brownian(delta / 2); }
    }

    void burn_in(u_long burn_in_length) {
        cout << "Burn-in for " << burn_in_length << " generations." << endl;
        for (u_long gen{1}; gen <= burn_in_length; gen++) {
            for (auto &exon : exons) {
                assert(!exon.haplotype_vector.empty());
                exon.forward(nuc_matrix, population_size, binomial_distribs.at(exon.nbr_sites), 0,
                    substitutions, non_syn_mutations, syn_mutations, false);
            }
            assert(check_consistency());
        }
        clear_events();
        cout << "Burn-in completed." << endl;
    }

    double theoretical_theta() const { return 4 * population_size * nuc_matrix.mutation_rate; };

    double predicted_dn_dn0(NucleotideRateMatrix const &rates, u_long const &pop_size) const {
        double dn{0.}, dn0{0.};
        for (auto const &exon : exons) {
            double exon_dn{0}, exon_d0{0};
            tie(exon_dn, exon_d0) =
                ::predicted_dn_dn0(exon.aa_fitness_profiles, rates, 4 * pop_size);
            dn += exon_dn;
            dn0 += exon_d0;
        }
        return dn / dn0;
    };

    double sequence_wise_predicted_dn_dn0(
        Population const &parent, NucleotideRateMatrix const &rates, u_long const &pop_size) const {
        double dn{0.}, dn0{0.};

        for (size_t i = 0; i < exons.size(); i++) {
            double exon_dn{0}, exon_dn0{0};
            tie(exon_dn, exon_dn0) = ::flow_dn_dn0(
                exons.at(i).aa_fitness_profiles, parent.exons.at(i).codon_seq, rates, 4 * pop_size);
            dn += exon_dn;
            dn0 += exon_dn0;
        }
        return dn / dn0;
    };

    double count_based_dn_dn0() const {
        double dn{0}, dn0{0}, time_total{0};
        for (auto const &substitution : substitutions) {
            if (substitution.is_non_synonymous()) { dn++; }
            time_total += substitution.time_between;
            dn0 += substitution.non_syn_mut_flow * substitution.time_between;
        }
        return time_total * dn / dn0;
    }

    double event_based_dn_dn0() const {
        double dn{0}, dn0{0};
        for (auto const &substitution : substitutions) {
            if (substitution.is_non_synonymous()) { dn++; }
            dn0 += substitution.non_syn_mut_flow;
        }
        return 2 * population_size * dn / dn0;
    }

    // Simulated omega from the substitutions
    double event_based_dn_ds() const {
        double dn{0}, ds{0}, dn0{0}, ds0{0};
        for (auto const &substitution : substitutions) {
            dn0 += substitution.non_syn_mut_flow;
            ds0 += substitution.syn_mut_flow;
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
            return dn / ds;
        }
    }

    // Simulated omega from the substitutions
    double count_based_dn_ds() const {
        double dn{0}, ds{0}, dn0{0}, ds0{0};
        for (auto const &substitution : substitutions) {
            dn0 += substitution.non_syn_mut_flow * substitution.time_between;
            ds0 += substitution.syn_mut_flow * substitution.time_between;
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

    void node_trace(
        string const &output_filename, Tree::NodeIndex node, Tree &tree, Population const *parent) {
        TimeVar t_start = timeNow();

        string node_name = tree.node_name(node);

        tree.set_tag(node, "population_size", to_string(log_multivariate.population_size()));
        tree.set_tag(node, "generation_time", d_to_string(log_multivariate.generation_time()));
        tree.set_tag(node, "mutation_rate_per_generation",
            d_to_string(log_multivariate.mutation_rate_per_generation()));

        if (tree.is_root(node)) { return; }

        auto sum_events = events();
        sum_events.add_to_tree(tree, node);
        auto geom_pop_size = static_cast<u_long>(piecewise_multivariate.GeometricPopSize());
        piecewise_multivariate.add_to_tree(tree, node, geom_pop_size);

        assert(parent != nullptr);
        tree.set_tag(node, "Branch_dNdN0_predicted",
            d_to_string(predicted_dn_dn0(nuc_matrix, geom_pop_size)));
        tree.set_tag(node, "Branch_dNdN0_sequence_wise_predicted",
            d_to_string(sequence_wise_predicted_dn_dn0(*parent, nuc_matrix, geom_pop_size)));
        tree.set_tag(node, "Branch_dNdN0_event_based", d_to_string(event_based_dn_dn0()));
        tree.set_tag(node, "Branch_dNdN0_count_based", d_to_string(count_based_dn_dn0()));
        tree.set_tag(node, "Branch_dNdS_event_based", d_to_string(event_based_dn_ds()));
        tree.set_tag(node, "Branch_dNdS_count_based", d_to_string(count_based_dn_ds()));

        for (auto const &exon : exons) {
            tracer_nodes.add("taxon_name", node_name);
            tracer_nodes.add("exon_id", exon.position);
            exon.events.add_to_trace(tracer_nodes);
        }

        if (!tree.is_leaf(node)) { return; }

        sample_one_individual();

        // If the node is a leaf, output the DNA sequence and name.
        write_sequence(output_filename, node_name, get_dna_str());

        // VCF file on the sample
        Polymorphism exome_poly(sample_size);

        string out;
        out += "##fileformat=VCFv4.0";
        out += "\n##source=SimuPoly";
        out += "\n##nodeName=" + node_name;
        out += "\n##sequenceSize=" + to_string(nbr_nucleotides());
        out += "\n##ploidyLevel=diploid";
        out += "\n##numberIndividuals=" + to_string(sample_size);
        out += "\n##numberGenotypes=" + to_string(2 * sample_size);
        out += "\n##reference=" + get_dna_str();
        out += "\n" + info;

        for (u_long indiv{1}; indiv <= sample_size; indiv++) {
            out += "\tId";
            out += to_string(indiv);
        }

        for (auto const &exon : exons) {
            Polymorphism exon_poly(sample_size);

            // Draw the sample of individuals
            vector<u_long> haplotypes_sample(2 * population_size, 0);
            u_long sum_copies = 0;
            for (u_long i_hap{0}; i_hap < exon.haplotype_vector.size(); i_hap++) {
                for (u_long copy_id{0}; copy_id < exon.haplotype_vector.at(i_hap).nbr_copies;
                     copy_id++) {
                    assert(sum_copies + copy_id < 2 * population_size);
                    haplotypes_sample[sum_copies + copy_id] = i_hap;
                }
                sum_copies += exon.haplotype_vector.at(i_hap).nbr_copies;
            }
            shuffle(haplotypes_sample.begin(), haplotypes_sample.end(), generator);
            haplotypes_sample.resize(2 * sample_size);

            for (u_long site{0}; site < exon.nbr_sites; site++) {
                map<tuple<char, char>, u_long> codon_from_to_copy{};
                for (auto const &i_hap : haplotypes_sample) {
                    if (exon.haplotype_vector.at(i_hap).diff_sites.count(site) > 0) {
                        char codon_from = exon.codon_seq.at(site);
                        char codon_to = exon.haplotype_vector.at(i_hap).diff_sites.at(site);
                        assert(codon_from != codon_to);
                        codon_from_to_copy[make_tuple(codon_from, codon_to)]++;
                    }
                }

                if (codon_from_to_copy.size() == 1) {
                    char codon_from = get<0>(codon_from_to_copy.begin()->first);
                    char codon_to = get<1>(codon_from_to_copy.begin()->first);
                    assert(codon_to != codon_from);
                    u_long alt_freq = codon_from_to_copy.begin()->second;
                    char nuc_pos{0};
                    while (nuc_pos < 3) {
                        if (Codon::codon_to_nuc(codon_from, nuc_pos) !=
                            Codon::codon_to_nuc(codon_to, nuc_pos)) {
                            break;
                        } else {
                            nuc_pos++;
                        }
                    }
                    assert(nuc_pos != 3);
                    string line{"\n"};
                    line += ".\t";
                    line += to_string(3 * (exon.position + site) + nuc_pos);
                    line += "\t.\t";
                    line += Codon::codon_to_nuc(codon_from, nuc_pos);
                    line += "\t";
                    line += Codon::codon_to_nuc(codon_to, nuc_pos);
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
                    line += to_string(nuc_pos);
                    line += ";ALTCOUNT=";
                    line += to_string(alt_freq);
                    line += ";SYN=";

                    if (is_synonymous(codon_from, codon_to)) {
                        exon_poly.syn_nbr++;
                        line += "TRUE";
                    } else {
                        exon_poly.non_syn_nbr++;
                        line += "FALSE";
                    }

                    line += "\tGT";
                    for (u_long indiv{0}; indiv < sample_size; indiv++) {
                        line += "\t";
                        for (u_long ploidy{0}; ploidy < 2; ploidy++) {
                            u_long i_hap = haplotypes_sample[indiv * 2 + ploidy];
                            char nuc{0};
                            if (exon.haplotype_vector.at(i_hap).diff_sites.count(site) > 0) {
                                char codon = exon.haplotype_vector.at(i_hap).diff_sites.at(site);
                                nuc = Codon::codon_to_nuc(codon, nuc_pos);
                            } else {
                                nuc = Codon::codon_to_nuc(exon.codon_seq.at(site), nuc_pos);
                            }
                            line += nuc;
                            if (ploidy == 0) { line += "|"; }
                        }
                    }
                    out += line;
                } else if (codon_from_to_copy.size() > 1) {
                    exon_poly.complex_sites++;
                }
            }

            // Theta pairwise computed on the sample
            for (u_long i{0}; i < 2 * sample_size; i++) {
                for (u_long j{i + 1}; j < 2 * sample_size; j++) {
                    if (i != j) {
                        u_long hap_i = haplotypes_sample.at(i);
                        u_long hap_j = haplotypes_sample.at(j);
                        auto it_first = exon.haplotype_vector.at(hap_i).diff_sites.begin();
                        auto end_first = exon.haplotype_vector.at(hap_i).diff_sites.end();
                        auto it_second = exon.haplotype_vector.at(hap_j).diff_sites.begin();
                        auto end_second = exon.haplotype_vector.at(hap_j).diff_sites.end();
                        while (it_first != end_first and it_second != end_second) {
                            char codon_from{-1};
                            char codon_to{-1};
                            if (it_second == end_second or (*it_first) < (*it_second)) {
                                codon_from = exon.codon_seq.at(it_first->first);
                                codon_to = it_first->second;
                                it_first++;
                            } else if (it_first == end_first or (*it_second) < (*it_first)) {
                                codon_from = exon.codon_seq.at(it_second->first);
                                codon_to = it_second->second;
                                it_second++;
                            } else if ((*it_first) == (*it_second)) {
                                codon_from = it_first->second;
                                codon_to = it_second->second;
                                it_first++;
                                it_second++;
                            }
                            if (is_synonymous(codon_from, codon_to)) {
                                exon_poly.pairwise_syn++;
                            } else {
                                exon_poly.pairwise_non_syn++;
                            }
                        }
                    }
                }
            }

            tie(exon_poly.non_syn_flow, exon_poly.syn_flow) = exon.flow(nuc_matrix);

            tracer_leaves.add("taxon_name", node_name);
            tracer_leaves.add("exon_id", exon.position);
            exon_poly.add_to_trace(tracer_leaves);

            exome_poly.add(exon_poly);
        }
        ofstream vcf;
        vcf.open(output_filename + "_" + node_name + ".vcf");
        vcf << out << endl;
        vcf.close();

        exome_poly.add_to_tree(tree, node);
        tree.set_tag(node, "theta_pred", d_to_string(theoretical_theta()));

        tracer_traits.add("TaxonName", node_name);
        tracer_traits.add("LogGenerationTime", log_multivariate.log_generation_time());
        tracer_traits.add("LogPopulationSize", log_multivariate.log_population_size());
        tracer_traits.add(
            "LogMutationRatePerGeneration", log_multivariate.log_mutation_rate_per_generation());

        timer.exportation += duration(timeNow() - t_start);
    }

    Events events() const {
        Events e;
        for (auto const &exon : exons) { e.add(exon.events); }
        return e;
    }

    void clear_events() {
        for (auto &exon : exons) { exon.events.clear(); }
        substitutions.clear();
        piecewise_multivariate.clear();
        non_syn_mutations = 0;
        syn_mutations = 0;
    }

    void sample_one_individual() {
        for (auto &exon : exons) { exon.sample_one_individual(); }
    }
    // Method returning the DNA string corresponding to the codon sequence.
    string get_dna_str() const {
        string dna_str{};
        dna_str.reserve(nbr_nucleotides());

        // For each site of the sequence.
        for (auto const &exon : exons) {
            for (u_long site{0}; site < exon.nbr_sites; site++) {
                // Assert there is no stop in the sequence.
                assert(Codon::codon_to_aa_array.at(exon.codon_seq.at(site)) != 20);

                // Translate the site to a triplet of DNA nucleotides
                array<char, 3> triplet = Codon::codon_to_triplet_array.at(exon.codon_seq.at(site));
                for (char nuc_pos{0}; nuc_pos < 3; nuc_pos++) {
                    dna_str += Codon::nucleotides.at(triplet.at(nuc_pos));
                }
            }
        }
        return dna_str;  // return the DNA sequence as a string.
    }

    static void timer_cout() {
        cout << Exon::nbr_fixations << " fixations" << endl;
        cout << Exon::nbr_mutations << " mutations" << endl;
        double total_time = timer.mutation + timer.selection + timer.extinction + timer.fixation +
                            timer.exportation + timer.correlation;
        cout << setprecision(3) << total_time / 1e9 << "s total time" << endl;
        cout << 100 * timer.mutation / total_time << "% of time spent in calculating mutation ("
             << timer.mutation / 1e9 << "s)" << endl;
        cout << 100 * timer.selection / total_time << "% of time spent in calculating selection ("
             << timer.selection / 1e9 << "s)" << endl;
        cout << 100 * timer.extinction / total_time << "% of time spent in calculating extinction ("
             << timer.extinction / 1e9 << "s)" << endl;
        cout << 100 * timer.fixation / total_time << "% of time spent in calculating fixation ("
             << timer.fixation / 1e9 << "s)" << endl;
        cout << 100 * timer.exportation / total_time << "% of time spent in exporting vcf ("
             << timer.exportation / 1e9 << "s)" << endl;
        cout << 100 * timer.correlation / total_time << "% of time spent in the correlation ("
             << timer.correlation / 1e9 << "s)" << endl;
    }
};

// Initialize static variables
u_long Exon::nbr_mutations = 0, Exon::nbr_fixations = 0;

class Process {
  private:
    static double years_computed;
    Tree &tree;
    vector<Population *> populations;  // Vector of sequence DNA.

  public:
    // Constructor
    Process(Tree &intree, Population &root_pop) : tree{intree}, populations() {
        populations.resize(tree.nb_nodes());
        populations[tree.root()] = &root_pop;
    }

    void run(string &output_filename) {
        populations.at(tree.root())->node_trace(output_filename, tree.root(), tree, nullptr);
        run_recursive(tree.root(), output_filename);
        ofstream nhx;
        nhx.open(output_filename + ".nhx");
        nhx << tree.as_string() << endl;
        nhx.close();
    }


    // Recursively iterate through the subtree.
    void run_recursive(Tree::NodeIndex node, string const &output_filename) {
        // Substitutions of the DNA sequence is generated.
        populations.at(node)->clear_events();

        if (!tree.is_root(node)) {
            populations.at(node)->run_and_trace(
                output_filename, node, tree, populations.at(tree.parent(node)));

            years_computed += tree.node_length(node);
            cout << years_computed << " years computed in total ("
                 << static_cast<int>(100 * years_computed / tree.total_length()) << "%) at node "
                 << tree.node_name(node) << " ("
                 << static_cast<int>(100 * tree.node_length(node) / tree.total_length()) << "%)."
                 << endl;
        }

        // Iterate through the direct children.
        for (auto &child : tree.children(node)) {
            populations.at(child) = new Population(*populations.at(node));
            run_recursive(child, output_filename);
        }
    }
};

// Initialize static variables
double Process::years_computed = 0.0;

class SimuPolyArgParse : public SimuArgParse {
  public:
    explicit SimuPolyArgParse(CmdLine &cmd) : SimuArgParse(cmd) {}

    TCLAP::ValueArg<u_long> pop_size{
        "n", "population_size", "Population size (at the root)", false, 500, "u_long", cmd};
    TCLAP::ValueArg<u_long> sample_size{
        "p", "sample_size", "Sample size (at the leaves)", false, 20, "u_long", cmd};
};

int main(int argc, char *argv[]) {
    CmdLine cmd{"SimuPoly", ' ', "0.1"};
    SimuPolyArgParse args(cmd);
    cmd.parse(argc, argv);

    u_long arg_seed = args.seed.getValue();
    cout << "Random generator seed: " << arg_seed << endl;
    generator.seed(arg_seed);

    string preferences_path{args.preferences_path.getValue()};
    string newick_path{args.newick_path.getValue()};
    string nuc_matrix_path{args.nuc_matrix_path.getValue()};
    string output_path{args.output_path.getValue()};
    string precision_path{args.precision_path.getValue()};
    bool fix_pop_size{args.fix_pop_size.getValue()};
    bool fix_mut_rate{args.fix_mut_rate.getValue()};
    bool fix_gen_time{args.fix_gen_time.getValue()};
    double mutation_rate_per_generation{args.mutation_rate_per_generation.getValue()};
    assert(mutation_rate_per_generation > 0.0);
    double root_age{args.root_age.getValue()};
    assert(root_age > 0.0);
    double generation_time{args.generation_time.getValue()};
    assert(generation_time > 0.0);
    assert(generation_time < root_age);
    double beta{args.beta.getValue()};
    assert(beta >= 0.0);
    bool branch_wise_correlation{args.branch_wise_correlation.getValue()};
    u_long pop_size{args.pop_size.getValue()};
    u_long sample_size{args.sample_size.getValue()};
    assert(sample_size <= pop_size);
    assert(sample_size > 0);

    vector<array<double, 20>> fitness_profiles =
        open_preferences(preferences_path, beta / (4 * pop_size));
    u_long nbr_sites = fitness_profiles.size();
    u_long exon_size{args.exons.getValue()};
    if (exon_size == 0) { exon_size = nbr_sites; }
    assert(0 <= exon_size and exon_size <= nbr_sites);

    Tree tree(newick_path);
    tree.set_root_age(root_age);

    u_long burn_in = 100 * pop_size;
    NucleotideRateMatrix nuc_matrix(nuc_matrix_path, mutation_rate_per_generation, true);

    LogMultivariate log_multivariate(pop_size, mutation_rate_per_generation, generation_time);
    CorrelationMatrix correlation_matrix(precision_path, fix_pop_size, fix_mut_rate, fix_gen_time);

    Trace parameters;
    parameters.add("seed", arg_seed);
    parameters.add("output_path", output_path);
    parameters.add("tree_path", newick_path);
    parameters.add("#tree_nodes", tree.nb_nodes());
    parameters.add("#tree_branches", tree.nb_branches());
    parameters.add("#tree_leaves", tree.nb_leaves());
    parameters.add("tree_ultrametric", tree.is_ultrametric());
    parameters.add("tree_min_distance_to_root_in_year", tree.min_distance_to_root());
    parameters.add("tree_max_distance_to_root_in_year", tree.max_distance_to_root());
    parameters.add("site_preferences_path", preferences_path);
    parameters.add("#codonsites", nbr_sites);
    parameters.add("#nucleotidesites", nbr_sites * 3);
    parameters.add("preferences_beta", beta);
    parameters.add("nucleotide_matrix_path", output_path);
    parameters.add("mutation_rate_per_generation", mutation_rate_per_generation);
    nuc_matrix.add_to_trace(parameters);
    parameters.add("generation_time_in_year", generation_time);
    parameters.add("#generations_burn_in", burn_in);
    parameters.add("population_size", pop_size);
    parameters.add("sample_size", sample_size);
    parameters.add("exon_size", exon_size);
    parameters.add("fix_pop_size", fix_pop_size);
    parameters.add("fix_mut_rate", fix_mut_rate);
    parameters.add("fix_gen_time", fix_gen_time);
    correlation_matrix.add_to_trace(parameters);
    parameters.write_tsv(output_path + ".parameters");


    init_alignments(output_path, tree.nb_leaves(), nbr_sites * 3);
    Eigen::SelfAdjointEigenSolver<EMatrix> eigen_solver(correlation_matrix);
    EMatrix transform_matrix =
        eigen_solver.eigenvectors() * eigen_solver.eigenvalues().cwiseSqrt().asDiagonal();

    Population root_population(fitness_profiles, sample_size, log_multivariate, exon_size,
        nuc_matrix, transform_matrix, branch_wise_correlation);
    root_population.burn_in(burn_in);

    Process simu_process(tree, root_population);
    simu_process.run(output_path);

    tracer_leaves.write_tsv(output_path + ".leaves");
    tracer_nodes.write_tsv(output_path + ".nodes");
    tracer_substitutions.write_tsv(output_path + ".substitutions");
    tracer_traits.write_tsv(output_path + ".traits");


    Population::timer_cout();

    cout << "Simulation computed." << endl;
    cout << nbr_sites * 3 * (mutation_rate_per_generation / generation_time) * tree.total_length()
         << " expected substitutions." << endl;
    cout << Exon::nbr_fixations << " simulated substitutions." << endl;
    cout << "Statistics summarized in: " << output_path + ".tsv" << endl;
    cout << "Fasta file in: " << output_path + ".fasta" << endl;
    cout << "Alignment (.ali) file in: " << output_path + ".ali" << endl;
    return 0;
}