#include <chrono>

#include "clonable_unique_ptr.hpp"
#include "codon.hpp"
#include "fitness.hpp"
#include "io.hpp"
#include "random.hpp"
#include "statistic.hpp"
#include "tree.hpp"

#define duration(a) std::chrono::duration_cast<std::chrono::nanoseconds>(a).count()
#define timeNow() std::chrono::high_resolution_clock::now()

typedef std::chrono::high_resolution_clock::time_point TimeVar;

static std::string info =
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
    return codonLexico.codon_to_aa.at(codon_from) == codonLexico.codon_to_aa.at(codon_to);
}

typedef std::binomial_distribution<u_long> BinomialDistr;
typedef std::unordered_map<u_long, BinomialDistr> MapBinomialDistr;

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
        u_long non_syn_mut_flow, u_long syn_mut_flow)
        : codon_from{codon_from},
          codon_to{codon_to},
          time_event{time_event},
          time_between{time_between},
          non_syn_mut_flow{non_syn_mut_flow},
          syn_mut_flow{syn_mut_flow} {}

    bool is_synonymous() const { return ::is_synonymous(codon_from, codon_to); }

    bool is_non_synonymous() const { return !::is_synonymous(codon_from, codon_to); }

    void add_to_trace(Trace &trace) const {
        trace.add("codon_from", codonLexico.codon_string(codon_from));
        trace.add("codon_to", codonLexico.codon_string(codon_to));
        trace.add("aa_from", codonLexico.codon_aa_string(codon_from));
        trace.add("aa_to", codonLexico.codon_aa_string(codon_to));
        trace.add("synonymous", is_synonymous());
        trace.add("time_event", time_event);
        trace.add("time_between", time_between);
    }

    char codon_from;
    char codon_to;
    double time_event;
    double time_between;
    u_long non_syn_mut_flow;
    u_long syn_mut_flow;
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

    explicit Polymorphism(u_long insample_size) : sample_size{insample_size} {
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
Trace tracer_fossils;

class Haplotype {
  public:
    // The nbr_copies of sites in the sequence (each position is a codon, thus the DNA sequence is 3
    // times greater).
    u_long nbr_copies{0};
    double sel_coeff{0.0};
    std::unordered_map<u_long, char> diff_sites;
    ClonableUniquePtr<FitnessState> fitness_state;

    explicit Haplotype() = default;

    // Constructor of Reference_seq.
    // size: the size of the DNA sequence.
    Haplotype(
        u_long const nbr_copies, double const sel_coeff, ClonableUniquePtr<FitnessState> f_state)
        : nbr_copies{nbr_copies}, sel_coeff{sel_coeff}, fitness_state{std::move(f_state)} {};

    bool check_consistency(u_long nbr_sites) const {
        for (auto &diff : diff_sites) {
            if (codonLexico.codon_to_aa.at(diff.second) == 20) {
                std::cerr << "The haplotype contains a stop codon" << std::endl;
                return false;
            }
        }
        return true;
    }

    struct GreaterThan {
        bool operator()(Haplotype const &left, Haplotype const &right) {
            return left.nbr_copies > right.nbr_copies;
        }

        bool operator()(Haplotype const &left, u_long right) { return left.nbr_copies > right; }

        bool operator()(u_long left, Haplotype const &right) { return left > right.nbr_copies; }
    };
};

// Class representing a genetically linked sequences (exons are unlinked between them)
class Exon {
  public:
    // Exon
    u_long position;

    // Reference sequence
    u_long nbr_sites;
    u_long nbr_nucleotides;
    std::vector<char> codon_seq;

    // The fitness profiles of amino-acids.
    ClonableUniquePtr<FitnessState> fitness_state;

    // Haplotypes
    std::vector<Haplotype> haplotype_vector;

    // Keep track of mutations and fixations
    Events events{};

    // Statics variables (shared by all instances)
    static u_long nbr_mutations, nbr_fixations;

    // Constructor
    Exon(std::unique_ptr<FitnessState> &f_state, u_long const &position,
        u_long const &population_size, NucleotideRateMatrix const &nuc_matrix)
        : position{position},
          nbr_sites{u_long(f_state->nbr_sites())},
          nbr_nucleotides{u_long(3 * f_state->nbr_sites())},
          codon_seq(f_state->nbr_sites(), 0),
          fitness_state{std::move(f_state)},
          haplotype_vector{} {
        // Draw codon from codon frequencies
        fitness_state->update(codon_seq, population_size);
        for (u_long site{0}; site < nbr_sites; site++) {
            std::array<double, 64> codon_freqs = fitness_state->codon_frequencies(
                fitness_state->aa_selection_coefficients(codon_seq, site, population_size),
                nuc_matrix, population_size, 0.0);
            std::discrete_distribution<short> freq_codon_distr(
                codon_freqs.begin(), codon_freqs.end());
            codon_seq[site] = freq_codon_distr(generator);
        }
        fitness_state->update(codon_seq, population_size);
        haplotype_vector.emplace_back(2 * population_size, 0.0, fitness_state);
        assert(check_consistency(population_size));
    }

    bool check_consistency(u_long const &population_size) const {
        if (haplotype_vector.size() > 2 * population_size) {
            std::cerr << "Too many haplotypes" << std::endl;
            return false;
        }
        if (haplotype_vector.empty()) {
            std::cerr << "No haplotype" << std::endl;
            return false;
        }

        u_long nbr_copies{0};
        for (auto &haplotype : haplotype_vector) {
            if (!haplotype.check_consistency(nbr_sites)) {
                std::cerr << "The haplotype is not consistent" << std::endl;
                return false;
            }
            nbr_copies += haplotype.nbr_copies;
        }
        if (nbr_copies != 2 * population_size) {
            std::cerr << "The number of copies is not equal to the population size." << std::endl;
            return false;
        }
        return true;
    }

    void forward(NucleotideRateMatrix const &nuc_matrix, u_long const &population_size,
        BinomialDistr &binomial_distrib, double time_current,
        std::vector<Substitution> &substitutions, u_long &non_syn_mutations, u_long &syn_mutations,
        bool burn_in) {
        mutation(nuc_matrix, binomial_distrib, non_syn_mutations, syn_mutations, burn_in,
            population_size);
        selection_and_drift(population_size);
        extinction();
        fixation(time_current, substitutions, non_syn_mutations, syn_mutations, burn_in,
            population_size);
    }

    u_long random_nuc_site() const {
        return std::uniform_int_distribution<u_long>(0, nbr_nucleotides - 1)(generator);
    }

    void mutation(NucleotideRateMatrix const &p, BinomialDistr &binomial_distrib,
        u_long &non_syn_mutations, u_long &syn_mutations, bool burn_in,
        double const &population_size) {
        TimeVar t_start = timeNow();

        // Randomly choose the number of mutations at this generation (for this exon)
        u_long binomial_draw = binomial_distrib(generator);

        // Early break if 0 mutation are drawn
        if (binomial_draw == 0) { return; }

        // Compute the distribution of haplotype frequency in the population
        std::discrete_distribution<u_long> haplotype_distri;
        if (haplotype_vector.size() != 1) {
            std::vector<u_long> nbr_copies(haplotype_vector.size(), 0);
            std::transform(haplotype_vector.begin(), haplotype_vector.end(), nbr_copies.begin(),
                [](Haplotype const &h) { return h.nbr_copies; });
            haplotype_distri =
                std::discrete_distribution<u_long>(nbr_copies.begin(), nbr_copies.end());
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
            std::array<char, 3> triplet_nuc = codonLexico.codon_to_triplet.at(codon_from);
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
            char codon_to = codonLexico.triplet_to_codon(
                triplet_nuc.at(0), triplet_nuc.at(1), triplet_nuc.at(2));
            if (codonLexico.codon_to_aa.at(codon_to) == 20) { continue; }

            Haplotype haplotype = haplotype_vector.at(hap_id);

            if (!is_synonymous(codon_from, codon_to)) {
                // Update the fitness of the new haplotype
                auto bk = codon_seq.at(codon_site);
                codon_seq[codon_site] = codon_from;
                haplotype.fitness_state->update(
                    codon_seq, haplotype.diff_sites, codon_site, codon_to, true, population_size);
                haplotype.sel_coeff =
                    fitness_state->selection_coefficient(*haplotype.fitness_state, burn_in);
                codon_seq[codon_site] = bk;
            }

            // Depending on whether the mutation is back to the reference sequence
            if (codon_to == codon_seq.at(codon_site)) {
                assert(haplotype.diff_sites.count(codon_site) > 0);
                haplotype.diff_sites.erase(codon_site);
            } else {
                haplotype.diff_sites[codon_site] = codon_to;
            }

            // Update the vector of haplotypes
            haplotype_vector.at(hap_id).nbr_copies--;
            haplotype.nbr_copies = 1;
            haplotype_vector.push_back(haplotype);

            // Track the created mutation
            if (burn_in) { continue; }
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
        std::vector<double> fitnesses(haplotype_vector.size(), 0);
        std::transform(haplotype_vector.begin(), haplotype_vector.end(), fitnesses.begin(),
            [](Haplotype const &h) { return (1.0 + h.sel_coeff) * h.nbr_copies; });
        double fit_tot = accumulate(fitnesses.begin(), fitnesses.end(), 0.0);

        // Random draws from the multinomial distribution
        for (std::size_t i_hap{0}; i_hap < haplotype_vector.size() - 1; i_hap++) {
            haplotype_vector.at(i_hap).nbr_copies = std::binomial_distribution<u_long>(
                children_tot, fitnesses.at(i_hap) / fit_tot)(generator);
            children_tot -= haplotype_vector.at(i_hap).nbr_copies;
            fit_tot -= fitnesses.at(i_hap);
        }
        assert(children_tot <= 2 * population_size);
        haplotype_vector.back().nbr_copies = children_tot;
        timer.selection += duration(timeNow() - t_start);
    }

    void extinction() {
        TimeVar t_start = timeNow();
        // Early break if only 1 haplotype
        if (haplotype_vector.size() == 1) { return; }

        // Sort the vector of haplotypes by number of copies
        if (!is_sorted(
                haplotype_vector.begin(), haplotype_vector.end(), Haplotype::GreaterThan{})) {
            sort(haplotype_vector.begin(), haplotype_vector.end(), Haplotype::GreaterThan{});
        }

        // Remove haplotypes with 0 copies
        if (haplotype_vector.back().nbr_copies == 0) {
            auto low_bound = lower_bound(
                haplotype_vector.begin(), haplotype_vector.end(), 0, Haplotype::GreaterThan{});
            if (low_bound != haplotype_vector.end()) {
                haplotype_vector.erase(low_bound, haplotype_vector.end());
            }
        }
        timer.extinction += duration(timeNow() - t_start);
    }

    void fixation(double time_current, std::vector<Substitution> &substitutions,
        u_long &non_syn_mutations, u_long &syn_mutations, bool burn_in,
        double const &population_size) {
        TimeVar t_start = timeNow();

        // ! Copy before the loop is necessary, since it can be updated during the loop
        auto diffs = haplotype_vector.begin()->diff_sites;

        // For each site (different to the reference) of the most common haplotype
        for (auto const &diff : diffs) {
            u_long site = diff.first;
            char codon_to = diff.second;
            bool polymorphic_site = false;
            std::set<char> poly_codons{codon_to};
            for (std::size_t hap_id{1}; hap_id < haplotype_vector.size(); hap_id++) {
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
                std::unordered_map<char, u_long> codon_copies{};
                for (char codon : poly_codons) { codon_copies[codon] = 0; }
                for (auto const &haplotype : haplotype_vector) {
                    codon_copies.at(haplotype.diff_sites.at(site)) += haplotype.nbr_copies;
                }
                codon_to = std::max_element(
                    codon_copies.begin(), codon_copies.end(), [](auto const &p1, auto const &p2) {
                        return p1.second < p2.second;
                    })->first;
            }
            if (!is_synonymous(codon_from, codon_to)) {
                fitness_state->update(codon_seq, site, codon_to, burn_in, population_size);
                for (auto &haplotype : haplotype_vector) {
                    haplotype.sel_coeff =
                        fitness_state->selection_coefficient(*fitness_state, true);
                }
            }
            codon_seq[site] = codon_to;
            // Update the vector of haplotypes
            for (auto &haplotype : haplotype_vector) {
                if (haplotype.diff_sites.at(site) == codon_to) {
                    haplotype.diff_sites.erase(site);
                } else {
                    assert(poly_codons.size() > 1);
                    assert(poly_codons.find(haplotype.diff_sites.at(site)) != poly_codons.end());
                }
            }

            // Track the created substitution
            if (burn_in) { continue; }
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
            auto sub = Substitution(
                codon_from, codon_to, time_current, t_between, non_syn_mutations, syn_mutations);
            non_syn_mutations = 0;
            syn_mutations = 0;
            sub.add_to_trace(tracer_substitutions);
            substitutions.push_back(sub);
        }
        timer.fixation += duration(timeNow() - t_start);
    }

    void sample_one_individual(double const &pop_size) {
        // Compute the distribution of haplotype frequency in the population
        std::vector<u_long> nbr_copies(haplotype_vector.size(), 0);
        std::transform(haplotype_vector.begin(), haplotype_vector.end(), nbr_copies.begin(),
            [](Haplotype const &h) { return h.nbr_copies; });
        u_long rand_hap =
            std::discrete_distribution<u_long>(nbr_copies.begin(), nbr_copies.end())(generator);

        for (auto const &diff : haplotype_vector[rand_hap].diff_sites) {
            codon_seq[diff.first] = diff.second;
        }
        fitness_state->update(codon_seq, pop_size);
        haplotype_vector.emplace_back(2 * pop_size, 0.0, fitness_state);
    }

    std::tuple<double, double> flow(NucleotideRateMatrix const &nuc_matrix) const {
        double non_syn_mut_flow{0.0}, syn_mut_flow{0.0};
        // For all site of the sequence.
        for (u_long site{0}; site < nbr_sites; site++) {
            // Codon original before substitution.
            char codon_from = codon_seq[site];

            // Array of neighbors of the original codon (codons differing by only 1 mutation).
            std::array<std::tuple<char, char, char>, 9> neighbors =
                codonLexico.codon_to_neighbors[codon_from];

            // For all possible neighbors.
            for (char neighbor{0}; neighbor < 9; neighbor++) {
                // Codon after mutation, Nucleotide original and Nucleotide after mutation.
                char codon_to{0}, n_from{0}, n_to{0};
                std::tie(codon_to, n_from, n_to) = neighbors[neighbor];

                // If the mutated amino-acid is a stop codon, the rate of fixation is 0.
                // Else, if the mutated and original amino-acids are non-synonymous, we compute the
                // rate of fixation. Note that, if the mutated and original amino-acids are
                // synonymous, the rate of fixation is 1.
                if (codonLexico.codon_to_aa[codon_to] != 20) {
                    if (codonLexico.codon_to_aa[codon_from] != codonLexico.codon_to_aa[codon_to]) {
                        non_syn_mut_flow += nuc_matrix.normalized_rate(n_from, n_to);
                    } else {
                        syn_mut_flow += nuc_matrix.normalized_rate(n_from, n_to);
                    }
                }
            }
        }
        return std::make_tuple(non_syn_mut_flow, syn_mut_flow);
    }
};

class Population {
  public:
    // TimeElapsed
    double time_from_root{0};
    // Exons
    std::vector<Exon> exons;

    u_long sample_size{0};
    LogMultivariate log_multivariate;
    double ornstein_uhlenbeck_sigma{};
    double ornstein_uhlenbeck_theta{};
    OrnsteinUhlenbeck ornstein_uhlenbeck;

    u_long population_size{0};
    double generation_time{0};
    NucleotideRateMatrix nuc_matrix;
    EMatrix const &transform_matrix;
    bool branch_wise{};

    std::vector<Substitution> substitutions{};
    u_long non_syn_mutations{0}, syn_mutations{0};
    PieceWiseMultivariate piecewise_multivariate{};

    mutable MapBinomialDistr binomial_distribs;

    Population(FitnessModel &seq_fitness, u_long sample_size, LogMultivariate &log_multi,
        NucleotideRateMatrix nucleotide_matrix, EMatrix const &transform_matrix, bool branch_wise,
        double ornstein_uhlenbeck_sigma, double ornstein_uhlenbeck_theta)
        : exons{},
          sample_size{sample_size},
          log_multivariate{log_multi},
          ornstein_uhlenbeck_sigma{ornstein_uhlenbeck_sigma},
          ornstein_uhlenbeck_theta{ornstein_uhlenbeck_theta},
          ornstein_uhlenbeck(ornstein_uhlenbeck_sigma, ornstein_uhlenbeck_theta, generator),
          population_size{log_multivariate.population_size()},
          generation_time{log_multivariate.generation_time()},
          nuc_matrix{std::move(nucleotide_matrix)},
          transform_matrix{transform_matrix},
          branch_wise{branch_wise} {
        u_long pos = 0;
        for (std::unique_ptr<FitnessState> &exon_seq_fitness : seq_fitness.fitness_states) {
            exons.emplace_back(exon_seq_fitness, pos, population_size, nuc_matrix);
            pos += exons.back().nbr_sites;
        }
        update_binomial_distribs();
    }

    void update_binomial_distribs() const {
        binomial_distribs.clear();
        for (auto const &exon : exons) {
            if (binomial_distribs.count(exon.nbr_sites) == 0) {
                binomial_distribs[exon.nbr_sites] =
                    std::binomial_distribution<u_long>(2 * population_size * exon.nbr_nucleotides,
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

    void set_from_aa_seq(std::string const &aa_seq) {
        assert(aa_seq.size() == nbr_sites());
        for (auto &exon : exons) {
            for (u_long site{0}; site < exon.nbr_sites; site++) {
                char aa_char = aa_seq.at(exon.position + site);
                char aa = codonLexico.aa_char_to_aa(aa_char);
                auto it =
                    std::find(codonLexico.codon_to_aa.begin(), codonLexico.codon_to_aa.end(), aa);
                assert(it != codonLexico.codon_to_aa.end());
                char codon_to = std::distance(codonLexico.codon_to_aa.begin(), it);
                exon.codon_seq[site] = codon_to;
            }
            exon.fitness_state->update(exon.codon_seq, population_size);
        }
    }

    bool check_consistency() const {
        for (auto const &exon : exons) {
            if (!exon.check_consistency(population_size)) {
                std::cerr << "The exon is not consistent" << std::endl;
                return false;
            }
        }
        for (std::size_t i = 1; i < substitutions.size(); ++i) {
            if (abs(substitutions.at(i).time_between -
                    (substitutions.at(i).time_event - substitutions.at(i - 1).time_event)) > 1e-6) {
                std::cerr << "The time between substitutions is not consistent" << std::endl;
                return false;
            }
            if (substitutions.at(i).time_between < 0.0) {
                std::cerr << "The substitutions are not ordered (in increasing time)" << std::endl;
                return false;
            }
        }
        return true;
    };

    EVector delta_log_multivariate(double distance) {
        TimeVar t_start = timeNow();
        EVector normal_vector = EVector::Zero(dimensions);
        for (int dim = 0; dim < dimensions; dim++) {
            normal_vector(dim) = normal_distrib(generator_brownian);
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
            std::cerr << "The population size (" << population_size
                      << ") became lower than sample size (" << sample_size << ")" << std::endl;
            population_size = 2 * sample_size - population_size;
            log_multivariate.set_population_size(population_size);
        }
        if (ornstein_uhlenbeck_sigma > 0) {
            ornstein_uhlenbeck.Next();
            population_size = std::max(
                static_cast<u_long>(population_size * ornstein_uhlenbeck.GetExpVal()), sample_size);
        }
        update_binomial_distribs();
        timer.correlation += duration(timeNow() - t_start);
    }

    void run_and_trace(std::string const &output_filename, Tree::NodeIndex node, Tree &tree,
        std::unique_ptr<Population> const &parent) {
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
                    time_current, substitutions, non_syn_mutations, syn_mutations, false);
            }
            assert(check_consistency());
            time_current += generation_time;
            time_from_root += generation_time;
        }
        node_trace(output_filename, node, tree, parent);
        if (branch_wise) { update_brownian(delta / 2); }
    }

    void burn_in(u_long burn_in_length) {
        std::cout << "Burn-in for " << burn_in_length << " generations." << std::endl;
        for (u_long gen{1}; gen <= burn_in_length; gen++) {
            for (auto &exon : exons) {
                assert(!exon.haplotype_vector.empty());
                exon.forward(nuc_matrix, population_size, binomial_distribs.at(exon.nbr_sites), 0,
                    substitutions, non_syn_mutations, syn_mutations, true);
            }
            assert(check_consistency());
        }
        clear_events();
        std::cout << "Burn-in completed." << std::endl;
    }

    double theoretical_theta() const { return 4 * population_size * nuc_matrix.mutation_rate; };

    double predicted_dn_dn0(NucleotideRateMatrix const &rates, u_long const &pop_size) const {
        double dn{0.}, dn0{0.};
        for (auto const &exon : exons) {
            double exon_dn{0}, exon_d0{0};
            std::tie(exon_dn, exon_d0) =
                exon.fitness_state->predicted_dn_dn0(exon.codon_seq, rates, pop_size, 0.0);
            dn += exon_dn;
            dn0 += exon_d0;
        }
        return dn / dn0;
    };

    double sequence_wise_predicted_dn_dn0(
        Population const &parent, NucleotideRateMatrix const &rates, u_long const &pop_size) const {
        double dn{0.}, dn0{0.};

        for (std::size_t i = 0; i < exons.size(); i++) {
            double exon_dn{0}, exon_dn0{0};
            std::tie(exon_dn, exon_dn0) =
                exons[i].fitness_state->flow_dn_dn0(parent.exons.at(i).codon_seq, rates, pop_size, 0.0);
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
            std::cerr << "There is no synonymous substitutions, dN/dS can't be computed!"
                      << std::endl;
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
            std::cerr << "There is no synonymous substitutions, dN/dS can't be computed!"
                      << std::endl;
            return .0;
        } else {
            return (dn * ds0) / (dn0 * ds);
        }
    }

    void node_trace(std::string const &output_filename, Tree::NodeIndex node, Tree &tree,
        std::unique_ptr<Population> const &parent) {
        TimeVar t_start = timeNow();

        std::string node_name = tree.node_name(node);

        tree.set_tag(node, "population_size", std::to_string(log_multivariate.population_size()));
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

        if (!tree.is_leaf(node)) {
            tracer_fossils.add("NodeName", node_name);
            double age = tree.max_distance_to_root() - time_from_root;
            tracer_fossils.add("Age", age);
            tracer_fossils.add("LowerBound", age * 0.9);
            tracer_fossils.add("UpperBound", age * 1.1);
            return;
        }

        // VCF file on the sample
        Polymorphism exome_poly(sample_size);

        std::string out;
        out += "##fileformat=VCFv4.0";
        out += "\n##source=SimuPoly";
        out += "\n##nodeName=" + node_name;
        out += "\n##sequenceSize=" + std::to_string(nbr_nucleotides());
        out += "\n##ploidyLevel=diploid";
        out += "\n##numberIndividuals=" + std::to_string(sample_size);
        out += "\n##numberGenotypes=" + std::to_string(2 * sample_size);
        out += "\n##reference=" + get_dna_str();
        out += "\n" + info;

        for (u_long indiv{1}; indiv <= sample_size; indiv++) {
            out += "\tId";
            out += std::to_string(indiv);
        }

        for (auto const &exon : exons) {
            Polymorphism exon_poly(sample_size);

            // Draw the sample of individuals
            std::vector<u_long> haplotypes_sample(2 * population_size, 0);
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
                std::map<std::tuple<char, char>, u_long> codon_from_to_copy{};
                for (auto const &i_hap : haplotypes_sample) {
                    if (exon.haplotype_vector.at(i_hap).diff_sites.count(site) > 0) {
                        char codon_from = exon.codon_seq.at(site);
                        char codon_to = exon.haplotype_vector.at(i_hap).diff_sites.at(site);
                        assert(codon_from != codon_to);
                        codon_from_to_copy[std::make_tuple(codon_from, codon_to)]++;
                    }
                }

                if (codon_from_to_copy.size() == 1) {
                    char codon_from = std::get<0>(codon_from_to_copy.begin()->first);
                    char codon_to = std::get<1>(codon_from_to_copy.begin()->first);
                    assert(codon_to != codon_from);
                    u_long alt_freq = codon_from_to_copy.begin()->second;
                    char nuc_pos{0};
                    while (nuc_pos < 3) {
                        if (codonLexico.codon_to_nuc(codon_from, nuc_pos) !=
                            codonLexico.codon_to_nuc(codon_to, nuc_pos)) {
                            break;
                        } else {
                            nuc_pos++;
                        }
                    }
                    assert(nuc_pos != 3);
                    std::string line{"\n"};
                    line += ".\t";
                    line += std::to_string(3 * (exon.position + site) + nuc_pos);
                    line += "\t.\t";
                    line += codonLexico.codon_to_nuc(codon_from, nuc_pos);
                    line += "\t";
                    line += codonLexico.codon_to_nuc(codon_to, nuc_pos);
                    line += "\t100\t";
                    if (alt_freq == 2 * sample_size) {
                        line += "s100";
                    } else if (alt_freq > sample_size) {
                        line += "s50";
                    } else {
                        line += "PASS";
                    }
                    line += "\tREFCODON=";
                    line += codonLexico.codon_string(codon_from);
                    line += ";ALTCODON=";
                    line += codonLexico.codon_string(codon_to);
                    line += ";REFAA=";
                    line += codonLexico.codon_aa_string(codon_from);
                    line += ";ALTAA=";
                    line += codonLexico.codon_aa_string(codon_to);
                    line += ";POSITION=";
                    line += std::to_string(nuc_pos);
                    line += ";ALTCOUNT=";
                    line += std::to_string(alt_freq);
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
                                nuc = codonLexico.codon_to_nuc(codon, nuc_pos);
                            } else {
                                nuc = codonLexico.codon_to_nuc(exon.codon_seq.at(site), nuc_pos);
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

            std::tie(exon_poly.non_syn_flow, exon_poly.syn_flow) = exon.flow(nuc_matrix);

            tracer_leaves.add("taxon_name", node_name);
            tracer_leaves.add("exon_id", exon.position);
            exon_poly.add_to_trace(tracer_leaves);

            exome_poly.add(exon_poly);
        }
        std::ofstream vcf;
        vcf.open(output_filename + "_" + node_name + ".vcf");
        vcf << out << std::endl;
        vcf.close();

        exome_poly.add_to_tree(tree, node);
        tree.set_tag(node, "theta_pred", d_to_string(theoretical_theta()));

        tracer_traits.add("TaxonName", node_name);
        tracer_traits.add("LogGenerationTime", log_multivariate.log_generation_time());
        tracer_traits.add("LogPopulationSize", log_multivariate.log_population_size());
        tracer_traits.add(
            "LogMutationRatePerGeneration", log_multivariate.log_mutation_rate_per_generation());

        // If the node is a leaf, output the DNA sequence and name.
        sample_one_individual(population_size);
        write_sequence(output_filename, node_name, get_dna_str());

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

    void sample_one_individual(double const &pop_size) {
        for (auto &exon : exons) { exon.sample_one_individual(pop_size); }
    }

    // Method returning the DNA std::string corresponding to the codon sequence.
    std::string get_dna_str() const {
        std::string dna_str{};
        dna_str.reserve(nbr_nucleotides());

        // For each site of the sequence.
        for (auto const &exon : exons) {
            for (u_long site{0}; site < exon.nbr_sites; site++) {
                // Assert there is no stop in the sequence.
                assert(codonLexico.codon_to_aa.at(exon.codon_seq.at(site)) != 20);

                // Translate the site to a triplet of DNA nucleotides
                std::array<char, 3> triplet =
                    codonLexico.codon_to_triplet.at(exon.codon_seq.at(site));
                for (char nuc_pos{0}; nuc_pos < 3; nuc_pos++) {
                    dna_str += codonLexico.nucleotides.at(triplet.at(nuc_pos));
                }
            }
        }
        return dna_str;  // return the DNA sequence as a std::string.
    }

    static void timer_cout() {
        std::cout << Exon::nbr_fixations << " fixations" << std::endl;
        std::cout << Exon::nbr_mutations << " mutations" << std::endl;
        double total_time = timer.mutation + timer.selection + timer.extinction + timer.fixation +
                            timer.exportation + timer.correlation;
        std::cout << std::setprecision(3) << total_time / 1e9 << "s total time" << std::endl;
        std::cout << 100 * timer.mutation / total_time
                  << "% of time spent in calculating mutation (" << timer.mutation / 1e9 << "s)"
                  << std::endl;
        std::cout << 100 * timer.selection / total_time
                  << "% of time spent in calculating selection (" << timer.selection / 1e9 << "s)"
                  << std::endl;
        std::cout << 100 * timer.extinction / total_time
                  << "% of time spent in calculating extinction (" << timer.extinction / 1e9 << "s)"
                  << std::endl;
        std::cout << 100 * timer.fixation / total_time
                  << "% of time spent in calculating fixation (" << timer.fixation / 1e9 << "s)"
                  << std::endl;
        std::cout << 100 * timer.exportation / total_time << "% of time spent in exporting vcf ("
                  << timer.exportation / 1e9 << "s)" << std::endl;
        std::cout << 100 * timer.correlation / total_time << "% of time spent in the correlation ("
                  << timer.correlation / 1e9 << "s)" << std::endl;
    }
};

// Initialize static variables
u_long Exon::nbr_mutations = 0, Exon::nbr_fixations = 0;

class Process {
  private:
    static double years_computed;
    Tree &tree;
    std::vector<std::unique_ptr<Population>> populations;  // Vector of sequence DNA.

  public:
    // Constructor
    Process(Tree &intree, std::unique_ptr<Population> &root_pop) : tree{intree}, populations() {
        populations.resize(tree.nb_nodes());
        populations[tree.root()] = std::move(root_pop);
    }

    void run(std::string &output_filename) {
        populations.at(tree.root())->node_trace(output_filename, tree.root(), tree, nullptr);
        run_recursive(tree.root(), output_filename);
        std::ofstream nhx;
        nhx.open(output_filename + ".nhx");
        nhx << tree.as_string() << std::endl;
        nhx.close();
    }


    // Recursively iterate through the subtree.
    void run_recursive(Tree::NodeIndex node, std::string const &output_filename) {
        // Substitutions of the DNA sequence is generated.
        populations.at(node)->clear_events();

        if (!tree.is_root(node)) {
            populations.at(node)->run_and_trace(
                output_filename, node, tree, populations.at(tree.parent(node)));

            years_computed += tree.node_length(node);
            std::cout << years_computed << " years computed in total ("
                      << static_cast<int>(100 * years_computed / tree.total_length())
                      << "%) at node " << tree.node_name(node) << " ("
                      << static_cast<int>(100 * tree.node_length(node) / tree.total_length())
                      << "%)." << std::endl;
        }

        // Iterate through the direct children.
        for (auto &child : tree.children(node)) {
            populations.at(child) = std::make_unique<Population>(*populations.at(node));
            run_recursive(child, output_filename);
        }
    }
};

// Initialize static variables
double Process::years_computed = 0.0;
