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

struct Change {
    u_long site;
    char codon_from{-1};
    char codon_to{-1};
};

bool operator<(const Change &lhs, const Change &rhs) { return lhs.site < rhs.site; }

bool operator==(const Change &lhs, const Change &rhs) { return lhs.site == rhs.site; }

bool is_synonymous(const Change &change) {
    return Codon::codon_to_aa_array[change.codon_from] == Codon::codon_to_aa_array[change.codon_to];
};

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
        trace.add("Branch_dS", this->ds());
        trace.add("Branch_dN", this->dn());
        trace.add("Branch_dNdS", this->dnds());
        trace.add("Branch_#mutations", this->mutations());
        trace.add("Branch_#fixations", this->fixations());
        trace.add("Branch_#non_synonymous_mutations", this->non_syn_mut);
        trace.add("Branch_#synonymous_mutations", this->syn_mut);
        trace.add("Branch_#non_synonymous_fixations", this->non_syn_fix);
        trace.add("Branch_#synonymous_fixations", this->syn_fix);
    }

    void add_to_tree(Tree &tree, Tree::NodeIndex node) const {
        tree.set_tag(node, "Branch_dS", d_to_string(this->ds()));
        tree.set_tag(node, "Branch_dN", d_to_string(this->dn()));
        tree.set_tag(node, "Branch_dNdS", d_to_string(this->dnds()));
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
    char codon_from{-1};
    char codon_to{-1};
    double time_event{0};
    double time_between_event{0};
    u_long non_syn_mut_flow{0};
    u_long syn_mut_flow{0};

    explicit Substitution() = default;

    bool is_synonymous() const {
        return (codon_from != codon_to and
                Codon::codon_to_aa_array[codon_from] == Codon::codon_to_aa_array[codon_to]);
    }

    bool is_non_synonymous() const {
        return (codon_from != codon_to and
                Codon::codon_to_aa_array[codon_from] != Codon::codon_to_aa_array[codon_to]);
    }

    bool is_dummy() const { return codon_from == codon_to; }
};

class Polymorphism {
  private:
    u_long sample_size;
    u_long seq_size;
    double inv_harmonic_sum;
    double inv_choice;

  public:
    u_long pairwise_syn{0};
    u_long pairwise_non_syn{0};
    u_long non_syn_nbr{0};
    u_long syn_nbr{0};
    u_long complex_sites{0};

    explicit Polymorphism(u_long sample_size_, u_long seq_size_)
        : sample_size{sample_size_}, seq_size{seq_size_} {
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
    }

    void add_to_trace(Trace &trace) const {
        trace.add("theta_watterson_non_syn", this->theta_watterson_non_syn());
        trace.add("theta_watterson_syn", this->theta_watterson_syn());
        trace.add("theta_watterson", this->theta_watterson());
        trace.add("theta_pairwise_non_syn", this->theta_pairwise_non_syn());
        trace.add("theta_pairwise_syn", this->theta_pairwise_syn());
        trace.add("theta_pairwise", this->theta_pairwise());
    }

    void add_to_tree(Tree &tree, Tree::NodeIndex node) const {
        tree.set_tag(node, "theta_watterson_non_syn", d_to_string(this->theta_watterson_non_syn()));
        tree.set_tag(node, "theta_watterson_syn", d_to_string(this->theta_watterson_syn()));
        tree.set_tag(node, "theta_watterson", d_to_string(this->theta_watterson()));
        tree.set_tag(node, "theta_pairwise_non_syn", d_to_string(this->theta_pairwise_non_syn()));
        tree.set_tag(node, "theta_pairwise_syn", d_to_string(this->theta_pairwise_syn()));
        tree.set_tag(node, "theta_pairwise", d_to_string(this->theta_pairwise()));
    }

  private:
    double theta_watterson_non_syn() const { return non_syn_nbr * inv_harmonic_sum / seq_size; }

    double theta_watterson_syn() const { return syn_nbr * inv_harmonic_sum / seq_size; }

    double theta_watterson() const { return (non_syn_nbr + syn_nbr) * inv_harmonic_sum / seq_size; }

    double theta_pairwise_non_syn() const { return pairwise_non_syn * inv_choice / seq_size; }

    double theta_pairwise_syn() const { return pairwise_syn * inv_choice / seq_size; }

    double theta_pairwise() const {
        return (pairwise_non_syn + pairwise_syn) * inv_choice / seq_size;
    }
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

class Haplotype {
  public:
    // The nbr_copies of sites in the sequence (each position is a codon, thus the DNA sequence is 3
    // times greater).
    u_long nbr_copies{0};
    double fitness{0.0};
    set<Change> set_change;

    explicit Haplotype() = default;

    // Constructor of Reference_seq.
    // size: the size of the DNA sequence.
    Haplotype(u_long const nbr_copies, double const fitness)
        : nbr_copies{nbr_copies}, fitness{fitness} {};

    bool check_consistency(u_long nbr_sites) const {
        for (auto &change : set_change) {
            if (Codon::codon_to_aa_array[change.codon_to] == 20) { return false; }
            if (Codon::codon_to_aa_array[change.codon_from] == 20) { return false; }
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
    u_long const &population_size;
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

    explicit Exon() = default;

    // Constructor
    explicit Exon(vector<array<double, 20>> const &fitness_profiles, const u_long &position,
        u_long &population_size, NucleotideRateMatrix const &nuc_matrix)
        : population_size{population_size},
          position{position},
          aa_fitness_profiles{fitness_profiles},
          nbr_sites{u_long(fitness_profiles.size())},
          nbr_nucleotides{u_long(3 * fitness_profiles.size())},
          codon_seq(fitness_profiles.size(), 0),
          haplotype_vector{} {
        // Draw codon from codon frequencies
        for (u_long site{0}; site < nbr_sites; site++) {
            array<double, 64> codon_freqs =
                codon_frequencies(aa_fitness_profiles[site], nuc_matrix, 4 * population_size);
            discrete_distribution<char> freq_codon_distr(codon_freqs.begin(), codon_freqs.end());
            codon_seq[site] = freq_codon_distr(generator);
        }

        haplotype_vector.emplace_back(Haplotype(2 * population_size, 0.0));
        assert(check_consistency());
    }

    bool check_consistency() const {
        if (haplotype_vector.size() > 2 * population_size) { return false; }
        if (haplotype_vector.empty()) { return false; }

        u_long nbr_copies{0};
        for (auto &haplotype : haplotype_vector) {
            if (!haplotype.check_consistency(nbr_sites)) { return false; }
            nbr_copies += haplotype.nbr_copies;
        }
        return nbr_copies == 2 * population_size;
    }

    void forward(
        NucleotideRateMatrix const &nuc_matrix, double time_current, vector<Substitution> &subs) {
        mutation(nuc_matrix, subs);
        selection_and_drift();
        extinction();
        fixation(time_current, subs);
    }

    discrete_distribution<u_long> haplotype_freq_distr() const {
        vector<u_long> nbr_copies(haplotype_vector.size(), 0);
        std::transform(haplotype_vector.begin(), haplotype_vector.end(), nbr_copies.begin(),
            [](const Haplotype &h) { return h.nbr_copies; });

        discrete_distribution<u_long> haplotype_distr(nbr_copies.begin(), nbr_copies.end());
        return haplotype_distr;
    }

    void mutation(NucleotideRateMatrix const &p, vector<Substitution> &subs) {
        TimeVar t_start = timeNow();

        binomial_distribution<u_long> binomial_distr(
            2 * population_size * nbr_nucleotides, p.max_sum_mutation_rates);
        u_long binomial_draw = binomial_distr(generator);

        if (binomial_draw > 0) {
            set<tuple<u_long, u_long, u_long>> coordinates_set{};
            u_long nbr_draws{0};

            discrete_distribution<u_long> haplotype_distr{};
            if (haplotype_vector.size() > 1) { haplotype_distr = haplotype_freq_distr(); }

            while (nbr_draws < binomial_draw) {
                u_long haplotype_draw = 0;

                if (haplotype_vector.size() > 1) { haplotype_draw = haplotype_distr(generator); }

                uniform_int_distribution<u_long> copy_and_site_distr(
                    0, haplotype_vector[haplotype_draw].nbr_copies * nbr_nucleotides - 1);
                u_long copy_and_site_draw = copy_and_site_distr(generator);
                u_long nuc_site = copy_and_site_draw % nbr_nucleotides;
                u_long copy = copy_and_site_draw / nbr_nucleotides;

                auto coordinate = make_tuple(haplotype_draw, copy, nuc_site);
                if (coordinates_set.count(coordinate) == 0) {
                    coordinates_set.insert(coordinate);
                    nbr_draws++;
                }
            }

            auto iter = coordinates_set.begin();
            while (iter != coordinates_set.end()) {
                u_long i_hap = get<0>(*iter);
                u_long copy_id = get<1>(*iter);
                bool at_least_one_mutation{false};
                Haplotype haplotype{};
                while (true) {
                    u_long site = get<2>(*iter);
                    u_long codon_site = site / 3;
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
                                events.syn_mut++;
                                subs.back().syn_mut_flow++;
                            } else {
                                events.non_syn_mut++;
                                subs.back().non_syn_mut_flow++;
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
        timer.mutation += duration(timeNow() - t_start);
    }

    void selection_and_drift() {
        TimeVar t_start = timeNow();
        if (haplotype_vector.size() > 1) {
            vector<double> fitnesses(haplotype_vector.size(), 0);

            std::transform(haplotype_vector.begin(), haplotype_vector.end(), fitnesses.begin(),
                [](const Haplotype &h) { return (1.0 + h.fitness) * h.nbr_copies; });

            double fit_tot = accumulate(fitnesses.begin(), fitnesses.end(), 0.0);
            u_long children_tot = 2 * population_size;

            for (size_t i_hap{0}; i_hap < haplotype_vector.size() - 1; i_hap++) {
                haplotype_vector[i_hap].nbr_copies = binomial_distribution<u_long>(
                    children_tot, fitnesses[i_hap] / fit_tot)(generator);
                children_tot -= haplotype_vector[i_hap].nbr_copies;
                fit_tot -= fitnesses[i_hap];
            }
            assert(children_tot <= 2 * population_size);
            assert(children_tot >= 0);
            haplotype_vector[haplotype_vector.size() - 1].nbr_copies = children_tot;
        }
        timer.selection += duration(timeNow() - t_start);
    }

    void extinction() {
        TimeVar t_start = timeNow();

        // remove haplotypes with 0 copies
        if (haplotype_vector.size() > 1) {
            if (!is_sorted(
                    haplotype_vector.begin(), haplotype_vector.end(), Haplotype::GreaterThan)) {
                sort(haplotype_vector.begin(), haplotype_vector.end(), Haplotype::GreaterThan);
            }
            if (haplotype_vector.back().nbr_copies == 0) {
                auto low_bound = lower_bound(
                    haplotype_vector.begin(), haplotype_vector.end(), 0, Haplotype::GreaterThan);
                if (low_bound != haplotype_vector.end()) {
                    haplotype_vector.erase(low_bound, haplotype_vector.end());
                }
            }
        }

        timer.extinction += duration(timeNow() - t_start);
    }

    void fixation(double time_current, vector<Substitution> &subs) {
        TimeVar t_start = timeNow();
        set<Change> current_change_set = haplotype_vector[0].set_change;
        u_long i_hap{1};
        while (!current_change_set.empty() and i_hap < haplotype_vector.size()) {
            set<Change> intersect;
            set_intersection(current_change_set.begin(), current_change_set.end(),
                haplotype_vector[i_hap].set_change.begin(),
                haplotype_vector[i_hap].set_change.end(), inserter(intersect, intersect.begin()));
            current_change_set = intersect;
            i_hap++;
        }
        for (auto const &change : current_change_set) {
            size_t count = 0;
            for (auto const &haplotype : haplotype_vector) {
                count += haplotype.set_change.count(change);
            }

            if (count == haplotype_vector.size()) {
                for (auto &haplotype : haplotype_vector) { haplotype.set_change.erase(change); }
                assert(change.site - position >= 0);
                assert(change.site - position < codon_seq.size());
                codon_seq[change.site - position] = change.codon_to;
                if (is_synonymous(change)) {
                    events.syn_fix++;
                } else {
                    events.non_syn_fix++;
                }
                nbr_fixations++;
                subs.back().codon_from = change.codon_from;
                subs.back().codon_to = change.codon_to;
                subs.back().time_event = time_current;
                if (subs.size() > 1) {
                    subs.back().time_between_event = time_current - subs.rbegin()[1].time_event;
                } else {
                    subs.back().time_between_event = time_current;
                }
                subs.emplace_back(Substitution());
            }
        }
        timer.fixation += duration(timeNow() - t_start);
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

    vector<Substitution> substitutions;

    explicit Population() = default;

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

            exons.emplace_back(Exon(exon_profiles, begin_exon, population_size, nuc_matrix));
        }
        assert(nbr_sites() == fitness_profiles.size());
        cout << exons.size() << " exons created." << endl;
        if (dv.rem != 0) { cout << "Last exon is " << dv.rem << " sites long." << endl; }
        substitutions.emplace_back(Substitution());
    }

    u_long nbr_sites() const {
        u_long sites = 0;
        for (auto const &exon : exons) { sites += exon.nbr_sites; }
        return sites;
    }

    u_long nbr_nucleotides() const { return 3 * nbr_sites(); }

    bool check_consistency() const {
        for (auto const &exon : exons) {
            if (!exon.check_consistency()) { return false; }
        }
        for (size_t i = 1; i < substitutions.size(); ++i) {
            if (abs(substitutions[i].time_between_event -
                    (substitutions[i].time_event - substitutions[i - 1].time_event)) > 1e-6) {
                return false;
            }
        }
        if (substitutions.empty()) { return false; }
        return count_if(substitutions.begin(), substitutions.end(),
                   [](Substitution const &s) { return s.is_dummy(); }) == 1;
    };

    EVector delta_log_multivariate(double distance) {
        TimeVar t_start = timeNow();
        EVector normal_vector = EVector::Zero(log_multivariate.dimensions);
        for (int dim = 0; dim < log_multivariate.dimensions; dim++) {
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
        nuc_matrix.set_mutation_rate(log_multivariate.mu());

        if (population_size < sample_size) {
            cerr << "The population size (" << population_size
                 << ") became lower than sample size (" << sample_size << ")" << endl;
            population_size = 2 * sample_size - population_size;
            log_multivariate.set_population_size(population_size);
        };
        timer.correlation += duration(timeNow() - t_start);
    }

    void run_forward(double t_max, Tree const &tree) {
        double time_current = 0.0;
        EVector delta;

        if (branch_wise) {
            delta = delta_log_multivariate(t_max / tree.max_distance_to_root());
            update_brownian(delta / 2);
        }
        while (time_current < t_max) {
            if (!branch_wise) {
                delta = delta_log_multivariate(generation_time / tree.max_distance_to_root());
                update_brownian(delta);
            }
            for (auto &exon : exons) {
                assert(!exon.haplotype_vector.empty());
                exon.forward(nuc_matrix, time_current, substitutions);
            }
            time_current += generation_time;
            time_from_root += generation_time;
        }
        substitutions.back().time_event = t_max;
        if (substitutions.size() > 1) {
            substitutions.back().time_between_event = t_max - substitutions.rbegin()[1].time_event;
        } else {
            substitutions.back().time_between_event = t_max;
        }
        if (branch_wise) { update_brownian(delta / 2); }
        assert(check_consistency());
    }

    void burn_in(u_long burn_in_length) {
        cout << "Burn-in for " << burn_in_length << " generations." << endl;
        for (u_long gen{1}; gen <= burn_in_length; gen++) {
            for (auto &exon : exons) {
                assert(!exon.haplotype_vector.empty());
                exon.forward(nuc_matrix, 0, substitutions);
            }
            assert(check_consistency());
        }
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
                exons[i].aa_fitness_profiles, parent.exons[i].codon_seq, rates, 4 * pop_size);
            dn += exon_dn;
            dn0 += exon_dn0;
        }
        return dn / dn0;
    };

    double count_based_dn_dn0() const {
        double dn{0}, dn0{0}, time_total{0};
        for (auto const &substitution : substitutions) {
            if (substitution.is_non_synonymous()) { dn++; }
            time_total += substitution.time_between_event;
            dn0 += substitution.non_syn_mut_flow * substitution.time_between_event;
        }
        return time_total * dn / dn0;
    }

    double event_based_dn_dn0() const {
        double dn{0}, dn0{0};
        for (auto const &substitution : substitutions) {
            if (substitution.is_non_synonymous()) { dn++; }
            dn0 += substitution.non_syn_mut_flow;
        }
        return dn / dn0;
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
            dn0 += substitution.non_syn_mut_flow * substitution.time_between_event;
            ds0 += substitution.syn_mut_flow * substitution.time_between_event;
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

    void node_trace(string const &output_filename, Tree::NodeIndex node, Tree &tree,
        Population const *parent) const {
        TimeVar t_start = timeNow();

        string node_name = tree.node_name(node);

        tree.set_tag(node, "population_size", to_string(population_size));
        tree.set_tag(node, "generation_time_in_year", d_to_string(generation_time));
        tree.set_tag(node, "mutation_rate_per_generation", d_to_string(nuc_matrix.mutation_rate));

        if (!tree.is_root(node)) {
            auto sum_events = events();
            sum_events.add_to_tree(tree, node);
            assert(parent != nullptr);
            LogMultivariate half = log_multivariate;
            half += parent->log_multivariate;
            half /= 2;
            NucleotideRateMatrix avg_nuc_matrix = nuc_matrix;
            avg_nuc_matrix.set_mutation_rate(half.mu());

            tree.set_tag(node, "Branch_population_size", to_string(half.population_size()));
            tree.set_tag(
                node, "Branch_generation_time_in_year", d_to_string(half.generation_time()));
            tree.set_tag(node, "Branch_mutation_rate_per_generation", d_to_string(half.mu()));
            tree.set_tag(node, "Branch_LogNe", d_to_string(log10(half.population_size())));
            tree.set_tag(node, "Branch_dNdN0_predicted",
                d_to_string(predicted_dn_dn0(avg_nuc_matrix, half.population_size())));
            tree.set_tag(node, "Branch_dNdN0_sequence_wise_predicted",
                d_to_string(sequence_wise_predicted_dn_dn0(
                    *parent, avg_nuc_matrix, half.population_size())));
            tree.set_tag(node, "Branch_dNdN0_event_based", d_to_string(event_based_dn_dn0()));
            tree.set_tag(node, "Branch_dNdN0_count_based", d_to_string(count_based_dn_dn0()));
            tree.set_tag(node, "Branch_dNdS_event_based", d_to_string(event_based_dn_ds()));
            tree.set_tag(node, "Branch_dNdS_count_based", d_to_string(count_based_dn_ds()));

            for (auto const &exon : exons) {
                tracer_nodes.add("taxon_name", node_name);
                tracer_nodes.add("exon_id", exon.position);
                exon.events.add_to_trace(tracer_nodes);
            }
        }

        if (tree.is_leaf(node)) {
            // If the node is a leaf, output the DNA sequence and name.
            write_sequence(output_filename, node_name, this->get_dna_str());

            double at_sites_obs{0};

            // VCF file on the sample
            Polymorphism exome_poly(sample_size, nbr_nucleotides());

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
                Polymorphism exon_poly(sample_size, exon.nbr_nucleotides);

                // Draw the sample of individuals
                vector<u_long> haplotypes_sample(2 * exon.population_size, 0);
                u_long sum_copies = 0;
                for (u_long i_hap{0}; i_hap < exon.haplotype_vector.size(); i_hap++) {
                    for (u_long copy_id{0}; copy_id < exon.haplotype_vector[i_hap].nbr_copies;
                         copy_id++) {
                        assert(sum_copies + copy_id < 2 * exon.population_size);
                        haplotypes_sample[sum_copies + copy_id] = i_hap;
                    }
                    sum_copies += exon.haplotype_vector[i_hap].nbr_copies;
                }
                shuffle(haplotypes_sample.begin(), haplotypes_sample.end(), generator);
                haplotypes_sample.resize(2 * sample_size);

                /*
                // DN/DS computed using one individual of the sample
                uniform_int_distribution<u_long> chosen_distr(0, 2 * sample_size - 1);
                u_long chosen = chosen_distr(generator);
                for (auto const &change :
                    exon.haplotype_vector[haplotypes_sample[chosen]].set_change) {
                    if (is_synonymous(change)) {
                        table.syn_fix++;
                    } else {
                        table.non_syn_fix++;
                    }
                }
                */

                for (u_long site{0}; site < exon.nbr_sites; site++) {
                    map<tuple<char, char>, u_long> codon_from_to_copy{};
                    for (auto const &i_hap : haplotypes_sample) {
                        char codon_to = exon.codon_seq[site];

                        auto it = exon.haplotype_vector[i_hap].set_change.find(
                            Change{exon.position + site});
                        if (it != exon.haplotype_vector[i_hap].set_change.end()) {
                            assert(it->site == site + exon.position);
                            char codon_from = it->codon_from;
                            codon_to = it->codon_to;
                            if (codon_to != codon_from) {
                                codon_from_to_copy[make_tuple(codon_from, codon_to)]++;
                            }
                        }

                        array<char, 3> triplet = Codon::codon_to_triplet_array[codon_to];
                        for (char position{0}; position < 3; position++) {
                            if (triplet[position] == 0 or triplet[position] == 3) {
                                at_sites_obs++;
                            }
                        }
                    }

                    if (codon_from_to_copy.size() == 1) {
                        char codon_from = get<0>(codon_from_to_copy.begin()->first);
                        char codon_to = get<1>(codon_from_to_copy.begin()->first);
                        if (codon_to != codon_from) {
                            u_long alt_freq = codon_from_to_copy.begin()->second;
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
                            line += to_string(3 * (exon.position + site) + position);
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
                                    auto it = exon.haplotype_vector[i_hap].set_change.find(
                                        Change{exon.position + site});
                                    char nuc{0};
                                    if (it != exon.haplotype_vector[i_hap].set_change.end()) {
                                        assert(it->site == exon.position + site);
                                        nuc = Codon::codon_to_nuc(it->codon_to, position);
                                    } else {
                                        nuc = Codon::codon_to_nuc(exon.codon_seq[site], position);
                                    }
                                    line += nuc;
                                    if (ploidy == 0) { line += "|"; }
                                }
                            }
                            out += line;
                        }
                    } else if (codon_from_to_copy.size() > 1) {
                        exon_poly.complex_sites++;
                    }
                }

                // Theta pairwise computed on the sample
                for (u_long i{0}; i < 2 * sample_size; i++) {
                    for (u_long j{i + 1}; j < 2 * sample_size; j++) {
                        if (i != j) {
                            u_long hap_i = haplotypes_sample[i];
                            u_long hap_j = haplotypes_sample[j];
                            auto it_first = exon.haplotype_vector[hap_i].set_change.begin();
                            auto end_first = exon.haplotype_vector[hap_i].set_change.end();
                            auto it_second = exon.haplotype_vector[hap_j].set_change.begin();
                            auto end_second = exon.haplotype_vector[hap_j].set_change.end();
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
                                    exon_poly.pairwise_syn++;
                                } else {
                                    exon_poly.pairwise_non_syn++;
                                }
                            }
                        }
                    }
                }


                /*
                exon_poly.mean_fitness +=
                accumulate(exon.haplotype_vector.begin(), exon.haplotype_vector.end(), 0.0,
                    [](double acc, Haplotype const &h) {
                        return acc + h.nbr_copies * h.fitness;
                    }) /
                exon.population_size;
                */

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
            tree.set_tag(node, "theta_pairwise_pred", d_to_string(theoretical_theta()));
            tree.set_tag(node, "theta_watterson_pred", d_to_string(theoretical_theta()));
        }

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
        substitutions.emplace_back(Substitution());
    }

    // Method returning the DNA string corresponding to the codon sequence.
    string get_dna_str() const {
        string dna_str{};
        dna_str.reserve(nbr_nucleotides());

        // For each site of the sequence.
        for (auto const &exon : exons) {
            for (u_long site{0}; site < exon.nbr_sites; site++) {
                // Assert there is no stop in the sequence.
                assert(Codon::codon_to_aa_array[exon.codon_seq[site]] != 20);

                // Translate the site to a triplet of DNA nucleotides
                array<char, 3> triplet = Codon::codon_to_triplet_array[exon.codon_seq[site]];
                for (char position{0}; position < 3; position++) {
                    dna_str += Codon::nucleotides[triplet[position]];
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
        cout << 100 * timer.exportation / total_time << "% of time spent in exportationing vcf ("
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
        populations[tree.root()]->node_trace(output_filename, tree.root(), tree, nullptr);
        run_recursive(tree.root(), output_filename);
        ofstream nhx;
        nhx.open(output_filename + ".nhx");
        nhx << tree.as_string() << endl;
        nhx.close();
    }


    // Recursively iterate through the subtree.
    void run_recursive(Tree::NodeIndex node, string const &output_filename) {
        // Substitutions of the DNA sequence is generated.

        if (!tree.is_root(node)) {
            populations[node]->run_forward(tree.node_length(node), tree);

            years_computed += tree.node_length(node);
            cout << years_computed << " years computed ("
                 << static_cast<int>(100 * years_computed / tree.total_length()) << "%)." << endl;

            populations[node]->node_trace(
                output_filename, node, tree, populations[tree.parent(node)]);
        }

        // Iterate through the direct children.
        for (auto &child : tree.children(node)) {
            populations[child] = new Population(*populations[node]);
            populations[child]->clear_events();
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
        "n", "pop_size", "Population size (at the root)", false, 500, "u_long", cmd};
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
    string correlation_path{args.correlation_path.getValue()};
    double mu{args.mu.getValue()};
    assert(mu > 0.0);
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
    NucleotideRateMatrix nuc_matrix(nuc_matrix_path, mu, true);

    LogMultivariate log_multivariate(pop_size, generation_time, mu);
    CorrelationMatrix correlation_matrix(correlation_path);

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
    parameters.add("mutation_rate_per_generation", mu);
    nuc_matrix.add_to_trace(parameters);
    parameters.add("generation_time_in_year", generation_time);
    parameters.add("#generations_burn_in", burn_in);
    parameters.add("population_size", pop_size);
    parameters.add("sample_size", sample_size);
    parameters.add("exon_size", exon_size);
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

    Population::timer_cout();

    cout << "Simulation computed." << endl;
    cout << nbr_sites * 3 * (mu / generation_time) * tree.total_length()
         << " expected substitutions." << endl;
    cout << Exon::nbr_fixations << " simulated substitutions." << endl;
    cout << "Statistics summarized in: " << output_path + ".tsv" << endl;
    cout << "Fasta file in: " << output_path + ".fasta" << endl;
    cout << "Alignment (.ali) file in: " << output_path + ".ali" << endl;
    return 0;
}