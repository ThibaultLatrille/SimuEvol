#pragma once

#include <queue>
#include "clonable_unique_ptr.hpp"
#include "fitness.hpp"
#include "random.hpp"

std::unordered_map<std::string, double> pfix_map{};

class DFE {
  public:
    u_long deleterious{0};
    u_long advantageous{0};
    double sub_deleterious{0};
    double sub_advantageous{0};
    double sel_coeff_subs{0};

    DFE() = default;

    void operator+=(const DFE &other) {
        this->deleterious += other.deleterious;
        this->advantageous += other.advantageous;
        this->sub_deleterious += other.sub_deleterious;
        this->sub_advantageous += other.sub_advantageous;
        this->sel_coeff_subs += other.sel_coeff_subs;
    }

    void add(Trace &trace) const {
        trace.add("advantageous_over_deleterious",
            static_cast<double>(advantageous) / static_cast<double>(deleterious));
        trace.add("sub_advantageous_minus_deleterious", sub_advantageous - sub_deleterious);
        trace.add("sub_advantageous_over_deleterious", sub_advantageous / sub_deleterious);
        trace.add("selection_coefficient_substitutions", sel_coeff_subs);
    }
};

class Substitution {
  public:
    u_long site;
    char codon_from;
    char codon_to;
    char n_from;
    char n_to;
    double time_event;
    double time_between_event;
    double non_syn_sub_flow;
    double non_syn_mut_flow;
    double non_syn_sub_flow_WS;
    double non_syn_mut_flow_WS;
    double non_syn_sub_flow_SW;
    double non_syn_mut_flow_SW;
    double syn_mut_flow;
    double dndn0;
    DFE dfe;

    Substitution(double time_event, double time_between_event, double non_syn_sub_flow,
        double non_syn_mut_flow, double non_syn_sub_flow_WS, double non_syn_mut_flow_WS,
        double non_syn_sub_flow_SW, double non_syn_mut_flow_SW, double syn_mut_flow, double dndn0,
        char codon_from = -1, char codon_to = -1, char n_from = -1, char n_to = -1, u_long site = 0,
        DFE const &dfe = {})
        : site{site},
          codon_from{codon_from},
          codon_to{codon_to},
          n_from{n_from},
          n_to{n_to},
          time_event{time_event},
          time_between_event{time_between_event},
          non_syn_sub_flow{non_syn_sub_flow},
          non_syn_mut_flow{non_syn_mut_flow},
          non_syn_sub_flow_WS{non_syn_sub_flow_WS},
          non_syn_mut_flow_WS{non_syn_mut_flow_WS},
          non_syn_sub_flow_SW{non_syn_sub_flow_SW},
          non_syn_mut_flow_SW{non_syn_mut_flow_SW},
          syn_mut_flow{syn_mut_flow},
          dndn0{dndn0},
          dfe{dfe} {}

    bool is_dummy() const { return codon_from == codon_to; }

    bool is_synonymous() const {
        return (!is_dummy() and
                codonLexico.codon_to_aa[codon_from] == codonLexico.codon_to_aa[codon_to]);
    }

    bool is_non_synonymous() const {
        return (!is_dummy() and
                codonLexico.codon_to_aa[codon_from] != codonLexico.codon_to_aa[codon_to]);
    }
};

// Class representing DNA sequence.
class Exon {
  public:
    // The number of sites in the exon (each position is a codon, thus the DNA sequence is 3
    // times greater).
    u_long nbr_sites;
    u_long position;

    // The sequence of codons.
    std::vector<char> codon_seq;

    // The fitness landscape (cannot be modified), which map genotype->phenotype->fitness.
    ClonableUniquePtr<FitnessState> fitness_state;

    std::queue<Substitution> substitutions{};

    // Constructor of Exon.
    // size: the size of the DNA sequence.
    Exon(std::unique_ptr<FitnessState> &f_state, u_long const &position, double const &pop_size)
        : nbr_sites{f_state->nbr_sites()},
          position{position},
          codon_seq(nbr_sites, 0),
          fitness_state{std::move(f_state)} {
        assert(substitutions.empty());
        fitness_state->update(codon_seq, pop_size);
    }

    bool operator==(Exon const &other) const { return *fitness_state == *other.fitness_state; };

    bool operator!=(Exon const &other) const { return !(*this == other); };

    // Method computing the next substitution event to occur, and the time for it to happen.
    // This method takes a time as input. If there is enough time given, the substitution event is
    // computed. This method also returns the time (given as input) decremented by the time needed
    // for the substitution event to occur. If there is not enough time given, no substitution event
    // is computed and the method returns 0. time_left: The time available for a substitution event
    // to occur.
    double next_substitution(NucleotideRateMatrix const &nuc_matrix, double beta, double gBGC,
        double time_start, double time_end, bool burn_in) {
        // Number of possible substitutions is 9 times the number of sites (3 substitutions for each
        // 3 possible positions).
        DFE dfe{};
        u_long nbr_substitutions{9 * nbr_sites};

        // Vector of substitution rates.
        std::vector<double> substitution_rates(nbr_substitutions, 0);

        // Sum of substitution rates.
        double total_substitution_rates{0.};

        double non_syn_sub_flow_WS{0.0}, non_syn_mut_flow_WS{0.0}, non_syn_sub_flow_SW{0.0},
            non_syn_mut_flow_SW{0.0}, non_syn_sub_flow{0.0}, non_syn_mut_flow{0.0},
            syn_mut_flow{0.0}, dndn0{0.0};
        u_long ns{0};

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
                if (codonLexico.codon_to_aa[codon_to] == 20) { continue; }
                double rate_mutation = nuc_matrix(n_from, n_to);

                double s{0.0};
                bool w2s = weak_nuc.count(n_from) != 0 and strong_nuc.count(n_to) != 0;
                bool s2w = strong_nuc.count(n_from) != 0 and weak_nuc.count(n_to) != 0;
                if (w2s) {
                    s += gBGC;
                } else if (s2w) {
                    s -= gBGC;
                }
                bool non_syn_mut =
                    codonLexico.codon_to_aa[codon_from] != codonLexico.codon_to_aa[codon_to];
                if (non_syn_mut) {
                    s += fitness_state->selection_coefficient(
                        codon_seq, site, codon_to, burn_in, beta);
                }
                double pfix = fitness_state->fixation_rate(beta, s);
                if (!std::isfinite(pfix)) { continue; }
                double rate_substitution = rate_mutation * pfix;

                if (non_syn_mut) {
                    non_syn_mut_flow += rate_mutation;
                    non_syn_sub_flow += rate_substitution;
                    dndn0 += pfix;
                    if (w2s) {
                        non_syn_mut_flow_WS += rate_mutation;
                        non_syn_sub_flow_WS += rate_substitution;
                    } else if (s2w) {
                        non_syn_mut_flow_SW += rate_mutation;
                        non_syn_sub_flow_SW += rate_substitution;
                    }
                    ns++;

                    if (!burn_in) {
                        if (pfix > 1) {
                            dfe.advantageous++;
                            dfe.sub_advantageous += pfix;
                        } else {
                            dfe.deleterious++;
                            dfe.sub_deleterious += pfix;
                        }
                        dfe.sel_coeff_subs += s * pfix;
                        distribution_map["DFE"].add(4 * beta * s);
                    }
                } else {
                    syn_mut_flow += rate_mutation;
                }

                substitution_rates[9 * site + neighbor] = rate_substitution;
                // Increment the sum of substitution rates
                total_substitution_rates += rate_substitution;
            }
        }

        dndn0 /= static_cast<double>(ns);

        // Increment the time by drawing from an exponential distribution (mean equal to the inverse
        // sum of substitution rates).
        double time_draw =
            std::exponential_distribution<double>(total_substitution_rates)(generator);

        // Substitute the sequence if the time is positive, else there is no substitution but the
        // time left is set to 0.
        assert(total_substitution_rates != 0.0);
        if (time_start + time_draw <= time_end) {
            std::discrete_distribution<u_long> substitution_distr(
                substitution_rates.begin(), substitution_rates.end());

            u_long index = substitution_distr(generator);
            u_long site = index / 9;
            char codom_from = codon_seq[site];

            // Array of neighbors of the original codon (codons differing by only 1 mutation).
            std::array<std::tuple<char, char, char>, 9> neighbors =
                codonLexico.codon_to_neighbors[codom_from];

            // Codon after mutation, Nucleotide original and Nucleotide after mutation.
            char codon_to{0}, n_from{0}, n_to{0};

            std::tie(codon_to, n_from, n_to) = neighbors[index % 9];

            assert(non_syn_mut_flow < 10e10);
            assert(syn_mut_flow < 10e10);

            if (!burn_in) {
                substitutions.emplace(time_start + time_draw, time_draw, non_syn_sub_flow,
                    non_syn_mut_flow, non_syn_sub_flow_WS, non_syn_mut_flow_WS, non_syn_sub_flow_SW,
                    non_syn_mut_flow_SW, syn_mut_flow, dndn0, codom_from, codon_to, n_from, n_to,
                    site + position, dfe);
            }

            if (codonLexico.codon_to_aa[codon_seq[site]] != codonLexico.codon_to_aa[codon_to]) {
                fitness_state->update(codon_seq, site, codon_to, burn_in, beta);
            }
            codon_seq[site] = codon_to;
            time_start += time_draw;
        } else {
            if (!burn_in) {
                substitutions.emplace(time_end, time_end - time_start, non_syn_sub_flow,
                    non_syn_mut_flow, non_syn_sub_flow_WS, non_syn_mut_flow_WS, non_syn_sub_flow_SW,
                    non_syn_mut_flow_SW, syn_mut_flow, dndn0);
            }
            time_start = time_end;
        }
        return time_start;
    }

    // Method computing all substitution event occurring during a given time-frame.
    // t: time during which substitution events occur (typically branch length).
    // This method is o(n²) where n is the number of sites, but can take into account epistatic
    // effects
    void run_substitutions(NucleotideRateMatrix const &nuc_matrix, double beta, double gBGC,
        double t_start, double t_end) {
        while (t_start < t_end) {
            t_start = next_substitution(nuc_matrix, beta, gBGC, t_start, t_end, false);
        }
    }
};

static double time_grid_step;
Trace tracer_traits;
Trace tracer_fossils;
Trace tracer_substitutions;

class Sequence {
  public:
    // TimeElapsed
    double time_from_root{0};

    // Blocks
    std::vector<Exon> exons;

    LogMultivariate log_multivariate;
    double beta{0};
    double generation_time{0};
    double gBGC{0};
    NucleotideRateMatrix nuc_matrix;

    EMatrix const &transform_matrix;
    BiasMultivariate &bias_vector;
    bool branch_wise;
    std::vector<Substitution> interspersed_substitutions;
    PieceWiseMultivariate piecewise_multivariate{};

    Sequence(FitnessModel &seq_fitness, LogMultivariate &log_multi,
        NucleotideRateMatrix &nucleotide_matrix, EMatrix const &transform_matrix,
        BiasMultivariate &bias_vector, bool branch_wise)
        : exons{},
          log_multivariate{log_multi},
          beta{log_multivariate.beta()},
          generation_time{log_multivariate.generation_time()},
          gBGC{log_multivariate.gBGC()},
          nuc_matrix{nucleotide_matrix},
          transform_matrix{transform_matrix},
          bias_vector{bias_vector},
          branch_wise{branch_wise} {
        u_long pos = 0;
        for (std::unique_ptr<FitnessState> &exon_seq_fitness : seq_fitness.fitness_states) {
            exons.emplace_back(exon_seq_fitness, pos, beta);
            pos += exons.back().nbr_sites;
        }
    }

    bool operator==(Sequence const &other) const {
        if (nbr_exons() != other.nbr_exons()) { return false; }
        for (std::size_t i = 0; i < nbr_exons(); ++i) {
            if (exons.at(i) != other.exons.at(i)) { return false; }
        }
        return true;
    };

    u_long nbr_exons() const { return exons.size(); }

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
                char codon_to =
                    static_cast<char>(std::distance(codonLexico.codon_to_aa.begin(), it));
                exon.codon_seq[site] = codon_to;
            }
            exon.fitness_state->update(exon.codon_seq, beta);
        }
    }

    void set_from_dna_seq(std::string const &dna_string) {
        assert(dna_string.size() == nbr_nucleotides());
        for (auto &exon : exons) {
            for (u_long site{0}; site < exon.nbr_sites; site++) {
                char n1 = Codon::nuc_to_index[dna_string.at(3 * (exon.position + site))];
                char n2 = Codon::nuc_to_index[dna_string.at(3 * (exon.position + site) + 1)];
                char n3 = Codon::nuc_to_index[dna_string.at(3 * (exon.position + site) + 2)];
                char codon_to = Codon::triplet_to_codon(n1, n2, n3);
                exon.codon_seq[site] = codon_to;
            }
            exon.fitness_state->update(exon.codon_seq, beta);
        }
    }

    bool check_consistency() const {
        for (std::size_t i = 1; i < interspersed_substitutions.size(); i++) {
            if (abs(interspersed_substitutions.at(i).time_between_event -
                    (interspersed_substitutions.at(i).time_event -
                        interspersed_substitutions.at(i - 1).time_event)) > 1e-6) {
                std::cerr
                    << "The time between interspersed substitutions does not match the timing of "
                       "events."
                    << std::endl;
                return false;
            }
        }
        if (!is_sorted(interspersed_substitutions.begin(), interspersed_substitutions.end(),
                [](Substitution const &a, Substitution const &b) -> bool {
                    return a.time_event <= b.time_event;
                })) {
            std::cerr << "The interspersed substitutions are not sorted." << std::endl;
            return false;
        }
        return true;
    };

    void intersperse_exon_substitutions() {
        // Find time of the next substitutions
        // At the same time sum the mut flow and sub flow
        while (true) {
            auto exon_next_sub =
                min_element(exons.begin(), exons.end(), [](Exon const &a, Exon const &b) -> bool {
                    return a.substitutions.front().time_event <= b.substitutions.front().time_event;
                });
            auto sub = exon_next_sub->substitutions.front();
            if (sub.is_dummy()) {
                for (auto &exon : exons) {
                    assert(exon.substitutions.front().is_dummy());
                    assert(exon.substitutions.front().time_event == sub.time_event);
                    exon.substitutions.pop();
                    assert(exon.substitutions.empty());
                }
                break;
            } else {
                if (!interspersed_substitutions.empty()) {
                    sub.time_between_event =
                        sub.time_event - interspersed_substitutions.back().time_event;
                }
                sub.non_syn_sub_flow = 0, sub.non_syn_mut_flow = 0;
                sub.non_syn_sub_flow_WS = 0, sub.non_syn_mut_flow_WS = 0;
                sub.non_syn_sub_flow_SW = 0, sub.non_syn_mut_flow_SW = 0;
                sub.syn_mut_flow = 0;
                sub.dfe = DFE();
                for (auto const &exon : exons) {
                    Substitution const &p = exon.substitutions.front();
                    sub.non_syn_sub_flow += p.non_syn_sub_flow;
                    sub.non_syn_mut_flow += p.non_syn_mut_flow;
                    sub.non_syn_sub_flow_WS += p.non_syn_sub_flow_WS;
                    sub.non_syn_mut_flow_WS += p.non_syn_mut_flow_WS;
                    sub.non_syn_sub_flow_SW += p.non_syn_sub_flow_SW;
                    sub.non_syn_mut_flow_SW += p.non_syn_mut_flow_SW;
                    sub.syn_mut_flow += p.syn_mut_flow;
                    sub.dndn0 += p.dndn0 * static_cast<double>(exon.nbr_sites);
                    sub.dfe += p.dfe;
                }
                sub.dndn0 /= static_cast<double>(nbr_sites());
                interspersed_substitutions.push_back(sub);
                exon_next_sub->substitutions.pop();
            }
        }
    };

    // Set the the DNA sequence to the mutation-selection equilibrium.
    void at_equilibrium(u_long nbr_rounds = 1, double init_pop_size = -1, double init_gBGC = -1) {
        // For all site of the sequence.
        if (init_pop_size < 0) { init_pop_size = beta; }
        if (init_gBGC < 0) { init_gBGC = gBGC; }
        for (u_long i = 0; i < nbr_rounds; i++) {
            for (auto &exon : exons) {
                for (u_long site{0}; site < exon.nbr_sites; site++) {
                    std::array<double, 64> codon_freqs = exon.fitness_state->codon_frequencies(
                        exon.fitness_state->aa_selection_coefficients(
                            exon.codon_seq, site, init_pop_size),
                        nuc_matrix, init_pop_size, gBGC);
                    std::discrete_distribution<char> freq_codon_distr(
                        codon_freqs.begin(), codon_freqs.end());
                    char chosen_codon = freq_codon_distr(generator);
                    if (codonLexico.codon_to_aa[chosen_codon] !=
                        codonLexico.codon_to_aa[exon.codon_seq[site]]) {
                        exon.fitness_state->update(
                            exon.codon_seq, site, chosen_codon, true, init_pop_size);
                    }
                    exon.codon_seq[site] = chosen_codon;
                }
            }
        }
    }

    void burn_in(int nbr_sub) {
        for (int i = 0; i < nbr_sub; i++) {
            for (auto &exon : exons) {
                exon.next_substitution(
                    nuc_matrix, beta, gBGC, 0.0, std::numeric_limits<double>::infinity(), true);
            }
        }
        clear();
        std::cout << get_aa_str() << std::endl;
    }

    EVector delta_log_multivariate(double distance) const {
        EVector normal_vector = EVector::Zero(dimensions);
        for (int dim = 0; dim < dimensions; dim++) {
            normal_vector(dim) = normal_distrib(generator_brownian);
        }
        return sqrt(distance) * (transform_matrix * normal_vector);
    }

    void update_brownian(EVector const &delta) {
        log_multivariate += delta;

        beta = log_multivariate.beta();
        generation_time = log_multivariate.generation_time();
        gBGC = log_multivariate.gBGC();
        nuc_matrix.set_mutation_rate(
            log_multivariate.mutation_rate_per_generation() / generation_time);
    }

    void run_forward(double t_max, Tree const &tree) {
        if (bias_vector.step_wise_pop_size) {
            piecewise_multivariate.AddMultivariate(t_max, log_multivariate);
            for (auto &exon : exons) { exon.run_substitutions(nuc_matrix, beta, gBGC, 0.0, t_max); }
            intersperse_exon_substitutions();
            double d = t_max / tree.max_distance_to_root();
            update_brownian(delta_log_multivariate(d) + d * bias_vector);
            bias_vector(dim_population_size) = -bias_vector(dim_population_size);
            time_from_root += t_max;
            assert(check_consistency());
        } else if (branch_wise) {
            double d = t_max / tree.max_distance_to_root();
            EVector delta = delta_log_multivariate(d) + d * bias_vector;
            update_brownian(delta / 2);
            piecewise_multivariate.AddMultivariate(t_max, log_multivariate);
            for (auto &exon : exons) { exon.run_substitutions(nuc_matrix, beta, gBGC, 0.0, t_max); }
            intersperse_exon_substitutions();
            update_brownian(delta / 2);
            time_from_root += t_max;
            assert(check_consistency());
        } else {
            double time_current = 0.0;
            double step_in_year = time_grid_step;
            while (time_current < t_max) {
                if (time_current + step_in_year > t_max) { step_in_year = t_max - time_current; }
                double d = step_in_year / tree.max_distance_to_root();
                update_brownian(delta_log_multivariate(d) + d * bias_vector);
                piecewise_multivariate.AddMultivariate(step_in_year, log_multivariate);
                for (auto &exon : exons) {
                    exon.run_substitutions(
                        nuc_matrix, beta, gBGC, time_current, time_current + step_in_year);
                }
                intersperse_exon_substitutions();
                time_current += step_in_year;
                time_from_root += step_in_year;
                assert(check_consistency());
            }
        }
    }

    double predicted_dn_dn0(NucleotideRateMatrix const &rates, double const &relative_pop,
        double const &relative_gBGC) const {
        double dn{0.}, dn0{0.};
        for (auto const &exon : exons) {
            double exon_dn{0}, exon_d0{0};
            std::tie(exon_dn, exon_d0) = exon.fitness_state->predicted_dn_dn0(
                exon.codon_seq, rates, relative_pop, relative_gBGC);
            dn += exon_dn;
            dn0 += exon_d0;
        }
        return dn / dn0;
    };

    double sequence_wise_predicted_dn_dn0(Sequence const &parent, NucleotideRateMatrix const &rates,
        double const &relative_pop, double const &relative_gBGC) const {
        double dn{0.}, dn0{0.};

        for (std::size_t i = 0; i < nbr_exons(); i++) {
            double exon_dn{0}, exon_dn0{0};
            std::tie(exon_dn, exon_dn0) = exons[i].fitness_state->flow_dn_dn0(
                parent.exons[i].codon_seq, rates, relative_pop, relative_gBGC);
            dn += exon_dn;
            dn0 += exon_dn0;
        }
        return dn / dn0;
    };

    double count_based_dn_dn0() const {
        double dn{0}, dn0{0};
        for (auto const &substitution : interspersed_substitutions) {
            if (substitution.is_non_synonymous()) { dn++; }
            dn0 += substitution.non_syn_mut_flow * substitution.time_between_event;
        }
        if (dn0 == .0) {
            return .0;
        } else {
            return dn / dn0;
        }
    }

    double flow_based_dn_dn0_WS() const {
        double dn{0}, dn0{0};
        for (auto const &substitution : interspersed_substitutions) {
            dn0 += substitution.non_syn_mut_flow_WS * substitution.time_between_event;
            dn += substitution.non_syn_sub_flow_WS * substitution.time_between_event;
        }
        if (dn0 == .0) {
            return .0;
        } else {
            return dn / dn0;
        }
    }

    double flow_based_dn_dn0_SW() const {
        double dn{0}, dn0{0};
        for (auto const &substitution : interspersed_substitutions) {
            dn0 += substitution.non_syn_mut_flow_SW * substitution.time_between_event;
            dn += substitution.non_syn_sub_flow_SW * substitution.time_between_event;
        }
        if (dn0 == .0) {
            return .0;
        } else {
            return dn / dn0;
        }
    }

    double flow_based_dn_dn0() const {
        double dn{0}, dn0{0};
        for (auto const &substitution : interspersed_substitutions) {
            dn0 += substitution.non_syn_mut_flow * substitution.time_between_event;
            dn += substitution.non_syn_sub_flow * substitution.time_between_event;
        }
        if (dn0 == .0) {
            return .0;
        } else {
            return dn / dn0;
        }
    }

    // Simulated omega from the substitutions
    double event_based_dn_ds() const {
        double dn{0}, ds{0}, dn0{0}, ds0{0};
        for (auto const &substitution : interspersed_substitutions) {
            dn0 += substitution.non_syn_mut_flow;
            ds0 += substitution.syn_mut_flow;
            if (substitution.is_synonymous()) {
                ds++;
            } else if (substitution.is_non_synonymous()) {
                dn++;
            }
        }
        if (ds == .0) {
            std::cerr
                << "There is no synonymous substitutions, event based dN/dS can't be computed !"
                << std::endl;
            return .0;
        } else {
            return (dn * ds0) / (dn0 * ds);
        }
    }

    // Simulated omega from the substitutions
    double count_based_dn_ds() const {
        double dn{0}, ds{0}, dn0{0}, ds0{0};
        for (auto const &substitution : interspersed_substitutions) {
            dn0 += substitution.non_syn_mut_flow * substitution.time_between_event;
            ds0 += substitution.syn_mut_flow * substitution.time_between_event;
            if (substitution.is_synonymous()) {
                ds++;
            } else if (substitution.is_non_synonymous()) {
                dn++;
            }
        }
        if (ds == .0) {
            std::cerr
                << "There is no synonymous substitutions, count based dN/dS can't be computed !"
                << std::endl;
            return .0;
        } else {
            return (dn * ds0) / (dn0 * ds);
        }
    }

    void node_trace(std::string const &output_filename, Tree::NodeIndex node, Tree &tree,
        std::unique_ptr<Sequence> const &parent) const {
        std::string node_name = tree.node_name(node);

        tree.set_tag(node, "population_size", d_to_string(beta));
        tree.set_tag(node, "generation_time", d_to_string(generation_time));
        tree.set_tag(node, "gBGC", d_to_string(gBGC));
        tree.set_tag(node, "mutation_rate", d_to_string(nuc_matrix.mutation_rate));
        tree.set_tag(node, "mutation_rate_per_generation",
            d_to_string(log_multivariate.mutation_rate_per_generation()));

        if (!tree.is_leaf(node)) {
            tracer_fossils.add("NodeName", node_name);
            double age = tree.max_distance_to_root() - time_from_root;
            tracer_fossils.add("Age", age);
            tracer_fossils.add("LowerBound", age);
            tracer_fossils.add("UpperBound", age);
        }

        if (tree.is_root(node)) { return; }
        assert(parent != nullptr);
        double geom_pop_size = piecewise_multivariate.GeometricPopSize();
        piecewise_multivariate.add_to_tree(tree, node, geom_pop_size);
        double geom_gBGC = piecewise_multivariate.GeometricgBGC();

        tree.set_tag(node, "Branch_dNdN0_predicted",
            d_to_string(predicted_dn_dn0(nuc_matrix, geom_pop_size, geom_gBGC)));
        tree.set_tag(node, "Branch_dNdN0_sequence_wise_predicted",
            d_to_string(
                sequence_wise_predicted_dn_dn0(*parent, nuc_matrix, geom_pop_size, geom_gBGC)));

        double flow_dn_dn0 = flow_based_dn_dn0();
        double flow_dn_dn0_WS = flow_based_dn_dn0_WS();
        double flow_dn_dn0_SW = flow_based_dn_dn0_SW();
        double count_dn_dn0 = count_based_dn_dn0();
        double event_dn_ds = event_based_dn_ds();
        double count_dn_ds = count_based_dn_ds();
        pfix_map["dnd0_event_tot"] += flow_dn_dn0 * tree.node_length(node);
        pfix_map["dnd0_WS_event_tot"] += flow_dn_dn0_WS * tree.node_length(node);
        pfix_map["dnd0_SW_event_tot"] += flow_dn_dn0_SW * tree.node_length(node);
        pfix_map["dnd0_count_tot"] += count_dn_dn0 * tree.node_length(node);
        pfix_map["dnds_event_tot"] += event_dn_ds * tree.node_length(node);
        pfix_map["dnds_count_tot"] += count_dn_ds * tree.node_length(node);
        tree.set_tag(node, "Branch_dNdN0_flow_based", d_to_string(flow_dn_dn0));
        tree.set_tag(node, "Branch_dNdN0_WS_flow_based", d_to_string(flow_dn_dn0_WS));
        tree.set_tag(node, "Branch_dNdN0_SW_flow_based", d_to_string(flow_dn_dn0_SW));
        tree.set_tag(node, "Branch_dNdN0_count_based", d_to_string(count_dn_dn0));
        tree.set_tag(node, "Branch_dNdS_event_based", d_to_string(event_dn_ds));
        tree.set_tag(node, "Branch_dNdS_count_based", d_to_string(count_dn_ds));

        for (auto const &sub : interspersed_substitutions) {
            if (!sub.is_dummy()) {
                tracer_substitutions.add("NodeName", node_name);
                tracer_substitutions.add("EndTime", sub.time_event);
                tracer_substitutions.add(
                    "AbsoluteStartTime", (time_from_root - tree.node_length(node)) +
                                             (sub.time_event - sub.time_between_event));
                tracer_substitutions.add("NucFrom", Codon::nucleotide(sub.n_from));
                tracer_substitutions.add("NucTo", Codon::nucleotide(sub.n_to));
                tracer_substitutions.add("CodonFrom", codonLexico.codon_string(sub.codon_from));
                tracer_substitutions.add("CodonTo", codonLexico.codon_string(sub.codon_to));
                tracer_substitutions.add("AAFrom", codonLexico.codon_aa_string(sub.codon_from));
                tracer_substitutions.add("AATo", codonLexico.codon_aa_string(sub.codon_to));
                tracer_substitutions.add("Site", sub.site);
                tracer_substitutions.add("<dN>>", sub.non_syn_sub_flow);
                tracer_substitutions.add("<dN0>", sub.non_syn_mut_flow);
                tracer_substitutions.add("<dN>/<dN0>", sub.non_syn_sub_flow / sub.non_syn_mut_flow);
                tracer_substitutions.add("<dN>>_WS", sub.non_syn_sub_flow_WS);
                tracer_substitutions.add("<dN0>_WS", sub.non_syn_mut_flow_WS);
                tracer_substitutions.add(
                    "<dN>/<dN0>_WS", sub.non_syn_sub_flow_WS / sub.non_syn_mut_flow_WS);
                tracer_substitutions.add("<dS>", sub.syn_mut_flow);
                tracer_substitutions.add("<dN/dN0>", sub.dndn0);
            }
        }
        write_sequence(output_filename, node_name, this->get_dna_str(), true);
        if (!tree.is_leaf(node)) { return; }
        // If the node is a leaf, output the DNA sequence and name.
        write_sequence(output_filename, node_name, this->get_dna_str());

        tracer_traits.add("TaxonName", node_name);
        tracer_traits.add("LogPopulationSize", log_multivariate.log_population_size());
        tracer_traits.add(
            "LogMutationRatePerGeneration", log_multivariate.log_mutation_rate_per_generation());
        tracer_traits.add("LogGenerationTime", log_multivariate.log_generation_time());
        tracer_traits.add("gBGC", log_multivariate.gBGC());
    }

    // Method returning the DNA std::string corresponding to the codon sequence.
    std::string get_dna_str() const {
        std::string dna_str{};
        dna_str.reserve(nbr_nucleotides());

        // For each site of the sequence.
        for (auto const &exon : exons) {
            for (u_long site{0}; site < exon.nbr_sites; site++) {
                // Assert there is no stop in the sequence.
                assert(codonLexico.codon_to_aa[exon.codon_seq[site]] != 20);

                // Translate the site to a triplet of DNA nucleotides
                std::array<char, 3> triplet = codonLexico.codon_to_triplet[exon.codon_seq[site]];
                for (char position{0}; position < 3; position++) {
                    dna_str += Codon::nucleotide(triplet[position]);
                }
            }
        }
        return dna_str;  // return the DNA sequence as a string.
    }

    // Method returning the AA std::string corresponding to the codon sequence.
    std::string get_aa_str() const {
        std::string aa_str{};
        aa_str.reserve(nbr_sites());

        // For each site of the sequence.
        for (auto const &exon : exons) {
            for (u_long site{0}; site < exon.nbr_sites; site++) {
                // Assert there is no stop in the sequence.
                assert(codonLexico.codon_to_aa[exon.codon_seq[site]] != 20);
                // Translate the site to amino-acids
                aa_str += codonLexico.codon_aa_string(exon.codon_seq[site]);
            }
        }
        return aa_str;  // return the DNA sequence as a string.
    }

    void clear() {
        interspersed_substitutions.clear();
        piecewise_multivariate.clear();
        for (auto &exon : exons) {
            while (!exon.substitutions.empty()) { exon.substitutions.pop(); }
        }
    }
};

class Process {
  private:
    static double years_computed;
    Tree &tree;
    std::vector<std::unique_ptr<Sequence>> sequences;  // Vector of sequence DNA.

  public:
    // Constructor
    Process(Tree &intree, std::unique_ptr<Sequence> &root_seq) : tree{intree}, sequences() {
        sequences.resize(tree.nb_nodes());
        sequences[tree.root()] = std::move(root_seq);
    }

    void run(std::string const &output_filename, bool branch_length_in_dS_unit = false) {
        sequences[tree.root()]->node_trace(output_filename, tree.root(), tree, nullptr);
        run_recursive(tree.root(), output_filename);
        std::ofstream nhx;
        nhx.open(output_filename + ".nhx");
        nhx << tree.as_string(true, branch_length_in_dS_unit) << std::endl;
        nhx.close();
        nhx.open(output_filename + ".tree");
        nhx << tree.as_string(false, branch_length_in_dS_unit) << std::endl;
        nhx.close();
    }

    // Recursively iterate through the subtree.
    void run_recursive(Tree::NodeIndex node, std::string const &output_filename) {
        // Substitutions of the DNA sequence is generated.
        sequences.at(node)->clear();

        if (!tree.is_root(node)) {
            sequences.at(node)->run_forward(tree.node_length(node), tree);

            years_computed += tree.node_length(node);
            std::cout << years_computed << " years computed in total ("
                      << static_cast<int>(100 * years_computed / tree.total_length())
                      << "%) at node " << tree.node_name(node) << " ("
                      << static_cast<int>(100 * tree.node_length(node) / tree.total_length())
                      << "%)." << std::endl;

            sequences.at(node)->node_trace(
                output_filename, node, tree, sequences.at(tree.parent(node)));
        }

        // If the node is internal, iterate through the direct children.
        for (auto &child : tree.children(node)) {
            sequences.at(child) = std::make_unique<Sequence>(*sequences.at(node));
        }
        for (auto &child : tree.children(node)) {
            assert(*sequences.at(node) == *sequences.at(child));
            run_recursive(child, output_filename);
        }
    }

    // Recursively iterate through the subtree and count the number of substitutions.
    void summary(std::string const &output_path, double expected_subs) {
        u_long nbr_non_syn{0}, nbr_syn{0};

        DFE dfe{};
        for (std::unique_ptr<Sequence> const &seq : sequences) {
            for (auto const &s : seq->interspersed_substitutions) { dfe += s.dfe; }
            nbr_syn += count_if(seq->interspersed_substitutions.begin(),
                seq->interspersed_substitutions.end(),
                [](Substitution const &s) { return s.is_synonymous(); });
            nbr_non_syn += count_if(seq->interspersed_substitutions.begin(),
                seq->interspersed_substitutions.end(),
                [](Substitution const &s) { return s.is_non_synonymous(); });
        }

        // .txt output
        Trace trace;
        trace.add("#substitutions", nbr_syn + nbr_non_syn);
        double n_sites = static_cast<double>(sequences.begin()->get()->nbr_sites());
        trace.add("#substitutions_per_site", static_cast<double>(nbr_syn + nbr_non_syn) / n_sites);
        trace.add("#synonymous_substitutions", nbr_syn);
        trace.add("#non_synonymous_substitutions", nbr_non_syn);
        for (auto &v : pfix_map) {
            v.second /= tree.total_length();
            std::cout << v.first << " is :" << v.second << std::endl;
            trace.add(v.first, v.second);
        }

        dfe.add(trace);
        for (auto const &ss : FitnessState::summary_stats) {
            trace.add(ss.first + "-mean", ss.second.mean());
            trace.add(ss.first + "-std", ss.second.std());
            trace.add(ss.first + "-abs-mean", ss.second.abs_mean());
            trace.add(ss.first + "-abs-std", ss.second.abs_std());
        }
        trace.write_tsv(output_path);

        tracer_traits.write_tsv(output_path + ".traits");
        tracer_fossils.write_tsv(output_path + ".fossils");
        tracer_substitutions.write_tsv(output_path + ".substitutions");
        distribution_map.write_tsv(output_path);

        std::cout << "Simulation computed." << std::endl;
        std::cout << expected_subs << " expected substitutions." << std::endl;
        std::cout << nbr_syn + nbr_non_syn << " simulated substitutions." << std::endl;
        std::cout << nbr_non_syn << " simulated non-synonymous substitutions." << std::endl;
        std::cout << nbr_syn << " simulated synonymous substitutions." << std::endl;
        std::cout << "Statistics summarized in: " << output_path + ".tsv" << std::endl;
        std::cout << "Fasta file in: " << output_path + ".fasta" << std::endl;
        std::cout << "Alignment (.ali) file in: " << output_path + ".ali" << std::endl;
    }
};

// Initialize static variables
double Process::years_computed = 0.0;