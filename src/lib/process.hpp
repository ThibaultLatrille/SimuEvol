#pragma once

#include <queue>
#include "fitness.hpp"
#include "matrices.hpp"
#include "random.hpp"

using namespace std;

double dnds_count_tot{0}, dnds_event_tot{0}, dnd0_count_tot{0}, dnd0_event_tot{0};
double s_tot{0}, S_tot{0}, s_abs_tot{0}, S_abs_tot{0}, pfix_tot{0};

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
    double syn_mut_flow;

    Substitution(double time_event, double time_between_event, double non_syn_sub_flow,
                 double non_syn_mut_flow, double syn_mut_flow, char codon_from = -1, char codon_to = -1,
                 char n_from = -1, char n_to = -1, u_long site = 0)
            : site{site},
              codon_from{codon_from},
              codon_to{codon_to},
              n_from{n_from},
              n_to{n_to},
              time_event{time_event},
              time_between_event{time_between_event},
              non_syn_sub_flow{non_syn_sub_flow},
              non_syn_mut_flow{non_syn_mut_flow},
              syn_mut_flow{syn_mut_flow} {}

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
    vector<char> codon_seq;

    // The fitness landscape (cannot be modified), which map genotype->phenotype->fitness.
    FitnessLandscape const &fitness_landscape;

    queue<Substitution> substitutions{};

    // Constructor of Exon.
    // size: the size of the DNA sequence.
    explicit Exon(FitnessLandscape const &exon_fitness, u_long const &position)
            : nbr_sites{exon_fitness.nbr_sites()},
              position{position},
              codon_seq(nbr_sites, 0),
              fitness_landscape{exon_fitness} {
        assert(substitutions.empty());

    }

    // Method computing the next substitution event to occur, and the time for it to happen.
    // This method takes a time as input. If there is enough time given, the substitution event is
    // computed. This method also returns the time (given as input) decremented by the time needed
    // for the substitution event to occur. If there is not enough time given, no substitution event
    // is computed and the method returns 0. time_left: The time available for a substitution event
    // to occur.
    double next_substitution(
            NucleotideRateMatrix const &nuc_matrix, double beta, double time_start, double time_end, bool only_non_syn) {
        // Number of possible substitutions is 9 times the number of sites (3 substitutions for each
        // 3 possible positions).
        u_long nbr_substitutions{9 * nbr_sites};

        // Vector of substitution rates.
        vector<double> substitution_rates(nbr_substitutions, 0);

        // Sum of substitution rates.
        double total_substitution_rates{0.};

        double non_syn_sub_flow{0.0}, non_syn_mut_flow{0.0}, syn_mut_flow{0.0};

        // For all site of the sequence.
        for (u_long site{0}; site < nbr_sites; site++) {
            // Codon original before substitution.

            char codon_from = codon_seq[site];

            // Array of neighbors of the original codon (codons differing by only 1 mutation).
            array<tuple<char, char, char>, 9> neighbors =
                    codonLexico.codon_to_neighbors[codon_from];

            // For all possible neighbors.
            for (char neighbor{0}; neighbor < 9; neighbor++) {
                // Codon after mutation, Nucleotide original and Nucleotide after mutation.
                char codon_to{0}, n_from{0}, n_to{0};
                tie(codon_to, n_from, n_to) = neighbors[neighbor];

                // Assign the substitution rate given by the method substitution rate.
                // Rate of substitution initialized to 0 (deleterious mutation)
                double rate_substitution{0.};

                // If the mutated amino-acid is a stop codon, the rate of fixation is 0.
                // Else, if the mutated and original amino-acids are non-synonymous, we compute the
                // rate of fixation. Note that, if the mutated and original amino-acids are
                // synonymous, the rate of fixation is 1.
                if (codonLexico.codon_to_aa[codon_to] == 20) { continue; }

                rate_substitution = nuc_matrix(n_from, n_to);
                if (codonLexico.codon_to_aa[codon_from] != codonLexico.codon_to_aa[codon_to]) {
                    non_syn_mut_flow += rate_substitution;
                    rate_substitution *= Pfix(beta, fitness_landscape.selection_coefficient(codon_seq, site, codon_to));
                    non_syn_sub_flow += rate_substitution;
                } else {
                    if (only_non_syn) { continue; }
                    syn_mut_flow += rate_substitution;
                }

                substitution_rates[9 * site + neighbor] = rate_substitution;
                // Increment the sum of substitution rates
                total_substitution_rates += rate_substitution;
            }
        }

        // Increment the time by drawing from an exponential distribution (mean equal to the inverse
        // sum of substitution rates).
        double time_draw = exponential_distribution<double>(total_substitution_rates)(generator);

        // Substitute the sequence if the time is positive, else there is no substitution but the
        // time left is set to 0.
        assert(total_substitution_rates != 0.0);
        if (time_start + time_draw <= time_end) {
            discrete_distribution<u_long> substitution_distr(
                    substitution_rates.begin(), substitution_rates.end());

            u_long index = substitution_distr(generator);
            u_long site = index / 9;
            char codom_from = codon_seq[site];

            // Array of neighbors of the original codon (codons differing by only 1 mutation).
            array<tuple<char, char, char>, 9> neighbors =
                    codonLexico.codon_to_neighbors[codom_from];

            // Codon after mutation, Nucleotide original and Nucleotide after mutation.
            char codon_to{0}, n_from{0}, n_to{0};

            tie(codon_to, n_from, n_to) = neighbors[index % 9];

            assert(non_syn_mut_flow < 10e10);
            assert(non_syn_sub_flow < 10e10);
            assert(syn_mut_flow < 10e10);

            substitutions.emplace(time_start + time_draw, time_draw, non_syn_sub_flow,
                                  non_syn_mut_flow, syn_mut_flow, codom_from, codon_to, n_from, n_to, site);

            if (codonLexico.codon_to_aa[codon_seq[site]] != codonLexico.codon_to_aa[codon_to]) {
                double s = fitness_landscape.selection_coefficient(codon_seq, site, codon_to);
                if (!only_non_syn) {
                    s_tot += s;
                    s_abs_tot += abs(s);
                    S_tot += 4 * s * beta;
                    S_abs_tot += abs(4 * s * beta);
                    pfix_tot += Pfix(beta, s);
                }
            }
            codon_seq[site] = codon_to;
            time_start += time_draw;
        } else {
            substitutions.emplace(
                    time_end, time_end - time_start, non_syn_sub_flow, non_syn_mut_flow, syn_mut_flow);
            time_start = time_end;
        }
        return time_start;
    }

    // Method computing all substitution event occurring during a given time-frame.
    // t: time during which substitution events occur (typically branch length).
    // This method is o(n²) where n is the number of sites, but can take into account epistatic
    // effects
    void run_substitutions(
            NucleotideRateMatrix const &nuc_matrix, double beta, double t_start, double t_end) {
        while (t_start < t_end) { t_start = next_substitution(nuc_matrix, beta, t_start, t_end, false); }
    }
};

static double time_grid_step;
Trace tracer_traits;
Trace tracer_fossils;
Trace tracer_substitutions;
Trace tracer_sequences;

class Sequence {
public:
    // TimeElapsed
    double time_from_root{0};

    // Blocks
    vector<Exon> exons;

    LogMultivariate log_multivariate;
    double beta{0};
    double generation_time{0};
    NucleotideRateMatrix nuc_matrix;

    EMatrix const &transform_matrix;
    bool branch_wise;
    vector<Substitution> interspersed_substitutions;
    PieceWiseMultivariate piecewise_multivariate{};

    Sequence(SequenceFitnessLandscape &seq_fitness, LogMultivariate &log_multi, NucleotideRateMatrix nucleotide_matrix,
             EMatrix const &transform_matrix, bool branch_wise)
            : exons{},
              log_multivariate{log_multi},
              beta{log_multivariate.beta()},
              generation_time{log_multivariate.generation_time()},
              nuc_matrix{move(nucleotide_matrix)},
              transform_matrix{transform_matrix},
              branch_wise{branch_wise} {
        u_long pos = 0;
        for (auto &exon_seq_fitness: seq_fitness) {
            exons.emplace_back(Exon(*exon_seq_fitness, pos));
            pos += exon_seq_fitness->nbr_sites();
        }
    }

    u_long nbr_sites() const {
        u_long sites = 0;
        for (auto const &exon : exons) { sites += exon.nbr_sites; }
        return sites;
    }

    u_long nbr_nucleotides() const { return 3 * nbr_sites(); }

    bool check_consistency() const {
        for (size_t i = 1; i < interspersed_substitutions.size(); i++) {
            if (abs(interspersed_substitutions.at(i).time_between_event -
                    (interspersed_substitutions.at(i).time_event -
                     interspersed_substitutions.at(i - 1).time_event)) > 1e-6) {
                cerr << "The time between interspersed substitutions does not match the timing of "
                        "events."
                     << endl;
                return false;
            }
        }
        if (!is_sorted(interspersed_substitutions.begin(), interspersed_substitutions.end(),
                       [](const Substitution &a, const Substitution &b) -> bool {
                           return a.time_event <= b.time_event;
                       })) {
            cerr << "The interspersed substitutions are not sorted." << endl;
            return false;
        }
        return true;
    };

    void intersperse_exon_substitutions() {
        // Find time of the next substitutions
        // At the same time sum the mut flow and sub flow
        while (true) {
            auto exon_next_sub =
                    min_element(exons.begin(), exons.end(), [](const Exon &a, const Exon &b) -> bool {
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
                sub.non_syn_sub_flow = 0;
                sub.non_syn_mut_flow = 0;
                sub.syn_mut_flow = 0;
                for (auto const &exon : exons) {
                    sub.non_syn_sub_flow += exon.substitutions.front().non_syn_sub_flow;
                    sub.non_syn_mut_flow += exon.substitutions.front().non_syn_mut_flow;
                    sub.syn_mut_flow += exon.substitutions.front().syn_mut_flow;
                }
                interspersed_substitutions.push_back(sub);
                exon_next_sub->substitutions.pop();
            }
        }
    };

    // Set the the DNA sequence to the mutation-selection equilibrium.
    void at_equilibrium(u_long nbr_rounds = 1) {
        // For all site of the sequence.
        for (u_long i = 0; i < nbr_rounds; i++) {
            for (auto &exon : exons) {
                for (u_long site{0}; site < exon.nbr_sites; site++) {
                    array<double, 64> codon_freqs =
                            codon_frequencies(exon.fitness_landscape.aa_selection_coefficients(exon.codon_seq, site), nuc_matrix, beta);
                    discrete_distribution<char> freq_codon_distr(
                            codon_freqs.begin(), codon_freqs.end());
                    char chosen_codon = freq_codon_distr(generator);
                    exon.codon_seq[site] = chosen_codon;
                }
            }
        }
    }

    void burn_in(int nbr_sub) {
        for (int i = 0; i < nbr_sub; i++) {
            for (auto &exon : exons) {
                exon.next_substitution(
                        nuc_matrix, beta, 0.0, numeric_limits<double>::infinity(), true);
            }
        }
        clear();
        cout << get_aa_str() << endl;
    }

    EVector delta_log_multivariate(double distance) const {
        EVector normal_vector = EVector::Zero(dimensions);
        for (int dim = 0; dim < dimensions; dim++) {
            normal_vector(dim) = normal_distrib(generator);
        }
        return sqrt(distance) * (transform_matrix * normal_vector);
    }

    void update_brownian(EVector const &delta) {
        log_multivariate += delta;

        beta = log_multivariate.beta();
        generation_time = log_multivariate.generation_time();
        nuc_matrix.set_mutation_rate(
                log_multivariate.mutation_rate_per_generation() / generation_time);
    }

    void run_forward(double t_max, Tree const &tree) {
        if (branch_wise) {
            EVector delta = delta_log_multivariate(t_max / tree.max_distance_to_root());
            update_brownian(delta / 2);
            piecewise_multivariate.AddMultivariate(t_max, log_multivariate);
            for (auto &exon : exons) { exon.run_substitutions(nuc_matrix, beta, 0.0, t_max); }
            intersperse_exon_substitutions();
            update_brownian(delta / 2);
            time_from_root += t_max;
            assert(check_consistency());
        } else {
            double time_current = 0.0;
            double step_in_year = time_grid_step;
            while (time_current < t_max) {
                if (time_current + step_in_year > t_max) { step_in_year = t_max - time_current; }
                update_brownian(delta_log_multivariate(step_in_year / tree.max_distance_to_root()));
                piecewise_multivariate.AddMultivariate(step_in_year, log_multivariate);
                for (auto &exon : exons) {
                    exon.run_substitutions(
                            nuc_matrix, beta, time_current, time_current + step_in_year);
                }
                intersperse_exon_substitutions();
                time_current += step_in_year;
                time_from_root += step_in_year;
                assert(check_consistency());
            }
        }
    }

    double predicted_dn_dn0(NucleotideRateMatrix const &rates, double const &relative_pop) const {
        double dn{0.}, dn0{0.};
        for (auto const &exon : exons) {
            double exon_dn{0}, exon_d0{0};
            tie(exon_dn, exon_d0) = exon.fitness_landscape.predicted_dn_dn0(exon.codon_seq, rates, relative_pop);
            dn += exon_dn;
            dn0 += exon_d0;
        }
        return dn / dn0;
    };

    double sequence_wise_predicted_dn_dn0(Sequence const &parent, NucleotideRateMatrix const &rates,
                                          double const &relative_pop) const {
        double dn{0.}, dn0{0.};

        for (size_t i = 0; i < exons.size(); i++) {
            double exon_dn{0}, exon_dn0{0};
            tie(exon_dn, exon_dn0) = exons[i].fitness_landscape.flow_dn_dn0(parent.exons[i].codon_seq, rates, relative_pop);
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
            cerr << "There is no synonymous substitutions, dN/dS can't be computed!" << endl;
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
            cerr << "There is no synonymous substitutions, dN/dS can't be computed!" << endl;
            return .0;
        } else {
            return (dn * ds0) / (dn0 * ds);
        }
    }

    void node_trace(string const &output_filename, Tree::NodeIndex node, Tree &tree,
                    Sequence const *parent) const {
        string node_name = tree.node_name(node);

        tree.set_tag(node, "population_size", d_to_string(beta));
        tree.set_tag(node, "generation_time", d_to_string(generation_time));
        tree.set_tag(node, "mutation_rate", d_to_string(nuc_matrix.mutation_rate));
        tree.set_tag(node, "mutation_rate_per_generation",
                     d_to_string(log_multivariate.mutation_rate_per_generation()));

        if (tree.is_root(node)) { return; }
        assert(parent != nullptr);
        double geom_pop_size = piecewise_multivariate.GeometricPopSize();
        piecewise_multivariate.add_to_tree(tree, node, geom_pop_size);

        tree.set_tag(node, "Branch_dNdN0_predicted",
                     d_to_string(predicted_dn_dn0(nuc_matrix, geom_pop_size)));
        tree.set_tag(node, "Branch_dNdN0_sequence_wise_predicted",
                     d_to_string(sequence_wise_predicted_dn_dn0(*parent, nuc_matrix, geom_pop_size)));

        double flow_dn_dn0 = flow_based_dn_dn0();
        double count_dn_dn0 = count_based_dn_dn0();
        double event_dn_ds = event_based_dn_ds();
        double count_dn_ds = count_based_dn_ds();
        dnd0_event_tot += flow_dn_dn0 * tree.node_length(node);
        dnd0_count_tot += count_dn_dn0 * tree.node_length(node);
        dnds_event_tot += event_dn_ds * tree.node_length(node);
        dnds_count_tot += count_dn_ds * tree.node_length(node);
        tree.set_tag(node, "Branch_dNdN0_flow_based", d_to_string(flow_dn_dn0));
        tree.set_tag(node, "Branch_dNdN0_count_based", d_to_string(count_dn_dn0));
        tree.set_tag(node, "Branch_dNdS_event_based", d_to_string(event_dn_ds));
        tree.set_tag(node, "Branch_dNdS_count_based", d_to_string(count_dn_ds));

        for (auto const &sub : interspersed_substitutions) {
            if (!sub.is_dummy()) {
                tracer_substitutions.add("NodeName", node_name);
                tracer_substitutions.add("Time", sub.time_event);
                tracer_substitutions.add("NucFrom", codonLexico.nucleotides[sub.n_from]);
                tracer_substitutions.add("NucTo", codonLexico.nucleotides[sub.n_to]);
                tracer_substitutions.add("CodonFrom", codonLexico.codon_string(sub.codon_from));
                tracer_substitutions.add("CodonTo", codonLexico.codon_string(sub.codon_to));
                tracer_substitutions.add("AAFrom", codonLexico.codon_aa_string(sub.codon_from));
                tracer_substitutions.add("AATo", codonLexico.codon_aa_string(sub.codon_to));
                tracer_substitutions.add("Site", sub.site);
            }
        }
        tracer_sequences.add("NodeName", node_name);
        tracer_sequences.add("CodonSequence", get_dna_str());
        tracer_sequences.add("AASequence", get_aa_str());

        if (!tree.is_leaf(node)) {
            tracer_fossils.add("NodeName", node_name);
            double age = tree.max_distance_to_root() - time_from_root;
            tracer_fossils.add("Age", age);
            tracer_fossils.add("LowerBound", age * 0.9);
            tracer_fossils.add("UpperBound", age * 1.1);
            return;
        }
        // If the node is a leaf, output the DNA sequence and name.
        write_sequence(output_filename, node_name, this->get_dna_str());

        tracer_traits.add("TaxonName", node_name);
        tracer_traits.add("LogGenerationTime", log_multivariate.log_generation_time());
        tracer_traits.add("LogPopulationSize", log_multivariate.log_population_size());
        tracer_traits.add(
                "LogMutationRatePerGeneration", log_multivariate.log_mutation_rate_per_generation());
    }

    // Method returning the DNA string corresponding to the codon sequence.
    string get_dna_str() const {
        string dna_str{};
        dna_str.reserve(nbr_nucleotides());

        // For each site of the sequence.
        for (auto const &exon : exons) {
            for (u_long site{0}; site < exon.nbr_sites; site++) {
                // Assert there is no stop in the sequence.
                assert(codonLexico.codon_to_aa[exon.codon_seq[site]] != 20);

                // Translate the site to a triplet of DNA nucleotides
                array<char, 3> triplet = codonLexico.codon_to_triplet[exon.codon_seq[site]];
                for (char position{0}; position < 3; position++) {
                    dna_str += codonLexico.nucleotides[triplet[position]];
                }
            }
        }
        return dna_str;  // return the DNA sequence as a string.
    }

    // Method returning the AA string corresponding to the codon sequence.
    string get_aa_str() const {
        string aa_str{};
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

    void write_matrices(string const &output_filename) {
        // For all site of the sequence.
        Trace trace;
        for (auto const &exon : exons) {
            for (u_long site{0}; site < exon.nbr_sites; site++) {
                // Codon original before substitution.
                trace.add("site", exon.position + site + 1);
                for (char codon_from{0}; codon_from < 64; codon_from++) {
                    // For all possible neighbors.
                    for (auto &neighbor : codonLexico.codon_to_neighbors[codon_from]) {
                        // Codon after mutation, Nucleotide original and Nucleotide after
                        // mutation.
                        char codon_to{0}, n_from{0}, n_to{0};
                        tie(codon_to, n_from, n_to) = neighbor;

                        if (codonLexico.codon_to_aa[codon_to] != 20) {
                            string key = "q_" + codonLexico.codon_string(codon_from) + "_" +
                                         codonLexico.codon_string(codon_to);
                            double s = exon.fitness_landscape.selection_coefficient(exon.codon_seq, site, codon_to);
                            double fix = nuc_matrix.normalized_rate(n_from, n_to) * Pfix(beta, s);
                            trace.add(key, fix);
                        }
                    }
                }
            }
        }
        trace.write_tsv(output_filename + ".matrices");
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
    vector<Sequence *> sequences;  // Vector of sequence DNA.

public:
    // Constructor
    Process(Tree &intree, Sequence &root_seq) : tree{intree}, sequences() {
        sequences.resize(tree.nb_nodes());
        sequences[tree.root()] = &root_seq;
    }

    void run(string const &output_filename) {
        sequences[tree.root()]->node_trace(output_filename, tree.root(), tree, nullptr);
        run_recursive(tree.root(), output_filename);
        ofstream nhx;
        nhx.open(output_filename + ".nhx");
        nhx << tree.as_string() << endl;
        nhx.close();
    }

    // Recursively iterate through the subtree.
    void run_recursive(Tree::NodeIndex node, string const &output_filename) {
        // Substitutions of the DNA sequence is generated.
        sequences.at(node)->clear();

        if (!tree.is_root(node)) {
            sequences.at(node)->run_forward(tree.node_length(node), tree);

            years_computed += tree.node_length(node);
            cout << years_computed << " years computed in total ("
                 << static_cast<int>(100 * years_computed / tree.total_length()) << "%) at node "
                 << tree.node_name(node) << " ("
                 << static_cast<int>(100 * tree.node_length(node) / tree.total_length()) << "%)."
                 << endl;

            sequences.at(node)->node_trace(
                    output_filename, node, tree, sequences[tree.parent(node)]);
        }

        // If the node is internal, iterate through the direct children.
        for (auto &child : tree.children(node)) {
            sequences.at(child) = new Sequence(*sequences[node]);
            run_recursive(child, output_filename);
        }
    }

    // Recursively iterate through the subtree and count the number of substitutions.
    tuple<long, long> nbr_substitutions() {
        long nbr_non_syn{0}, nbr_syn{0};

        for (auto const seq : sequences) {
            nbr_syn += count_if(seq->interspersed_substitutions.begin(),
                                seq->interspersed_substitutions.end(),
                                [](Substitution const &s) { return s.is_synonymous(); });
            nbr_non_syn += count_if(seq->interspersed_substitutions.begin(),
                                    seq->interspersed_substitutions.end(),
                                    [](Substitution const &s) { return s.is_non_synonymous(); });
        }

        return make_tuple(nbr_non_syn, nbr_syn);
    }
};

// Initialize static variables
double Process::years_computed = 0.0;