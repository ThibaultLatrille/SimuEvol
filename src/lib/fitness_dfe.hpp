#pragma once

#include <array>
#include <vector>
#include "fitness.hpp"
#include "random.hpp"
#include "tclap/CmdLine.h"

class DfeParameters final : public FitnessLandscape {
  public:
    u_long exon_size;
    bool reflected;
    double mean;
    double shape;
    double scale;

    explicit DfeParameters(u_long exon_size, bool reflected, double mean, double shape)
        : exon_size{exon_size}, reflected{reflected}, mean{mean}, shape{shape} {
        additive_phenotype = true;
        if (shape != 0.0) { scale = mean / shape; }
    }

    double selection_coefficient(u_long site, char aa_to, double const &pop_size) const {
        auto bk_seed = generator();
        generator.seed(site * 20 + aa_to);
        double s{-mean};
        if (shape != 0.0) { s = -std::gamma_distribution<double>(shape, scale)(generator); }
        if (reflected) {
            // Estimating the distribution of fitness effects from DNA sequence data: Implications
            // for the molecular clock. Gwenaël Piganeau and Adam Eyre-Walker, PNAS, 2003.
            // https://doi.org/10.1073/pnas.1833064100
            double S = -4 * pop_size * s;
            double x = exp(S) / (exp(S) + 1);
            assert(0.5 <= x && x <= 1.0);
            double draw = std::uniform_real_distribution<double>(0.0, 1.0)(generator);
            if (draw > x) {
                s = -s;
                assert(s >= 0);
            }
        }
        generator.seed(bk_seed);
        return s;
    }

    u_long nbr_sites() const override { return exon_size; }
};

class DfeState final : public FitnessState {
  private:
    DfeParameters const &f;

  public:
    std::unique_ptr<FitnessState> clone() const override {
        return std::make_unique<DfeState>(*this);
    };

    explicit DfeState(DfeParameters const &f) : FitnessState(f), f{f} {}

    bool operator==(FitnessState const &other) const override {
        return log_fitness == other.log_fitness;
    };

    u_long nbr_sites() const override { return f.nbr_sites(); }

    void update(std::vector<char> const &codon_seq, double const &pop_size) override {
        log_fitness = 0;
        for (u_long i = 0; i < nbr_sites(); ++i) {
            log_fitness +=
                f.selection_coefficient(i, codonLexico.codon_to_aa[codon_seq.at(i)], pop_size);
        }
    }

    void update(std::vector<char> const &codon_seq, u_long site, char codon_to, bool burn_in,
        double const &pop_size) override {
        log_fitness += f.selection_coefficient(site, codonLexico.codon_to_aa[codon_to], pop_size);
        if (!burn_in) { summary_stats["sub-log-fitness"].add(log_fitness); }
    }

    double selection_coefficient(std::vector<char> const &codon_seq, u_long site, char codon_to,
        bool burn_in, double const &pop_size) const override {
        double s = f.selection_coefficient(site, codonLexico.codon_to_aa[codon_to], pop_size);
        if (!burn_in) {
            summary_stats["mut-s"].add(s);
            summary_stats["mut-log-fitness"].add(log_fitness + s);
        }
        return s;
    };

    std::array<double, 20> aa_selection_coefficients(
        std::vector<char> const &codon_seq, u_long site, double const &pop_size) const override {
        std::array<double, 20> profiles{};
        for (char aa = 0; aa < 20; ++aa) {
            profiles[aa] = f.selection_coefficient(site, aa, pop_size);
        }
        return profiles;
    }
};

std::unordered_map<std::string, SummaryStatistic> FitnessState::summary_stats = {
    {"mut-s", SummaryStatistic()}, {"mut-log-fitness", SummaryStatistic()},
    {"sub-log-fitness", SummaryStatistic()}};

class DfeArgParse {
  protected:
    TCLAP::CmdLine &cmd;

  public:
    explicit DfeArgParse(TCLAP::CmdLine &cmd) : cmd{cmd} {}
    TCLAP::SwitchArg gamma_reflected{
        "", "gamma_reflected", "True if reflected gamma.", cmd, false};
    TCLAP::ValueArg<double> gamma_mean{
        "", "gamma_mean", "Mean value of selection coefficient effect", false, 1.0, "double", cmd};
    TCLAP::ValueArg<double> gamma_shape{
        "", "gamma_shape", "Shape of the gamma distribution", false, 1.0, "double", cmd};
    TCLAP::ValueArg<u_long> nbr_exons{"", "nbr_exons", "Number of exons", false, 30, "u_long", cmd};


    void add_to_trace(Trace &trace) {
        assert(gamma_mean.getValue() >= 0);
        assert(gamma_shape.getValue() >= 0);
        trace.add("#nbr_exons", nbr_exons.getValue());
        trace.add("gamma_distribution_reflected", gamma_reflected.getValue());
        trace.add("gamma_distribution_mean", gamma_mean.getValue());
        trace.add("gamma_distribution_shape", gamma_shape.getValue());
    }
};

class SequenceDfeModel : public FitnessModel {
  public:
    SequenceDfeModel(u_long &exon_size, DfeArgParse &args) : FitnessModel() {
        fitness_landscapes.reserve(args.nbr_exons.getValue());
        fitness_states.reserve(args.nbr_exons.getValue());
        for (u_long exon{0}; exon < args.nbr_exons.getValue(); exon++) {
            fitness_landscapes.emplace_back(
                std::make_unique<DfeParameters>(exon_size, args.gamma_reflected.getValue(),
                    args.gamma_mean.getValue(), args.gamma_shape.getValue()));
            fitness_states.emplace_back(std::make_unique<DfeState>(
                *dynamic_cast<DfeParameters *>(fitness_landscapes.at(exon).get())));
        }
        assert(nbr_sites() == exon_size * args.nbr_exons.getValue());
        std::cout << fitness_landscapes.size() << " exons created." << std::endl;
    }
};
