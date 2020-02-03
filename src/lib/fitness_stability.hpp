#pragma once

#include "fitness.hpp"
#include "random.hpp"

class StabilityLandscape final : public FitnessLandscape {
  public:
    std::vector<char> optimal_aa_seq;
    double dG_min{-10};
    double ddG{0.1};

    explicit StabilityLandscape(int nbr_sites, double dG_min, double ddG)
        : optimal_aa_seq(nbr_sites, 0), dG_min{dG_min}, ddG{ddG} {
        additive_phenotype = true;
        auto d = std::uniform_int_distribution<char>(0, 19);
        for (auto &site : optimal_aa_seq) { site = d(generator); }
    }

    u_long nbr_sites() const override { return optimal_aa_seq.size(); }
};


double fitness(double deltaG) {
    double factor = exp(-deltaG);
    return (factor / (1 + factor));
}

double logfitness(double deltaG) { return log(fitness(deltaG)); }

double selection_coefficient_ddG(double deltaG, double deltaGMutant) {
    double fm = fitness(deltaGMutant), f = fitness(deltaG);
    return (fm - f) / f;
}

class StabilityState final : public FitnessState {
  private:
    StabilityLandscape const &f;

    unsigned distance{0};
    double dG;

  public:
    std::unique_ptr<FitnessState> clone() const override {
        return std::make_unique<StabilityState>(*this);
    };

    explicit StabilityState(StabilityLandscape const &f) : FitnessState(f), f{f}, dG{f.dG_min} {}

    bool operator==(FitnessState const &other) const override {
        auto p = dynamic_cast<StabilityState const *>(&other);
        return log_fitness == other.log_fitness and (distance == p->distance) and (dG == p->dG);
    }

    u_long nbr_sites() const override { return f.nbr_sites(); }

    void update(std::vector<char> const &codon_seq, double const &pop_size) override {
        distance = 0;
        dG = f.dG_min;
        for (size_t site = 0; site < codon_seq.size(); ++site) {
            if (codonLexico.codon_to_aa.at(codon_seq.at(site)) != f.optimal_aa_seq.at(site)) {
                distance++;
                dG += f.ddG;
            }
        }
        log_fitness = logfitness(dG);
    }

    void update(std::vector<char> const &codon_seq, u_long site, char codon_to, bool burn_in,
        double const &pop_size) override {
        if (codonLexico.codon_to_aa.at(codon_seq.at(site)) == f.optimal_aa_seq.at(site)) {
            // We were optimal
            distance++;
            dG += f.ddG;
        } else if (codonLexico.codon_to_aa.at(codon_to) == f.optimal_aa_seq.at(site)) {
            // We were not optimal, and we are going optimal
            distance--;
            dG -= f.ddG;
        }
        log_fitness = logfitness(dG);
        if (!burn_in) {
            summary_stats["sub-log-fitness"].add(log_fitness);
            summary_stats["sub-distance/#sites"].add(static_cast<double>(distance) / nbr_sites());
            summary_stats["sub-distance"].add(distance);
            summary_stats["sub-ΔG"].add(dG);
        }
    }

    double selection_coefficient(std::vector<char> const &codon_seq, u_long site, char codon_to,
        bool burn_in, double const &pop_size) const override {
        double mutantdG = dG;
        if (codonLexico.codon_to_aa[codon_seq[site]] == codonLexico.codon_to_aa[codon_to]) {
            return 0.0;
        }
        if (codonLexico.codon_to_aa.at(codon_seq.at(site)) == f.optimal_aa_seq.at(site)) {
            // We were optimal
            mutantdG += f.ddG;
        } else if (codonLexico.codon_to_aa.at(codon_to) == f.optimal_aa_seq.at(site)) {
            // We were not optimal, and we are going optimal
            mutantdG -= f.ddG;
        }
        double s = selection_coefficient_ddG(dG, mutantdG);
        if (!burn_in) {
            summary_stats["mut-s"].add(s);
            summary_stats["mut-ΔG"].add(mutantdG);
            summary_stats["mut-ΔΔG"].add(mutantdG - dG);
        }
        return s;
    }
};

std::unordered_map<std::string, SummaryStatistic> FitnessState::summary_stats = {
    {"mut-s", SummaryStatistic()}, {"mut-ΔG", SummaryStatistic()}, {"mut-ΔΔG", SummaryStatistic()},
    {"sub-ΔG", SummaryStatistic()}, {"sub-log-fitness", SummaryStatistic()},
    {"sub-distance/#sites", SummaryStatistic()}, {"sub-distance", SummaryStatistic()}};

class StabilityArgParse {
  protected:
    TCLAP::CmdLine &cmd;

  public:
    explicit StabilityArgParse(TCLAP::CmdLine &cmd) : cmd{cmd} {}
    TCLAP::ValueArg<u_long> nbr_exons{
        "", "nbr_exons", "Number of exons in the protein", false, 5000, "u_long", cmd};
    TCLAP::ValueArg<double> dg_min{
        "", "dg_min", "ΔG minimum for the optimal sequence", false, -10, "double", cmd};
    TCLAP::ValueArg<double> ddg{"", "ddg",
        "ΔΔG for each mutation going away from the optimal sequence", false, 0.1, "double", cmd};
    void add_to_trace(Trace &trace) {
        trace.add("#nbr_exons", nbr_exons.getValue());
        trace.add("ΔG-minimum", dg_min.getValue());
        trace.add("ΔΔG-mutation", ddg.getValue());
    }
};

class StabilityModel : public FitnessModel {
  public:
    StabilityModel(u_long exon_size, StabilityArgParse &args) : FitnessModel() {
        fitness_landscapes.reserve(args.nbr_exons.getValue());
        fitness_states.reserve(args.nbr_exons.getValue());
        StabilityLandscape landscape(exon_size, args.dg_min.getValue(), args.ddg.getValue());

        for (u_long exon{0}; exon < args.nbr_exons.getValue(); exon++) {
            fitness_landscapes.emplace_back(std::make_unique<StabilityLandscape>(landscape));
            fitness_states.emplace_back(std::make_unique<StabilityState>(
                *dynamic_cast<StabilityLandscape *>(fitness_landscapes.at(exon).get())));
        }
        assert(nbr_sites() == exon_size * args.nbr_exons.getValue());
        std::cout << fitness_landscapes.size() << " exons created." << std::endl;
    }

    std::string aa_seq() {
        std::string s{};
        for (u_long exon{0}; exon < nbr_exons(); exon++) {
            auto tmp = dynamic_cast<StabilityLandscape *>(fitness_landscapes.at(exon).get())
                           ->optimal_aa_seq;
            std::for_each(tmp.begin(), tmp.end(),
                [](char &c) -> char { return c = codonLexico.amino_acids.at(c); });
            s += std::string(tmp.begin(), tmp.end());
        }
        return s;
    }
};