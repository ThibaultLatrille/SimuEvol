#pragma once

#include "fitness.hpp"
#include "random.hpp"

// composition
const std::array<double, 20> grantham_com = {0, 2.75, 1.38, 0.92, 0, 0.74, 0.58, 0, 0.33, 0, 0,
    1.33, 0.39, 0.89, 0.65, 1.42, 0.71, 0, 0.13, 0.2};

// polarity
const std::array<double, 20> grantham_pol = {8.1, 5.5, 13.0, 12.3, 5.2, 9.0, 10.4, 5.2, 11.3, 4.9,
    5.7, 11.6, 8.0, 10.5, 10.5, 9.2, 8.6, 5.9, 5.4, 6.2};

// volume
const std::array<double, 20> grantham_vol = {
    31, 55, 54, 83, 132, 3, 96, 111, 119, 111, 105, 56, 32.5, 85, 124, 32, 61, 84, 170, 136};

// associated weights
const double grantham_wcom = 1.833;  // 0.459 in Beaulieu et al
const double grantham_wpol = 0.1018;
const double grantham_wvol = 0.000399;

static double grantham_dist(char a, char b) {
    double tcom = grantham_com[b] - grantham_com[a];
    double tpol = grantham_pol[b] - grantham_pol[a];
    double tvol = grantham_vol[b] - grantham_vol[a];
    double d = sqrt(
        grantham_wcom * tcom * tcom + grantham_wpol * tpol * tpol + grantham_wvol * tvol * tvol);
    return d;
}

class StabilityLandscape final : public FitnessLandscape {
  private:
    std::vector<double> site_ddg;
    std::array<std::array<double, 20>, 20> DistanceMatrix{};

  public:
    std::vector<char> optimal_aa_seq;
    double alpha;
    double beta;
    double expression_level;
    bool grantham;

    explicit StabilityLandscape(int nbr_sites, double alpha, double gamma, double gamma_std,
        double beta, bool grantham, double expression_level)
        : site_ddg(nbr_sites, gamma),
          optimal_aa_seq(nbr_sites, 0),
          alpha{alpha},
          beta{beta},
          expression_level{expression_level},
          grantham{grantham} {
        additive_phenotype = true;
        auto bk_seed = generator();
        generator.seed(568);
        auto d = std::uniform_int_distribution<char>(0, 19);
        for (auto &site : optimal_aa_seq) { site = d(generator); }
        if (gamma_std != 0.0) {
            auto distr = std::gamma_distribution<double>(
                gamma * gamma / (gamma_std * gamma_std), gamma_std * gamma_std / gamma);
            for (auto &ddg : site_ddg) { ddg = distr(generator); }
        }
        generator.seed(bk_seed);
        if (grantham) {
            for (char i = 0; i < 20; ++i) {
                for (char j = 0; j < 20; ++j) {
                    DistanceMatrix[i][j] = grantham_dist(i, j);
                    DistanceMatrix[j][i] = DistanceMatrix[i][j];
                }
            }
            std::cout << "Grantham distance is activated" << std::endl;
        } else {
            std::cout << "Grantham distance is not activated" << std::endl;
        }
    }

    double ddg(size_t site, char codon_to) const {
        if (optimal_aa_seq[site] == codonLexico.codon_to_aa[codon_to]) { return 0.0; }
        auto out = site_ddg[site];
        if (grantham) {
            out *= DistanceMatrix[optimal_aa_seq[site]][codonLexico.codon_to_aa[codon_to]];
        }
        return out;
    }

    u_long nbr_sites() const override { return optimal_aa_seq.size(); }
};

class StabilityState final : public FitnessState {
  private:
    StabilityLandscape const &f;

    double x{0};
    double dG;

  public:
    std::unique_ptr<FitnessState> clone() const override {
        return std::make_unique<StabilityState>(*this);
    };

    explicit StabilityState(StabilityLandscape const &f) : FitnessState(f), f{f}, dG{f.alpha} {}

    bool operator==(FitnessState const &other) const override {
        auto p = dynamic_cast<StabilityState const *>(&other);
        return log_fitness == other.log_fitness and (x == p->x) and (dG == p->dG);
    }

    u_long nbr_sites() const override { return f.nbr_sites(); }

    double fitness(double deltaG) const {
        double factor = exp(-f.beta * deltaG);
        return (factor / (1 + factor));
    }

    double logfitness(double deltaG) const { return f.expression_level * log(fitness(deltaG)); }

    double selection_coefficient_ddG(double deltaG, double deltaGMutant) const {
        double fm = fitness(deltaGMutant), fr = fitness(deltaG);
        return f.expression_level * (fm - fr) / fr;
    }

    void update(std::vector<char> const &codon_seq, double const &pop_size) override {
        x = 0;
        dG = f.alpha;
        for (std::size_t site = 0; site < codon_seq.size(); ++site) {
            dG += f.ddg(site, codon_seq[site]);
            if (codonLexico.codon_to_aa.at(codon_seq.at(site)) != f.optimal_aa_seq.at(site)) {
                x += 1.0 / f.nbr_sites();
            }
        }
        log_fitness = logfitness(dG);
    }

    void update(std::vector<char> const &codon_seq, u_long site, char codon_to, bool burn_in,
        double const &pop_size) override {
        if (codonLexico.codon_to_aa.at(codon_seq.at(site)) == f.optimal_aa_seq.at(site)) {
            // We were optimal
            x += 1.0 / f.nbr_sites();
        } else if (codonLexico.codon_to_aa.at(codon_to) == f.optimal_aa_seq.at(site)) {
            // We were not optimal, and we are going optimal
            x -= 1.0 / f.nbr_sites();
        }
        dG += f.ddg(site, codon_to) - f.ddg(site, codon_seq.at(site));
        log_fitness = logfitness(dG);
        if (!burn_in) {
            summary_stats["sub-log-fitness"].add(log_fitness);
            summary_stats["sub-x"].add(x);
            summary_stats["sub-distance"].add(x * nbr_sites());
            summary_stats["sub-ΔG"].add(dG);
        }
    }

    double selection_coefficient(std::vector<char> const &codon_seq, u_long site, char codon_to,
        bool burn_in, double const &pop_size) const override {
        double mutantdG = dG;
        if (codonLexico.codon_to_aa[codon_seq[site]] == codonLexico.codon_to_aa[codon_to]) {
            return 0.0;
        }
        mutantdG += f.ddg(site, codon_to) - f.ddg(site, codon_seq.at(site));
        double s = selection_coefficient_ddG(dG, mutantdG);
        if (!burn_in) {
            summary_stats["mut-s"].add(s);
            summary_stats["mut-ΔG"].add(mutantdG);
            summary_stats["mut-ΔΔG"].add(mutantdG - dG);
            distribution_map["mut-ΔΔG"].add(mutantdG - dG);
        }
        return s;
    }
};

std::unordered_map<std::string, SummaryStatistic> FitnessState::summary_stats = {};

class StabilityArgParse {
  protected:
    TCLAP::CmdLine &cmd;

  public:
    explicit StabilityArgParse(TCLAP::CmdLine &cmd) : cmd{cmd} {}
    TCLAP::ValueArg<u_long> nbr_exons{
        "", "nbr_exons", "Number of exons in the protein", false, 10, "u_long", cmd};
    TCLAP::ValueArg<double> alpha{
        "", "alpha", "ΔG minimum for the optimal sequence", false, -118, "double", cmd};
    TCLAP::ValueArg<double> gamma{
        "", "gamma", "Increase in ΔG for destabilizying mutations", false, 1.0, "double", cmd};
    TCLAP::ValueArg<double> expression_level{
        "", "expression_level", "Expression level", false, 1.0, "double", cmd};
    TCLAP::ValueArg<double> gamma_std{
        "", "gamma_std", "Standard deviation for gamma", false, 0, "double", cmd};
    TCLAP::ValueArg<double> beta{"", "beta", "beta=1/kT in mol/kcal", false, 1.686, "double", cmd};
    TCLAP::SwitchArg grantham{"", "grantham", "Use grantham distance", cmd, false};

    void add_to_trace(Trace &trace) {
        trace.add("#nbr_exons", nbr_exons.getValue());
        trace.add("expression_level", expression_level.getValue());
        trace.add("alpha", alpha.getValue());
        trace.add("beta", beta.getValue());
        trace.add("gamma", gamma.getValue());
        trace.add("gamma_std", gamma_std.getValue());
        trace.add("grantham", grantham.getValue());
    }
};

class StabilityModel : public FitnessModel {
  public:
    StabilityModel(u_long exon_size, StabilityArgParse &args) : FitnessModel() {
        fitness_landscapes.reserve(args.nbr_exons.getValue());
        fitness_states.reserve(args.nbr_exons.getValue());
        StabilityLandscape landscape(exon_size, args.alpha.getValue(), args.gamma.getValue(),
            args.gamma_std.getValue(), args.beta.getValue(), args.grantham.getValue(),
            args.expression_level.getValue());

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