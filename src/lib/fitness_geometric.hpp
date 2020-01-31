#pragma once

#include "fitness.hpp"
#include "random.hpp"


double distance(std::vector<double> const &v) {
    return std::sqrt(std::accumulate(
        v.begin(), v.end(), 0.0, [](double a, double const &b) { return a + b * b; }));
}

std::vector<double> spherical_coord(u_long n, double radius) {
    std::vector<double> coord(n, 0);
    for (u_long i = 0; i < n; ++i) { coord[i] = normal_distrib(generator); }
    double norm = radius / distance(coord);
    for (u_long i = 0; i < n; ++i) { coord[i] *= norm; }
    return coord;
}

class GeometricLandscape final : public FitnessLandscape {
  public:
    // The fitness profiles of amino-acids.
    std::vector<std::array<std::vector<double>, 20>> sites_aa_phenotypes;
    u_long complexity;
    double peakness;
    double epistasis;
    double radius;
    double radius_std;

    GeometricLandscape(u_long nbr_sites, u_long const &complexity, double const &peakness,
        double const &epistasis, double const &radius, double const &radius_std)
        : complexity{complexity},
          peakness{peakness},
          epistasis{epistasis},
          radius{radius},
          radius_std{radius_std} {
        additive_phenotype = true;
        sites_aa_phenotypes.resize(nbr_sites);
        double r{radius};
        if (radius_std != 0.0) {
            r = std::gamma_distribution<double>((radius / radius_std) * (radius / radius_std),
                radius_std * radius_std / radius)(generator);
        }
        for (auto &site_aa_phenotypes : sites_aa_phenotypes) {
            for (auto &aa_phenotype : site_aa_phenotypes) {
                aa_phenotype = spherical_coord(complexity, r);
            }
        }
    }

    u_long nbr_sites() const override { return sites_aa_phenotypes.size(); }
};

class GeometricState final : public FitnessState {
  private:
    GeometricLandscape const &f;
    std::vector<double> phenotype;

  public:
    std::unique_ptr<FitnessState> clone() const override {
        return std::make_unique<GeometricState>(*this);
    };

    explicit GeometricState(GeometricLandscape const &f) : FitnessState(f), f{f} {
        phenotype.resize(f.complexity, 0);
    }

    bool operator==(FitnessState const &other) const override {
        return phenotype == dynamic_cast<GeometricState const *>(&other)->phenotype;
    };

    u_long nbr_sites() const override { return f.sites_aa_phenotypes.size(); }

    void update(std::vector<char> const &codon_seq, double const &pop_size) override {
        phenotype.resize(f.complexity, 0);
        for (u_long dim = 0; dim < f.complexity; ++dim) {
            for (u_long i = 0; i < nbr_sites(); ++i) {
                phenotype[dim] += f.sites_aa_phenotypes.at(i)
                                      .at(codonLexico.codon_to_aa.at(codon_seq.at(i)))
                                      .at(dim);
            }
        }
    }

    void update(std::vector<char> const &codon_seq, u_long site, char codon_to, bool burn_in,
        double const &pop_size) override {
        for (u_long dim = 0; dim < f.complexity; ++dim) {
            phenotype[dim] -= f.sites_aa_phenotypes.at(site)
                                  .at(codonLexico.codon_to_aa.at(codon_seq.at(site)))
                                  .at(dim);
            phenotype[dim] +=
                f.sites_aa_phenotypes.at(site).at(codonLexico.codon_to_aa.at(codon_to)).at(dim);
        }
        if (!burn_in) { summary_stats["sub-distance"].add(distance(phenotype)); }
    }

    double fitness(std::vector<double> const &v) const {
        return exp(-f.peakness * pow(distance(v), f.epistasis));
    }


    double selection_coefficient(std::vector<double> const &mutant_phenotype) const {
        double fp = fitness(phenotype);
        return (fitness(mutant_phenotype) - fp) / fp;
    }

    double selection_coefficient(FitnessState const &mutant, bool burn_in) const override {
        auto mutant_cast = dynamic_cast<GeometricState const *>(&mutant);
        double s = selection_coefficient(mutant_cast->phenotype);
        if (!burn_in) {
            summary_stats["mut-s"].add(s);
            summary_stats["mut-distance"].add(distance(mutant_cast->phenotype));
        }
        return s;
    };

    double selection_coefficient(std::vector<char> const &codon_seq, u_long site, char codon_to,
        bool burn_in, double const &pop_size) const override {
        auto mutant_phenotype = phenotype;
        for (u_long dim = 0; dim < f.complexity; ++dim) {
            mutant_phenotype[dim] -=
                f.sites_aa_phenotypes[site][codonLexico.codon_to_aa[codon_seq[site]]][dim];
            mutant_phenotype[dim] +=
                f.sites_aa_phenotypes[site][codonLexico.codon_to_aa[codon_to]][dim];
        }
        double s = selection_coefficient(mutant_phenotype);
        if (!burn_in) {
            summary_stats["mut-s"].add(s);
            summary_stats["mut-distance"].add(distance(mutant_phenotype));
        }
        return s;
    };
};

std::unordered_map<std::string, SummaryStatistic> FitnessState::summary_stats = {
    {"mut-s", SummaryStatistic()}, {"mut-distance", SummaryStatistic()},
    {"sub-distance", SummaryStatistic()}};

class GeometricArgParse {
  protected:
    TCLAP::CmdLine &cmd;

  public:
    explicit GeometricArgParse(TCLAP::CmdLine &cmd) : cmd{cmd} {}
    TCLAP::ValueArg<double> peakness{"", "peakness",
        "'alpha' parameter (peakness) of the fitness function (exp(-alpha*(d^beta))", false, 1e-3,
        "double", cmd};
    TCLAP::ValueArg<double> epistasis{"", "epistasis",
        "'beta' parameter (epistasis) of fitness function (exp(-alpha*(d^beta))", false, 2.0,
        "double", cmd};
    TCLAP::ValueArg<u_long> complexity{
        "", "complexity", "Complexity of the fitness landscape", false, 3, "u_long", cmd};
    TCLAP::ValueArg<double> radius{
        "", "radius", "Mean radius of mutations", false, 0.1, "double", cmd};
    TCLAP::ValueArg<double> radius_std{
        "", "radius_std", "Standard Deviation radius of mutations", false, 0.1, "double", cmd};
    TCLAP::ValueArg<u_long> nbr_exons{"", "nbr_exons", "Number of exons", false, 30, "u_long", cmd};

    void add_to_trace(Trace &trace) {
        assert(complexity.getValue() > 0);
        assert(peakness.getValue() > 0);
        assert(epistasis.getValue() > 0);
        assert(radius.getValue() > 0);
        assert(radius_std.getValue() >= 0);
        trace.add("radius", radius.getValue());
        trace.add("radius_std", radius_std.getValue());
        trace.add("peakness", peakness.getValue());
        trace.add("epistasis", epistasis.getValue());
        trace.add("#nbr_exons", nbr_exons.getValue());
        trace.add("complexity", complexity.getValue());
    }
};

class GeometricModel : public FitnessModel {
  public:
    GeometricModel(u_long &exon_size, GeometricArgParse &args) : FitnessModel() {
        fitness_landscapes.reserve(args.nbr_exons.getValue());
        fitness_states.reserve(args.nbr_exons.getValue());

        for (u_long exon{0}; exon < args.nbr_exons.getValue(); exon++) {
            fitness_landscapes.emplace_back(std::make_unique<GeometricLandscape>(exon_size,
                args.complexity.getValue(), args.peakness.getValue(), args.epistasis.getValue(),
                args.radius.getValue(), args.radius_std.getValue()));
            fitness_states.emplace_back(std::make_unique<GeometricState>(
                *dynamic_cast<GeometricLandscape *>(fitness_landscapes.at(exon).get())));
        }
        assert(nbr_sites() == exon_size * args.nbr_exons.getValue());
        std::cout << fitness_landscapes.size() << " exons created." << std::endl;
    }
};