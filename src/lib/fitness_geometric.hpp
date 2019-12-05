#pragma once

#include "fitness.hpp"
#include "random.hpp"


double distance(std::vector<double> const &v) {
    return std::sqrt(std::accumulate(
        v.begin(), v.end(), 0.0, [](double a, double const &b) { return a + b * b; }));
}

double fitness(std::vector<double> const &v, double a = 0.5, int q = 2) {
    return exp(-a * pow(distance(v), q));
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

    GeometricLandscape(u_long nbr_sites, u_long const &complexity) : complexity{complexity} {
        sites_aa_phenotypes.resize(nbr_sites);
        for (auto &site_aa_phenotypes : sites_aa_phenotypes) {
            for (auto &aa_phenotype : site_aa_phenotypes) {
                aa_phenotype = spherical_coord(complexity, 0.1);
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

    explicit GeometricState(GeometricLandscape const &f) : f{f} {
        phenotype.resize(f.complexity, 0);
    }

    bool operator==(FitnessState const &other) const override {
        return phenotype == dynamic_cast<GeometricState const *>(&other)->phenotype;
    };

    u_long nbr_sites() const override { return f.sites_aa_phenotypes.size(); }

    void update(std::vector<char> const &codon_seq) override {
        phenotype.resize(f.complexity, 0);
        for (u_long dim = 0; dim < f.complexity; ++dim) {
            for (u_long i = 0; i < nbr_sites(); ++i) {
                phenotype[dim] += f.sites_aa_phenotypes.at(i)
                                      .at(codonLexico.codon_to_aa.at(codon_seq.at(i)))
                                      .at(dim);
            }
        }
    }

    void update(
        std::vector<char> const &codon_seq, u_long site, char codon_to, bool burn_in) override {
        for (u_long dim = 0; dim < f.complexity; ++dim) {
            phenotype[dim] -= f.sites_aa_phenotypes.at(site)
                                  .at(codonLexico.codon_to_aa.at(codon_seq.at(site)))
                                  .at(dim);
            phenotype[dim] +=
                f.sites_aa_phenotypes.at(site).at(codonLexico.codon_to_aa.at(codon_to)).at(dim);
        }
        if (!burn_in) { summary_stats["sub-distance"].add(distance(phenotype)); }
    }

    double selection_coefficient(std::vector<char> const &codon_seq, u_long site, char codon_to,
        bool burn_in) const override {
        auto mutant_phenotype = phenotype;
        for (u_long dim = 0; dim < f.complexity; ++dim) {
            mutant_phenotype[dim] -=
                f.sites_aa_phenotypes[site][codonLexico.codon_to_aa[codon_seq[site]]][dim];
            mutant_phenotype[dim] +=
                f.sites_aa_phenotypes[site][codonLexico.codon_to_aa[codon_to]][dim];
        }
        double fp = fitness(phenotype);
        double fm = fitness(mutant_phenotype);
        double s = (fm - fp) / fp;
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

class GeometricModel : public FitnessModel {
  public:
    GeometricModel(u_long &exon_size, u_long &nbr_exons, double complexity) : FitnessModel() {
        fitness_landscapes.reserve(nbr_exons);
        fitness_states.reserve(nbr_exons);

        for (u_long exon{0}; exon < nbr_exons; exon++) {
            fitness_landscapes.emplace_back(
                std::make_unique<GeometricLandscape>(exon_size, complexity));
            fitness_states.emplace_back(std::make_unique<GeometricState>(
                *dynamic_cast<GeometricLandscape *>(fitness_landscapes.at(exon).get())));
        }
        assert(nbr_sites() == exon_size * nbr_exons);
        std::cout << fitness_landscapes.size() << " exons created." << std::endl;
    }
};