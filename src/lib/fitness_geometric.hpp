#pragma once

#include "fitness.hpp"
#include "random.hpp"


double distance(std::vector<double> const &v) {
    return std::sqrt(
            std::accumulate(v.begin(), v.end(), 0.0, [](double a, double const &b) { return a + b * b; }));
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

    explicit GeometricLandscape(u_long nbr_sites, const u_long &complexity) :
            complexity{complexity} {
        sites_aa_phenotypes.resize(nbr_sites);
        for (auto &site_aa_phenotypes : sites_aa_phenotypes) {
            for (auto &aa_phenotype : site_aa_phenotypes) {
                aa_phenotype = spherical_coord(complexity, 0.1);
            }
        }
    }

    u_long nbr_sites() const override { return sites_aa_phenotypes.size(); }

    double selection_coefficient(const std::vector<char> &codon_seq, u_long site, char codon_to) const override {
        std::vector<double> phenotype(complexity, 0);
        for (u_long dim = 0; dim < complexity; ++dim) {
            for (u_long i = 0; i < nbr_sites(); ++i) {
                phenotype[dim] += sites_aa_phenotypes.at(i)
                        .at(codonLexico.codon_to_aa.at(codon_seq.at(i)))
                        .at(dim);
            }
        }
        auto mutant_phenotype = phenotype;
        for (u_long dim = 0; dim < complexity; ++dim) {
            mutant_phenotype[dim] -=
                    sites_aa_phenotypes[site][codonLexico.codon_to_aa[codon_seq[site]]][dim];
            mutant_phenotype[dim] += sites_aa_phenotypes[site][codonLexico.codon_to_aa[codon_to]][dim];
        }
        double f = fitness(phenotype);
        double f1 = fitness(mutant_phenotype);
        return (f1 - f) / f;
    };

    std::array<double, 20> aa_selection_coefficients(const std::vector<char> &codon_seq, u_long site) const override {
        std::array<double, 20> aa_sel_coeffs{};
        for (char aa = 0; aa < 20; aa++) {
            auto it = std::find(
                    codonLexico.codon_to_aa.begin(), codonLexico.codon_to_aa.end(), aa);
            assert(it != codonLexico.codon_to_aa.end());
            char codon_to = std::distance(codonLexico.codon_to_aa.begin(), it);
            assert(codonLexico.codon_to_aa[codon_to] == aa);
            aa_sel_coeffs[aa] = selection_coefficient(codon_seq, site, codon_to);
        }
        return aa_sel_coeffs;
    }
};

class SequenceGeometricLandscape : public SequenceFitnessLandscape {
public:

    SequenceGeometricLandscape(u_long &exon_size, u_long &nbr_exons, double complexity) : SequenceFitnessLandscape() {
        this->reserve(nbr_exons);
        for (u_long exon{0}; exon < nbr_exons; exon++) {
            this->emplace_back(std::make_unique<GeometricLandscape>(exon_size, complexity));
        }
        assert(nbr_sites() == exon_size * nbr_exons);
        std::cout << this->size() << " exons created." << std::endl;
    }
};