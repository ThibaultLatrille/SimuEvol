#pragma once

#include "fitness.hpp"

double LOG_UNFOLDED_STATES = 368.413615;  // log(1.0E160)
double TEMPERATURE = 0.6;

std::array<std::string, 2> nativeProtein = {"1qhw", "A"};
std::vector<std::array<std::string, 2>> unfoldedProteinList = {{"1mty", "B"}, {"1n00", "A"},
    {"1gwu", "A"}, {"1kwf", "A"}, {"1iom", "A"}, {"1jfb", "A"}, {"1wer", " "}, {"1hz4", "A"},
    {"1uby", " "}, {"1t5j", "A"}, {"1dmh", "A"}, {"2bbv", "A"}, {"1ojj", "A"}, {"1nsz", "A"},
    {"1mkf", "A"}, {"1jj2", "B"}, {"1wkr", "A"}, {"1gyh", "A"}, {"3sil", " "}, {"1pby", "B"},
    {"1o88", "A"}, {"1odm", "A"}, {"1jub", "A"}, {"1ek6", "A"}, {"1oc7", "A"}, {"1jl5", "A"},
    {"1esd", " "}, {"1jil", "A"}, {"1t5o", "A"}, {"1umd", "A"}, {"1svm", "A"}, {"1l5o", "A"},
    {"1ga6", "A"}, {"1woh", "A"}, {"1wch", "A"}, {"1m4l", "A"}, {"1nd6", "A"}, {"1i4w", "A"},
    {"1o4s", "A"}, {"1jkm", "A"}, {"2mas", "A"}, {"1rkd", " "}, {"1e19", "A"}, {"1cnz", "A"},
    {"1qop", "B"}, {"1moq", " "}, {"1v6s", "A"}, {"1jix", "A"}, {"1o7j", "A"}, {"1pfk", "A"},
    {"1ir6", "A"}, {"1to6", "A"}, {"1qo0", "A"}, {"1sbp", " "}, {"1nbf", "A"}};

std::vector<std::vector<double>> pdbInteractionMatrix = {
    {-0.13, 0.43, 0.28, 0.12, 0.00, 0.08, 0.26, -0.07, 0.34, -0.22, -0.01, 0.14, 0.25, 0.03, 0.10,
        -0.06, -0.09, -0.09, 0.09, -0.10, 0.00},
    {0.43, 0.11, -0.14, -0.72, 0.24, -0.52, -0.74, -0.04, -0.12, 0.42, 0.35, 0.75, 0.31, 0.41,
        -0.38, 0.17, -0.35, -0.16, -0.25, 0.30, 0.00},
    {0.28, -0.14, -0.53, -0.30, 0.13, -0.25, -0.32, -0.14, -0.24, 0.53, 0.30, -0.33, 0.08, 0.18,
        -0.18, -0.14, -0.11, 0.06, -0.20, 0.50, 0.00},
    {0.12, -0.72, -0.30, 0.04, 0.03, -0.17, -0.15, -0.22, -0.39, 0.59, 0.67, -0.76, 0.65, 0.39,
        0.04, -0.31, -0.29, 0.24, 0.00, 0.58, 0.00},
    {0.00, 0.24, 0.13, 0.03, -0.20, 0.05, 0.69, -0.08, -0.19, 0.16, -0.08, 0.71, 0.19, -0.23, 0.00,
        -0.02, 0.19, 0.08, 0.04, 0.06, 0.00},
    {0.08, -0.52, -0.25, -0.17, 0.05, 0.29, -0.17, -0.06, -0.02, 0.36, 0.26, -0.38, 0.46, 0.49,
        -0.42, -0.14, -0.14, 0.08, -0.20, 0.24, 0.00},
    {0.26, -0.74, -0.32, -0.15, 0.69, -0.17, -0.03, 0.25, -0.45, 0.35, 0.43, -0.97, 0.44, 0.27,
        -0.10, -0.26, 0.00, 0.29, -0.10, 0.34, 0.00},
    {-0.07, -0.04, -0.14, -0.22, -0.08, -0.06, 0.25, -0.38, 0.20, 0.25, 0.23, 0.11, 0.19, 0.38,
        -0.11, -0.16, -0.26, 0.18, 0.14, 0.16, 0.00},
    {0.34, -0.12, -0.24, -0.39, -0.19, -0.02, -0.45, 0.20, -0.29, 0.49, 0.16, 0.22, 0.99, -0.16,
        -0.21, -0.05, -0.19, -0.12, -0.34, 0.19, 0.00},
    {-0.22, 0.42, 0.53, 0.59, 0.16, 0.36, 0.35, 0.25, 0.49, -0.22, -0.41, 0.36, -0.28, -0.19, 0.25,
        0.21, 0.14, 0.02, 0.11, -0.25, 0.00},
    {-0.01, 0.35, 0.30, 0.67, -0.08, 0.26, 0.43, 0.23, 0.16, -0.41, -0.27, 0.19, -0.20, -0.30, 0.42,
        0.25, 0.20, -0.09, 0.24, -0.29, 0.00},
    {0.14, 0.75, -0.33, -0.76, 0.71, -0.38, -0.97, 0.11, 0.22, 0.36, 0.19, 0.25, 0.00, 0.44, 0.11,
        -0.13, -0.09, 0.22, -0.21, 0.44, 0.00},
    {0.25, 0.31, 0.08, 0.65, 0.19, 0.46, 0.44, 0.19, 0.99, -0.28, -0.20, 0.00, 0.04, -0.42, -0.34,
        0.14, 0.19, -0.67, -0.13, -0.14, 0.00},
    {0.03, 0.41, 0.18, 0.39, -0.23, 0.49, 0.27, 0.38, -0.16, -0.19, -0.30, 0.44, -0.42, -0.44, 0.20,
        0.29, 0.31, -0.16, 0.00, -0.22, 0.00},
    {0.10, -0.38, -0.18, 0.04, 0.00, -0.42, -0.10, -0.11, -0.21, 0.25, 0.42, 0.11, -0.34, 0.20,
        0.26, 0.01, -0.07, -0.28, -0.33, 0.09, 0.00},
    {-0.06, 0.17, -0.14, -0.31, -0.02, -0.14, -0.26, -0.16, -0.05, 0.21, 0.25, -0.13, 0.14, 0.29,
        0.01, -0.20, -0.08, 0.34, 0.09, 0.18, 0.00},
    {-0.09, -0.35, -0.11, -0.29, 0.19, -0.14, 0.00, -0.26, -0.19, 0.14, 0.20, -0.09, 0.19, 0.31,
        -0.07, -0.08, 0.03, 0.22, 0.13, 0.25, 0.00},
    {-0.09, -0.16, 0.06, 0.24, 0.08, 0.08, 0.29, 0.18, -0.12, 0.02, -0.09, 0.22, -0.67, -0.16,
        -0.28, 0.34, 0.22, -0.12, -0.04, -0.07, 0.00},
    {0.09, -0.25, -0.20, 0.00, 0.04, -0.20, -0.10, 0.14, -0.34, 0.11, 0.24, -0.21, -0.13, 0.00,
        -0.33, 0.09, 0.13, -0.04, -0.06, 0.02, 0.00},
    {-0.10, 0.30, 0.50, 0.58, 0.06, 0.24, 0.34, 0.16, 0.19, -0.25, -0.29, 0.44, -0.14, -0.22, 0.09,
        0.18, 0.25, -0.07, 0.02, -0.29, 0.00},
    {0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00,
        0.00, 0.00, 0.00, 0.00, 0.00}};

Codon codonPDB = Codon("ARNDCQEGHILKMFPSTWYVX");

double getCodonInteraction(char codonA, char codonB) {
    return pdbInteractionMatrix[codonPDB.codon_to_aa[codonA]][codonPDB.codon_to_aa[codonB]];
}

void removeSpaces(std::string &str) {
    str.erase(remove_if(str.begin(), str.end(), ::isspace), str.end());
}

class Structure {
  private:
    std::vector<std::vector<double>> distanceMatrix;
    std::vector<std::array<double, 3>> alphaCoordinateVector;
    std::vector<std::array<double, 3>> betaCoordinateVector;
    std::string proteinName;
    std::string chainName;
    size_t size;
    double cut_off;
    std::vector<std::array<size_t, 2>> contactVector;
    std::vector<std::vector<size_t>> siteContactVector;


    /**
     * Make a list of all contacts made between sites in the structure
     */
    void makeContactVector() {
        distanceMatrix.resize(size, std::vector<double>(size, 0.0));
        siteContactVector.resize(size, {});
        for (size_t iRes = 0; iRes < size - 1; iRes++) {
            std::array<double, 3> iCoord = betaCoordinateVector.at(iRes);
            for (size_t jRes = iRes + 1; jRes < size; jRes++) {
                std::array<double, 3> jCoord = betaCoordinateVector.at(jRes);
                distanceMatrix[iRes][jRes] =
                    sqrt(((iCoord[0] - jCoord[0]) * (iCoord[0] - jCoord[0]) +
                          (iCoord[1] - jCoord[1]) * (iCoord[1] - jCoord[1]) +
                          (iCoord[2] - jCoord[2]) * (iCoord[2] - jCoord[2])));
                distanceMatrix[jRes][iRes] = distanceMatrix.at(iRes).at(jRes);
            }
        }
        for (size_t iRes = 0; iRes < size - 2; iRes++) {
            for (size_t jRes = iRes + 2; jRes < size; jRes++) {
                if (distanceMatrix.at(iRes).at(jRes) < cut_off) {
                    std::array<size_t, 2> contact{};
                    contact[0] = iRes;
                    contact[1] = jRes;
                    contactVector.push_back(contact);
                    siteContactVector.at(iRes).push_back(jRes);
                    siteContactVector.at(jRes).push_back(iRes);
                }
            }
        }
    }

    /**
     * Check for obvious problems in protein structure
     */
    bool allOK() const {
        bool ok = true;
        for (size_t iRes = 0; iRes < alphaCoordinateVector.size() - 1; iRes++) {
            size_t jRes = iRes + 1;
            std::array<double, 3> iCoord = alphaCoordinateVector.at(iRes);
            std::array<double, 3> jCoord = alphaCoordinateVector.at(jRes);
            double dist2 = ((iCoord[0] - jCoord[0]) * (iCoord[0] - jCoord[0]) +
                            (iCoord[1] - jCoord[1]) * (iCoord[1] - jCoord[1]) +
                            (iCoord[2] - jCoord[2]) * (iCoord[2] - jCoord[2]));
            if (dist2 > (4.25 * 4.25)) {
                std::cout << "Discontinuity in protein " << proteinName << chainName << std::endl;
                std::cout << "  Residues " << iRes << " and " << jRes << "\t" << sqrt(dist2)
                          << std::endl;
                ok = false;
            }
        }

        double rg = 0.0;
        for (size_t iRes = 0; iRes < alphaCoordinateVector.size(); iRes++) {
            for (size_t jRes = 0; jRes < alphaCoordinateVector.size(); jRes++) {
                std::array<double, 3> iCoord = alphaCoordinateVector.at(iRes);
                std::array<double, 3> jCoord = alphaCoordinateVector.at(jRes);
                double dist2 = ((iCoord[0] - jCoord[0]) * (iCoord[0] - jCoord[0]) +
                                (iCoord[1] - jCoord[1]) * (iCoord[1] - jCoord[1]) +
                                (iCoord[2] - jCoord[2]) * (iCoord[2] - jCoord[2]));
                rg += dist2;
            }
        }
        rg /= (2.0 * alphaCoordinateVector.size() * alphaCoordinateVector.size());
        if (rg > 650) {
            std::cout << "Problem with rg " << proteinName << chainName + "\t"
                      << alphaCoordinateVector.size() << "\t" << rg << std::endl;
            ok = false;
        }

        double nContacts = 0.0;
        for (size_t iRes = 0; iRes < betaCoordinateVector.size(); iRes++) {
            for (size_t jRes = 0; jRes < betaCoordinateVector.size(); jRes++) {
                std::array<double, 3> iCoord = betaCoordinateVector.at(iRes);
                std::array<double, 3> jCoord = betaCoordinateVector.at(jRes);
                double dist2 = ((iCoord[0] - jCoord[0]) * (iCoord[0] - jCoord[0]) +
                                (iCoord[1] - jCoord[1]) * (iCoord[1] - jCoord[1]) +
                                (iCoord[2] - jCoord[2]) * (iCoord[2] - jCoord[2]));
                if (dist2 < 36.0) { nContacts++; }
            }
        }

        nContacts /= betaCoordinateVector.size();
        if ((nContacts < 5) || (nContacts > 7)) {
            std::cout << "Problem with nContacts " << proteinName << chainName << "\t"
                      << alphaCoordinateVector.size() << "\t" << nContacts << std::endl;
            ok = false;
        }
        return ok;
    }

    /**
     * Read in data from pdb file
     */
    void readPDBFile(std::string const &pdb_folder) {
        std::ifstream input_stream(pdb_folder + "/" + proteinName + ".pdb");
        if (!input_stream) {
            std::cerr << "PDB file " << pdb_folder << "/" << proteinName << ".pdb doesn't exist"
                      << std::endl;
            exit(1);
        }

        std::string line;
        int prevAlpha = -999;
        int prevBeta = -999;
        while (getline(input_stream, line)) {
            std::string line_atom = line.substr(0, 4);
            std::string line_chain = line.substr(21, 1);
            if ((line_atom == "ATOM") && (line_chain == chainName) &&
                (line[16] == ' ' || line[16] == 'A')) {
                std::string atomType = line.substr(13, 2);
                std::string resType = line.substr(17, 3);
                if (atomType == "CA") {
                    proteinSeq.push_back(Codon::AAcode_3_to_1.at(resType));
                    std::string alpha_str = line.substr(22, 4);
                    removeSpaces(alpha_str);
                    int currentAlpha = std::stoi(alpha_str);
                    if (prevAlpha < -100) { prevAlpha = currentAlpha - 1; }
                    if (currentAlpha != (prevAlpha + 1)) {
                        std::cout << "Problem with Alpha " << proteinName << chainName
                                  << ", residue " << currentAlpha << "	" << prevAlpha << std::endl;
                        std::cout << line;
                        std::exit(1);
                    }
                    prevAlpha = currentAlpha;
                    std::array<double, 3> coords{};
                    coords[0] = std::stod(line.substr(30, 8));
                    coords[1] = std::stod(line.substr(38, 8));
                    coords[2] = std::stod(line.substr(46, 8));
                    alphaCoordinateVector.push_back(coords);
                    if (resType == "GLY") {
                        betaCoordinateVector.push_back(coords);
                        prevBeta++;
                    }
                } else if (atomType == "CB") {
                    std::string beta_str = line.substr(22, 4);
                    removeSpaces(beta_str);
                    int currentBeta = std::stoi(beta_str);
                    if (prevBeta < -100) { prevBeta = currentBeta - 1; }
                    if (currentBeta != (prevBeta + 1)) {
                        std::cout << "Problem with Beta " << proteinName << chainName
                                  << ", residue " << currentBeta << "	" << prevBeta << std::endl;
                        std::cout << line;
                        std::exit(1);
                    }
                    prevBeta = currentBeta;
                    std::array<double, 3> coords{};
                    coords[0] = std::stod(line.substr(30, 8));
                    coords[1] = std::stod(line.substr(38, 8));
                    coords[2] = std::stod(line.substr(46, 8));
                    betaCoordinateVector.push_back(coords);
                }
            }
        }
    }

  public:
    std::string proteinSeq;

    explicit Structure(std::string const &pdb_folder,
        std::array<std::string, 2> protein_chain_names, size_t size, double cut_off)
        : proteinName{protein_chain_names[0]},
          chainName{protein_chain_names[1]},
          size{size},
          cut_off{cut_off} {
        for (char c = 0; c < 64; ++c) {
            assert(codonPDB.codon_string(c) == codonLexico.codon_string(c));
        }
        readPDBFile(pdb_folder);
        if (!allOK()) { std::cout << "Yell!!!\n"; }
        makeContactVector();
        std::cout << proteinName << chainName << " with " << contactVector.size()
                  << " contacts and " << alphaCoordinateVector.size() << " residues." << std::endl;
    }

    /**
     * Computes energy of an amino acid sequence in the structure of the protein
     */
    double getEnergy(std::vector<char> const &codonSeq) const {
        double energy = 0.0;
        for (auto const &contact : contactVector) {
            energy += getCodonInteraction(codonSeq[contact[0]], codonSeq[contact[1]]);
        }
        return energy;
    }

    /**
     * Computes energy of an amino acid sequence in the structure of the protein
     */
    double getMutantEnergy(
        std::vector<char> const &codonSeq, size_t site, char codon_from, char codon_to) const {
        assert(codonSeq[site] == codon_from);
        double energy = 0.0;
        for (size_t const &contact : siteContactVector[site]) {
            energy += (getCodonInteraction(codon_to, codonSeq[contact]) -
                       getCodonInteraction(codon_from, codonSeq[contact]));
        }
        return energy;
    }
};

class StabilityLandscape final : public FitnessLandscape {
  public:
    Structure native;
    std::vector<Structure> unfoldedVector{};

    StabilityLandscape(std::string const &pdb_folder, int nbr_sites, double cut_off)
        : native(pdb_folder, nativeProtein, nbr_sites, cut_off) {
        for (auto const &unfolded : unfoldedProteinList) {
            unfoldedVector.emplace_back(Structure(pdb_folder, unfolded, nbr_sites, cut_off));
        }
    }

    u_long nbr_sites() const override { return native.proteinSeq.size(); }

    std::string seq() const { return native.proteinSeq; }
};

double fitness(double deltaG) {
    double factor = exp(-deltaG / TEMPERATURE);
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

    double nativeDeltaG = 0.0;
    double nativeEnergy = 0.0;
    std::vector<double> unfoldedEnergyVector;

  public:
    std::unique_ptr<FitnessState> clone() const override {
        return std::make_unique<StabilityState>(*this);
    };

    explicit StabilityState(StabilityLandscape const &f) : FitnessState(f), f{f} {
        nativeDeltaG = 0.0;
        nativeEnergy = 0.0;
        unfoldedEnergyVector.resize(f.unfoldedVector.size(), 0.0);
    }

    bool operator==(FitnessState const &other) const override {
        auto p = dynamic_cast<StabilityState const *>(&other);
        return log_fitness == other.log_fitness and (nativeDeltaG == p->nativeDeltaG) and
               (nativeEnergy == p->nativeEnergy) and
               (unfoldedEnergyVector == p->unfoldedEnergyVector);
    }

    u_long nbr_sites() const override { return f.nbr_sites(); }

    void update(std::vector<char> const &codon_seq, double const &pop_size) override {
        nativeEnergy = f.native.getEnergy(codon_seq);

        for (size_t i = 0; i < f.unfoldedVector.size(); i++) {
            unfoldedEnergyVector[i] = f.unfoldedVector[i].getEnergy(codon_seq);
        }
        double avgUnfoldedEnergy = avg(unfoldedEnergyVector);
        double sigma2 = variance(unfoldedEnergyVector, avgUnfoldedEnergy);

        nativeDeltaG = nativeEnergy - avgUnfoldedEnergy + (TEMPERATURE * LOG_UNFOLDED_STATES) +
                       (sigma2 / (2.0 * TEMPERATURE));
        log_fitness = logfitness(nativeDeltaG);
    }

    void update(std::vector<char> const &codon_seq, u_long site, char codon_to, bool burn_in,
        double const &pop_size) override {
        nativeEnergy += f.native.getMutantEnergy(codon_seq, site, codon_seq[site], codon_to);

        for (size_t i = 0; i < f.unfoldedVector.size(); i++) {
            unfoldedEnergyVector[i] +=
                f.unfoldedVector[i].getMutantEnergy(codon_seq, site, codon_seq[site], codon_to);
        }
        double avgUnfoldedEnergy = avg(unfoldedEnergyVector);
        double sigma2 = variance(unfoldedEnergyVector, avgUnfoldedEnergy);

        nativeDeltaG = nativeEnergy - avgUnfoldedEnergy + (TEMPERATURE * LOG_UNFOLDED_STATES) +
                       (sigma2 / (2.0 * TEMPERATURE));
        log_fitness = logfitness(nativeDeltaG);
        if (!burn_in) {
            summary_stats["sub-log-fitness"].add(log_fitness);
            summary_stats["sub-ΔG"].add(nativeDeltaG);
            summary_stats["sub-GFold"].add(nativeEnergy);
            summary_stats["sub-GUnfold"].add(avgUnfoldedEnergy);
        }
    }

    double selection_coefficient(std::vector<char> const &codon_seq, u_long site, char codon_to,
        bool burn_in, double const &pop_size) const override {
        if (codonLexico.codon_to_aa[codon_seq[site]] == codonLexico.codon_to_aa[codon_to]) {
            return 0.0;
        }
        double nativeMutantEnergy =
            nativeEnergy + f.native.getMutantEnergy(codon_seq, site, codon_seq[site], codon_to);

        std::vector<double> unfoldedMutantEnergyVector(f.unfoldedVector.size(), 0.0);
        for (size_t i = 0; i < f.unfoldedVector.size(); i++) {
            unfoldedMutantEnergyVector[i] =
                unfoldedEnergyVector[i] +
                f.unfoldedVector[i].getMutantEnergy(codon_seq, site, codon_seq[site], codon_to);
        }
        double avgUnfoldedMutantEnergy = avg(unfoldedMutantEnergyVector);
        double sigma2 = variance(unfoldedMutantEnergyVector, avgUnfoldedMutantEnergy);

        double mutantDeltaG = nativeMutantEnergy - avgUnfoldedMutantEnergy +
                              (TEMPERATURE * LOG_UNFOLDED_STATES) + (sigma2 / (2.0 * TEMPERATURE));

        double s = selection_coefficient_ddG(nativeDeltaG, mutantDeltaG);
        if (!burn_in) {
            summary_stats["mut-s"].add(s);
            summary_stats["mut-ΔG-native"].add(nativeDeltaG);
            summary_stats["mut-ΔG-mutant"].add(mutantDeltaG);
            summary_stats["mut-ΔΔG"].add(mutantDeltaG - nativeDeltaG);
        }
        return s;
    }
};

std::unordered_map<std::string, SummaryStatistic> FitnessState::summary_stats = {
    {"mut-s", SummaryStatistic()}, {"mut-ΔG-native", SummaryStatistic()},
    {"mut-ΔG-mutant", SummaryStatistic()}, {"mut-ΔΔG", SummaryStatistic()},
    {"sub-ΔG", SummaryStatistic()}, {"sub-log-fitness", SummaryStatistic()},
    {"sub-GFold", SummaryStatistic()}, {"sub-GUnfold", SummaryStatistic()}};

class FoldingArgParse {
  protected:
    TCLAP::CmdLine &cmd;

  public:
    explicit FoldingArgParse(TCLAP::CmdLine &cmd) : cmd{cmd} {}
    TCLAP::ValueArg<std::string> pdb_folder{"", "pdb_folder",
        "The folder containing the .pdb files", false, "data/pdbfiles", "string", cmd};
    TCLAP::ValueArg<u_long> nbr_exons{
        "", "nbr_exons", "Number of exons in the protein", false, 12, "u_long", cmd};
    TCLAP::ValueArg<double> cut_off{"", "cut_off",
        "The distance (in angstrom) to determine if 2 sites are in contact", false, 7.0, "double",
        cmd};
    void add_to_trace(Trace &trace) {
        trace.add("#nbr_exons", nbr_exons.getValue());
        trace.add("cut_off", cut_off.getValue());
        trace.add("pdb_folder", pdb_folder.getValue());
    }
};

class StabilityModel : public FitnessModel {
  public:
    StabilityModel(u_long exon_size, FoldingArgParse &args) : FitnessModel() {
        fitness_landscapes.reserve(args.nbr_exons.getValue());
        fitness_states.reserve(args.nbr_exons.getValue());
        StabilityLandscape landscape(
            args.pdb_folder.getValue(), exon_size, args.cut_off.getValue());

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
            std::string tmp =
                dynamic_cast<StabilityLandscape *>(fitness_landscapes.at(exon).get())->seq();
            s += std::string(tmp.begin(), tmp.end());
        }
        return s;
    }
};