#include <assert.h>
#include <vector>
#include <array>
#include <map>
#include <algorithm>

// Definitions:
// Nucleotide: a char from 0 to 3 (included) encoding one of the nucleotide (ATGC).
// triplet: an array of 3 nucleotides (first, second and third position).
// codon: a char from 0 to 63 (included) encoding for a triplet.
// neighboring codons: two codons differing by only 1 mutation.
// neighbor: a 3-tuple containing - in 1st position the codon after mutation,
//                                - in 2nd position the nucleotide before mutation,
//                                - in 3rd position the nucleotide after mutation.

namespace Codon {
    double epsilon = std::numeric_limits<double>::epsilon();

    // String of all nucleotides.
    std::string const nucleotides{"ACGT"};


    // Function to map a triplet of 3 nucleotides (1st, 2nd and 3rd positions) to the corresponding codon.
    // n_1 : the nucleotide in 1st position
    // n_2 : the nucleotide in 2nd position
    // n_3 : the nucleotide in 3rd position
    char triplet_to_codon(char n_1, char n_2, char n_3) {
        assert(0 <= n_1 and n_1 <= 3);
        assert(0 <= n_2 and n_2 <= 3);
        assert(0 <= n_3 and n_3 <= 3);
        int codon{n_1 * 16 + n_2 * 4 + n_3};
        assert(0 <= codon and codon <= 63);
        return char(codon);
    }

    // Function to build an array mapping each of the 64 codons to their respective triplet.
    std::array<std::array<char, 3>, 64> build_codon_to_triplet_array() {
        std::array<std::array<char, 3>, 64> codon_to_triplet_vec = {0};

        for (char n_1{0}; n_1 < 4; n_1++) {
            // n_1 is the first position of the triplet.
            for (char n_2{0}; n_2 < 4; n_2++) {
                // n_2 is the second position of the triplet.
                for (char n_3{0}; n_3 < 4; n_3++) {
                    // n_3 is the third position of the triplet.
                    std::array<char, 3> triplet = {n_1, n_2, n_3};
                    codon_to_triplet_vec[triplet_to_codon(n_1, n_2, n_3)] = triplet;
                }
            }
        }
        return codon_to_triplet_vec;
    }

    // Function to build an array which maps each of the 64 codons to the set of their 9 respective neighbors.
    // The set of neighbors contains the 9 codons which differ by only 1 nucleotide (3 positions * 3 possible mutations).
    // codon_to_triplet: array mapping each of the 64 codons to their respective triplet (see build_codon_to_triplet_array)
    std::array<std::array<std::tuple<char, char, char>, 9>, 64>
    build_codon_to_neighbors_array(std::array<std::array<char, 3>, 64> const &codon_to_triplet) {
        // The array which maps the 64 codons to the set of their 9 respective neighbors.
        std::array<std::array<std::tuple<char, char, char>, 9>, 64> codon_to_neighbors;

        // For each possible codon.
        for (char codon{0}; codon < 64; codon++) {

            // The triplet corresponding to this codon.
            std::array<char, 3> triplet_from = codon_to_triplet[codon];

            // For each position in the triplet (3 positions).
            for (char position{0}; position < 3; position++) {
                // The original nucleotide which is going to be mutated.
                char n_from = triplet_from[position];

                // For each possible mutation (3 possible mutations).
                for (char mutation{0}; mutation < 3; mutation++) {
                    // The mutated nucleotide.
                    char n_to = n_from;
                    n_to++;
                    n_to += mutation;
                    n_to %= 4;

                    // The mutated triplet.
                    std::array<char, 3> triplet_to = triplet_from;
                    triplet_to[position] = n_to;

                    // The mutated codon.
                    char codon_to = triplet_to_codon(triplet_to[0], triplet_to[1], triplet_to[2]);

                    // Assign the neighbor to the array of neighbors.
                    codon_to_neighbors[codon][3 * position + mutation] = std::make_tuple(codon_to, n_from, n_to);
                }
            }
        }
        return codon_to_neighbors;
    }

    // Array mapping each of the 64 codons to their respective triplet.
    auto codon_to_triplet_array = build_codon_to_triplet_array();

    // The set of neighbors contains the 9 codons which differ by only 1 nucleotide (3 positions * 3 possible mutations).
    // Array which maps each of the 64 codons to the set of their 9 respective neighbors.
    auto const codon_to_neighbors_array = build_codon_to_neighbors_array(codon_to_triplet_array);

    // String containing the 21 amino-acids (20 + stop) for translations from index to amino-acid character
    std::string const amino_acids{"ACDEFGHIKLMNPQRSTVWYX"};

    // Map from codons to amino-acids
    std::map<std::string, char> const triplet_str_to_aa_char{{"GCA", 'A'},
                                                             {"GAA", 'E'},
                                                             {"ACT", 'T'},
                                                             {"CAT", 'H'},
                                                             {"ACG", 'T'},
                                                             {"GGT", 'G'},
                                                             {"GCG", 'A'},
                                                             {"GAG", 'E'},
                                                             {"CGC", 'R'},
                                                             {"TGA", 'X'},
                                                             {"CGG", 'R'},
                                                             {"GCC", 'A'},
                                                             {"TGC", 'C'},
                                                             {"GAC", 'D'},
                                                             {"CAA", 'Q'},
                                                             {"CGT", 'R'},
                                                             {"GAT", 'D'},
                                                             {"TCA", 'S'},
                                                             {"CAC", 'H'},
                                                             {"ATC", 'I'},
                                                             {"CGA", 'R'},
                                                             {"ATA", 'I'},
                                                             {"GCT", 'A'},
                                                             {"CAG", 'Q'},
                                                             {"TGG", 'W'},
                                                             {"GGC", 'G'},
                                                             {"TTC", 'F'},
                                                             {"CCA", 'P'},
                                                             {"ACC", 'T'},
                                                             {"TAC", 'Y'},
                                                             {"GTG", 'V'},
                                                             {"AAC", 'N'},
                                                             {"AAG", 'K'},
                                                             {"CCT", 'P'},
                                                             {"TCC", 'S'},
                                                             {"CCC", 'P'},
                                                             {"CTC", 'L'},
                                                             {"GTT", 'V'},
                                                             {"AGC", 'S'},
                                                             {"ATT", 'I'},
                                                             {"ACA", 'T'},
                                                             {"TTG", 'L'},
                                                             {"GTC", 'V'},
                                                             {"AGT", 'S'},
                                                             {"CTG", 'L'},
                                                             {"TCG", 'S'},
                                                             {"TAT", 'Y'},
                                                             {"TTT", 'F'},
                                                             {"AAT", 'N'},
                                                             {"CCG", 'P'},
                                                             {"TTA", 'L'},
                                                             {"TGT", 'C'},
                                                             {"GGA", 'G'},
                                                             {"CTA", 'L'},
                                                             {"AAA", 'K'},
                                                             {"GGG", 'G'},
                                                             {"ATG", 'M'},
                                                             {"GTA", 'V'},
                                                             {"TCT", 'S'},
                                                             {"AGA", 'R'},
                                                             {"TAA", 'X'},
                                                             {"TAG", 'X'},
                                                             {"AGG", 'R'},
                                                             {"CTT", 'L'}};

    // Function to map a particular amino-acid to its index in the string containing all the amino-acids
    // Equivalent to the inverse function of translating index to amino-acid character
    // aa_char: character of an amino-acid
    // amino_acids: string containing all the amino-acids
    char aa_char_to_aa(char const aa_char, std::string const &amino_acids) {
        return char(amino_acids.find(aa_char));
    }

    // Function to build an array mapping each of the 64 codons to their respective amino-acid.
    // triplet_str_to_aa_char: map from codons to amino-acids
    // amino_acids: string containing all the amino-acids
    // codon_to_triplet: array mapping each of the 64 codons to their respective triplet (see build_codon_to_triplet_array)
    std::array<char, 64> build_codon_to_aa_array(std::map<std::string, char> const &triplet_str_to_aa_char,
                                                 std::string const &amino_acids,
                                                 std::array<std::array<char, 3>, 64> const &codon_to_triplet) {
        // Array mapping each of the 64 codons to their respective amino-acid.
        std::array<char, 64> codon_to_aa = {0};

        // For each codon
        for (char codon{0}; codon < 64; codon++) {

            // Triplet corresponding to the codon
            std::array<char, 3> triplet = codon_to_triplet_array[codon];

            // String containing the translated triplet (e.g "ATG")
            std::string triplet_str(3, ' ');
            for (char position{0}; position < 3; position++) {
                triplet_str[position] = nucleotides[triplet[position]];
            }

            // Amino-acid corresponding to the translated string
            codon_to_aa[codon] = aa_char_to_aa(triplet_str_to_aa_char.at(triplet_str), amino_acids);
        }
        return codon_to_aa;
    }

    // Array mapping each of the 64 codons to their respective amino-acid.
    std::array<char, 64> const codon_to_aa_array = build_codon_to_aa_array(triplet_str_to_aa_char, amino_acids,
                                                                           codon_to_triplet_array);

    // Random generator engine with seed 0.
    double seed{0};
    std::default_random_engine re(seed);
}