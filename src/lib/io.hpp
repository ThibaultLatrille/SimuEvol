#pragma once

#include <fstream>
#include <iostream>
#include <numeric>
#include <sstream>
#include <unordered_map>
#include <vector>


std::string char_to_str(char const &_char) {
    std::string _str(1, _char);
    return _str;
}

void init_alignments(std::string const &output_path, u_long nb_leaves, u_long nb_sites) {
    // .ali format
    std::ofstream ali_file;
    ali_file.open(output_path + ".ali");
    ali_file << nb_leaves << " " << nb_sites << std::endl;
    ali_file.close();

    // .ali format
    std::ofstream fasta_file;
    fasta_file.open(output_path + ".fasta");
    fasta_file.close();
}

// Write sequence in ali and fasta files
void write_sequence(
    std::string const &output_filename, std::string const &name, std::string const &dna_str) {
    // .ali format
    std::ofstream ali_file;
    ali_file.open(output_filename + ".ali", std::ios_base::app);
    ali_file << name << " " << dna_str << std::endl;
    ali_file.close();

    // .fasta format
    std::ofstream fasta_file;
    fasta_file.open(output_filename + ".fasta", std::ios_base::app);
    fasta_file << ">" << name << std::endl << dna_str << std::endl;
    fasta_file.close();
}

std::string join(std::vector<std::string> const &v, char sep) {
    return std::accumulate(v.begin() + 1, v.end(), v[0],
        [sep](std::string const &acc, std::string const &b) { return acc + sep + b; });
}

std::string d_to_string(double val) {
    std::ostringstream so;
    so << std::scientific << val;
    return so.str();
}

class Trace {
  private:
    std::unordered_map<std::string, size_t> header_to_index;
    std::unordered_map<std::string, size_t> header_to_count;
    size_t nb_row{0};
    size_t nb_col{0};
    std::vector<std::string> header;
    std::vector<std::vector<std::string>> data;

  public:
    Trace() : header_to_index{}, header_to_count{}, data{} {}

    void add(std::string const &key, std::string const &val) {
        if (header_to_index.count(key) == 0) {
            header.push_back(key);
            header_to_count[key] = 0;
            header_to_index[key] = nb_col;
            nb_col++;
            for (auto &row : data) { row.resize(nb_col); }
        } else {
            header_to_count[key]++;
        }
        if (header_to_count[key] >= nb_row) {
            nb_row++;
            data.resize(nb_row);
            data.back().resize(nb_col);
        }
        data.at(header_to_count[key]).at(header_to_index[key]) = val;
    }

    void add(std::string const &key, bool val) { add(key, std::to_string(val)); }
    void add(std::string const &key, char val) { add(key, char_to_str(val)); }
    void add(std::string const &key, int val) { add(key, std::to_string(val)); }
    void add(std::string const &key, u_long val) { add(key, std::to_string(val)); }
    void add(std::string const &key, unsigned long long val) { add(key, std::to_string(val)); }
    void add(std::string const &key, long val) { add(key, std::to_string(val)); }
    void add(std::string const &key, unsigned val) { add(key, std::to_string(val)); }
    void add(std::string const &key, double val) { add(key, d_to_string(val)); }

    void write_tsv(std::string const &output_filename) {
        std::ofstream tsv_file;
        tsv_file.open(output_filename + ".tsv");
        tsv_file << join(header, '\t') << std::endl;
        for (auto const &row : data) {
            assert(row.size() == header_to_index.size());
            assert(row.size() == header_to_count.size());
            tsv_file << join(row, '\t') << std::endl;
        }
        tsv_file.close();
    }
};
