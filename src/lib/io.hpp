#pragma once

#include <algorithm>
#include <deque>
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

void init_alignments(
    std::string const &output_path, u_long nb_leaves, u_long nb_sites, bool fasta = false) {
    if (fasta) {
        // .fasta format
        std::ofstream fasta_file;
        fasta_file.open(output_path + ".fasta");
        fasta_file.close();
    } else {
        // .ali format
        std::ofstream ali_file;
        ali_file.open(output_path + ".ali");
        ali_file << nb_leaves << " " << nb_sites << std::endl;
        ali_file.close();
    }
}

// Write sequence in ali and fasta files
void write_sequence(std::string const &output_filename, std::string const &name,
    std::string const &dna_str, bool fasta = false) {
    if (fasta) {
        // .fasta format
        std::ofstream fasta_file;
        fasta_file.open(output_filename + ".fasta", std::ios_base::app);
        fasta_file << ">" << name << std::endl << dna_str << std::endl;
        fasta_file.close();
    } else {
        // .ali format
        std::ofstream ali_file;
        ali_file.open(output_filename + ".ali", std::ios_base::app);
        ali_file << name << " " << dna_str << std::endl;
        ali_file.close();
    }
}

std::string join(std::vector<std::string> const &v, char sep) {
    if (v.empty()) { return ""; }
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
    std::unordered_map<std::string, std::size_t> header_to_index;
    std::unordered_map<std::string, std::size_t> header_to_count;
    std::size_t nb_row{0};
    std::size_t nb_col{0};
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

class Distribution {
  private:
    double grid_step = 0.1;
    double grid_min = -10;
    double grid_max = 10;
    u_long bin_zero = 0;
    std::deque<u_long> bins{};
    std::deque<double> bounds{};


    static struct {
        bool operator()(const double &left, const double &right) { return left < right; }
    } LowerThan;

    void UpdateComputed(double s) {
        while (bounds.front() > s) {
            double new_s = bounds.front() - grid_step;
            bounds.emplace_front(new_s);
            bins.emplace_front(0);
        }
        while (bounds.back() < s) {
            double new_s = bounds.back() + grid_step;
            bounds.emplace_back(new_s);
            bins.emplace_back(0);
        }
    }


  public:
    Distribution() {
        double s = -1;
        while (s <= 1) {
            bounds.emplace_back(s);
            bins.emplace_back(0);
            s += grid_step;
        }
    }

    void add(double s) {
        if (s == 0.0) {
            bin_zero++;
            return;
        }
        if (s < grid_min) { s = grid_min; }
        if (s > grid_max) { s = grid_max; }
        if (bounds.front() > s or bounds.back() < s) { UpdateComputed(s); }

        // The index for the closest (lower and upper) selection coefficient for which
        // pre-computation is available
        auto it_up = std::upper_bound(bounds.begin(), bounds.end(), s, LowerThan);
        auto it_low = prev(it_up);
        auto id = it_low - bounds.begin();
        assert(bounds.at(it_low - bounds.begin()) <= s);
        assert(bounds.at(it_up - bounds.begin()) >= s);
        bins.at(id) += 1;
    }

    void write_tsv(std::string const &output_filename) {
        std::cout << "Writing " << output_filename << " with " << bounds.size() << " columns"
                  << std::endl;
        std::cout << "Min =" << bounds.front() << "; " << bins.front() << " counts" << std::endl;
        std::cout << "Max =" << bounds.back() << "; " << bins.back() << " counts" << std::endl;
        std::ofstream tsv_file;
        tsv_file.open(output_filename + "distrib.tsv");
        while (!bounds.empty()) {
            tsv_file << std::to_string(bounds.front()) << "\t";
            bounds.pop_front();
        }
        tsv_file << "BinNull" << std::endl;
        while (!bins.empty()) {
            tsv_file << std::to_string(bins.front()) << "\t";
            bins.pop_front();
        }
        tsv_file << bin_zero << std::endl;
        tsv_file.close();
        std::cout << "Written !" << std::endl;
    }
};

class DistributionMap : public std::unordered_map<std::string, Distribution> {
  public:
    DistributionMap() = default;

    void write_tsv(std::string const &output_filename) {
        for (auto &d : *this) { d.second.write_tsv(output_filename + d.first); }
    }
};

DistributionMap distribution_map{};