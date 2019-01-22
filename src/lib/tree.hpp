#pragma once

#include <cassert>
#include <set>
#include <vector>
#include "nhx-parser.hpp"

// a tree with both a vector of parents and a vector of children
class Tree {
  public:
    using NodeIndex = int;  // guaranteed to be [0... nb_nodes() - 1] U {no_parent}
    const NodeIndex no_parent{-1};
    using BranchIndex = int;  // guaranteed to be [0... nb_nodes() - 2]

    explicit Tree(std::string const& newick_path) {
        std::ifstream tree_stream{newick_path};
        NHXParser parser{tree_stream};
        const AnnotatedTree& input_tree = parser.get_tree();
        root_ = input_tree.root();
        for (NodeIndex node = 0; node < NodeIndex(input_tree.nb_nodes()); node++) {
            parent_.push_back(input_tree.parent(node));
            children_.emplace_back(
                input_tree.children(node).begin(), input_tree.children(node).end());
            name_.push_back(input_tree.tag(node, "name"));
            std::string length = input_tree.tag(node, "length");
            if (length.empty()) {
                length_.push_back(0.0);
            } else {
                length_.push_back(std::stod(length));
            }
            if (is_leaf(node)) { nb_leaves_++; }
        }
        std::vector<double> distances;
        for (NodeIndex node = 0; node < NodeIndex(nb_nodes()); node++) {
            if (is_leaf(node)) {
                double d = distance_to_root(node);
                std::cerr << d << std::endl;
                distances.push_back(d);
            }
        }
        min_distance_ = *std::min_element(distances.begin(), distances.end());
        max_distance_ = *std::max_element(distances.begin(), distances.end());
        if (std::abs(max_distance_ - min_distance_) > 1e-3) {
            std::cerr << "The tree is not ultrametric." << std::endl;
            std::cerr << "The maximum distance from leaf to root is " << max_distance_ << std::endl;
            std::cerr << "The minimum distance from leaf to root is " << min_distance_ << std::endl;
        } else {
            ultrametric_ = true;
            std::cout << "The tree is ultrametric." << std::endl;
        }
    }

    const std::set<NodeIndex>& children(NodeIndex node) const { return children_.at(node); }
    NodeIndex parent(NodeIndex node) const { return parent_.at(node); }
    std::string node_name(NodeIndex node) const { return name_[node]; }
    double node_length(NodeIndex node) const { return length_[node]; }
    double total_length() const {
        double tot = 0.0;
        for (NodeIndex i = 0; i < NodeIndex(nb_nodes()); i++) { tot += node_length(i); }
        return tot;
    }
    NodeIndex root() const { return root_; }
    std::size_t nb_nodes() const { return parent_.size(); }
    bool is_root(NodeIndex i) const { return i == root_; }
    bool is_leaf(NodeIndex i) const { return children_.at(i).size() == 0; }
    int nb_branches() const { return nb_nodes() - 1; }
    int nb_leaves() const { return nb_leaves_; }
    bool is_ultrametric() const { return ultrametric_; }
    double max_distance_to_root() const { return max_distance_; }
    double min_distance_to_root() const { return min_distance_; }

    //! find the longest path from this node to the farthest leaf.
    double distance_to_root(Tree::NodeIndex node) const {
        double d = node_length(node);
        while (!(is_root(node))) {
            node = parent(node);
            d += node_length(node);
        }
        return d;
    };

    void set_root_age(double root_age) {
        if (!ultrametric_) {
            std::cerr << "The longest distance between root to leaf is used to set the root age."
                      << std::endl;
        }
        root_age /= max_distance_;
        for (auto& length : length_) { length *= root_age; }
    };

    BranchIndex branch_index(NodeIndex i) const { return i - 1; }
    NodeIndex node_index(BranchIndex i) const { return i + 1; }

  private:
    std::vector<NodeIndex> parent_;
    std::vector<std::set<NodeIndex>> children_;
    std::vector<std::string> name_;
    std::vector<double> length_;
    NodeIndex root_{0};
    int nb_leaves_{0};
    bool ultrametric_{false};
    double max_distance_ = 0.0;
    double min_distance_ = 0.0;
};