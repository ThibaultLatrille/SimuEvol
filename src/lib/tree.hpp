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

    explicit Tree(const AnnotatedTree& input_tree) {
        root_ = input_tree.root();
        for (std::size_t i = 0; i < input_tree.nb_nodes(); i++) {
            parent_.push_back(input_tree.parent(i));
            children_.emplace_back(input_tree.children(i).begin(), input_tree.children(i).end());
            name_.push_back(input_tree.tag(i, "name"));
            std::string length = input_tree.tag(i, "length");
            if (length == ""){
                length_.push_back(0.0);
            } else {
                length_.push_back(std::stod(length));
            }
        }
    }

    const std::set<NodeIndex>& children(NodeIndex node) const { return children_.at(node); }
    NodeIndex parent(NodeIndex node) const { return parent_.at(node); }
    std::string node_name(NodeIndex node) const { return name_[node]; }
    double node_length(NodeIndex node) const { return length_[node]; }
    double total_length() const {
        double tot = 0.0;
        for (NodeIndex i = 0; i < nb_nodes(); i++) { tot += node_length(i); }
        return tot;
    }
    NodeIndex root() const { return root_; }
    std::size_t nb_nodes() const { return parent_.size(); }
    bool is_root(NodeIndex i) const { return i == root_; }
    bool is_leaf(NodeIndex i) const { return children_.at(i).size() == 0; }
    int nb_branches() const { return nb_nodes() - 1; }
    int nb_leaves() const {
        int leaves = 0;
        for (NodeIndex i = 0; i < nb_nodes(); i++) {
            if (is_leaf(i)) { leaves++; }
        }
        return leaves;
    }
    BranchIndex branch_index(NodeIndex i) const { return i - 1; }
    NodeIndex node_index(BranchIndex i) const { return i + 1; }

  private:
    std::vector<NodeIndex> parent_;
    std::vector<std::set<NodeIndex>> children_;
    std::vector<std::string> name_;
    std::vector<double> length_;
    NodeIndex root_{0};
};

std::unique_ptr<const Tree> make_from_parser(TreeParser& parser) {
    return std::unique_ptr<Tree>(new Tree(parser.get_tree()));
}
