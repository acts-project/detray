/** Detray library, part of the ACTS project (R&D line)
 *
 * (c) 2021 CERN for the benefit of the ACTS project
 *
 * Mozilla Public License Version 2.0
 */
#pragma once

#include <functional>
#include <iterator>
#include <type_traits>
#include <utility>

#include "detray/utils/enumerate.hpp"
#include "detray/utils/indexing.hpp"

namespace detray {

/** Default node in the hash tree.*/
template <typename hash_type>
struct hashed_node {

    hashed_node(hash_type hash)
        : _key(hash), _parent(-1), _left_child(-1), _right_child(-1) {}

    hash_type _key;
    dindex _parent, _left_child, _right_child;
    hash_type &key() { return _key; }

    void set_parent(dindex pi) { _parent = pi; }

    void set_children(dindex lc, dindex rc) {
        _left_child = lc;
        _right_child = rc;
    }
};

/** Default hash function makes use of default node in the hash tree.*/
template <typename value_type>
struct default_hash {
    using hash_type = std::size_t;
    using data_hasher = std::hash<value_type>;
    using node_hasher = std::hash<hash_type>;

    decltype(auto) operator()(const std::pair<value_type, value_type> &v) {
        return std::make_pair(data_hasher{}(v.first), data_hasher{}(v.second));
    }

    decltype(auto) operator()(const hash_type &left, const hash_type &right) {
        return node_hasher{}(left) + node_hasher{}(right);
    }
};

/**
 * @brief Builds a hash tree from the data of an input collection.
 *
 * This class provides a graph algorithm that walk along the volumes of a given
 * geomtery and uses the portals to check reachability between the volumes.
 *
 * @tparam hash_function the type hashing that can be called on the collection
 *                       data.
 * @tparam input collection the type of data collection to be hased
 * @tparam node_type how carries the hashes and links
 */
template <typename input_collection,
          typename data_type = typename input_collection::value_type,
          typename hash_function = default_hash<data_type>,
          // std::enable_if_t<std::is_invocable_v<hash_function, data_type>,
          // bool> = true, typename node_type =
          // hashed_node<std::result_of_t<hash_function(data_type)>>,
          typename node_type = hashed_node<std::size_t>,
          template <typename> class vector_type = dvector>
class hash_tree {

    public:
    /** No empty tree */
    hash_tree() = delete;

    /** Build from existing nodes and edges, which are provide by the geometry.
     *
     * @param volumes geometry volumes that become the graph nodes
     * @param portals geometry portals link volumes and become edges
     */
    hash_tree(const input_collection &data, const hash_function & /*hf*/ = {})
        : _hash(hash_function{}) {
        build(data);
    }

    /** Default destructor */
    ~hash_tree() = default;

    /** @return the root hash of the tree, which is always the last node in the
     *         node storage by way of construction. It is the fingeprint of the
     *         input data.
     */
    auto root() { return _tree.back().key(); }

    private:
    /** Go through the the sorted input data and recursively build the tree.
     */
    void build(const input_collection &input_data) {
        // Build leaves
        for (const auto &data : input_data) {
            // Construct leaves from input data type
            auto [left, right] = _hash(data);

            _tree.emplace_back(left);
            _tree.emplace_back(right);
        }
        // Size of the tree is already known (all iterators stay valid)
        _tree.reserve(2 * _tree.size() - 1);
        // Build next level
        build(_tree.begin(), _tree.size());
    }

    /** Build the hash tree recursively.
     *
     * @param first_child the beginning of the nodes for which to construct
     *                    the parents in this iteration.
     */
    template <typename iterator_type>
    void build(iterator_type &&first_child, std::size_t n_at_level) {
        auto last_child = first_child + n_at_level;
        // base case
        if (first_child == last_child) {
            return;
        }
        // Run over previous tree level to build the next level
        for (auto current_child = first_child; current_child != last_child;
             ++current_child) {
            auto parent_hash =
                _hash(current_child->key(), ++current_child->key());
            node_type parent = _tree.emplace_back(parent_hash);

            // Parent node index is at the back of the tree
            (current_child - 1)->set_parent(_tree.size() - 1);
            current_child->set_parent(_tree.size() - 1);

            // Set the indices as distances in the contiguous container
            dindex right_child_idx =
                std::distance(current_child, _tree.begin());
            parent.set_children(right_child_idx - 1, right_child_idx);
        }
        // begin next time where we ended this time
        build(last_child, static_cast<std::size_t>(0.5 * n_at_level));
    }

    /** How to encode tha node data */
    hash_function _hash;

    /** Tree nodes */
    vector_type<node_type> _tree = {};
};

}  // namespace detray