// Copyright (C) 2022 Roberto Rossini <roberros@uio.no>
//
// SPDX-License-Identifier: MIT

#pragma once

#include <cstdint>
#include <queue>
#include <string_view>
#include <utility>
#include <vector>

#include "coolerpp/coolerpp.hpp"

namespace coolerpp::utils {

enum class MergeStrategy { IN_MEMORY, PQUEUE };

/// Iterable of coolerpp::File or strings
template <typename Str>
void merge(Str first_file, Str last_file, std::string_view dest_uri,
           bool overwrite_if_exists = false, MergeStrategy strat = MergeStrategy::PQUEUE);

[[nodiscard]] bool equal(std::string_view uri1, std::string_view uri2,
                         bool ignore_attributes = true);
[[nodiscard]] bool equal(const File& clr1, const File& clr2, bool ignore_attributes = true);

namespace internal {

/// This class is basically a wrapper around a priority queue of objects of type Node
/// Node consist of a pixel and an index. The index represent from which iterator (i.e. file) the
/// pixel was read. This allows us to know from which iterator we should read the next pixel (i.e.
/// the same iterator from which the top pixel originated)
template <typename N>
class PixelMerger {
  struct Node {
    Pixel<N> pixel{};
    std::size_t i{};

    bool operator<(const Node& other) const noexcept;
    bool operator>(const Node& other) const noexcept;
    bool operator==(const Node& other) const noexcept;
    bool operator!=(const Node& other) const noexcept;
  };

  std::vector<Pixel<N>> _buffer{};
  std::priority_queue<Node, std::vector<Node>, std::greater<>> _pqueue{};
  using PixelIt = decltype(std::declval<File>().begin<N>());

  std::vector<PixelIt> _heads{};
  std::vector<PixelIt> _tails{};

 public:
  PixelMerger() = delete;
  explicit PixelMerger(const std::vector<File>& input_coolers);
  template <typename FileIt>
  PixelMerger(FileIt first_file, FileIt last_file);
  void merge(File& clr, MergeStrategy strat, std::size_t queue_capacity = 2'000'000,
             bool quiet = true);

 private:
  void replace_top_node(std::size_t i);
  [[nodiscard]] Pixel<N> next();
  void merge_pqueue(File& clr, std::size_t queue_capacity = 2'000'000, bool quiet = true);
  void merge_in_memory(File& clr, bool quiet = true);
};
}  // namespace internal
}  // namespace coolerpp::utils

#include "../../utils_equal_impl.hpp"
#include "../../utils_merge_impl.hpp"
