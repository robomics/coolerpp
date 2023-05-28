// Copyright (C) 2023 Roberto Rossini <roberros@uio.no>
//
// SPDX-License-Identifier: MIT

#include <fmt/format.h>

#include <algorithm>
#include <cassert>
#include <cstdint>
#include <queue>
#include <stdexcept>
#include <string_view>
#include <type_traits>
#include <utility>
#include <vector>

#include "coolerpp/coolerpp.hpp"
#include "coolerpp_tools/config.hpp"
#include "coolerpp_tools/tools.hpp"

namespace coolerpp::tools {

[[nodiscard]] static std::vector<File> open_coolers(const std::vector<std::string>& uris) {
  std::vector<File> clrs(uris.size());
  std::uint32_t bin_size = 0;
  std::transform(uris.begin(), uris.end(), clrs.begin(), [&](const auto& uri) {
    auto clr = File::open_read_only(uri);
    if (bin_size == 0) {
      bin_size = clr.bin_size();
    } else if (bin_size != clr.bin_size()) {
      throw std::runtime_error("TODO");
    }
    return clr;
  });
  return clrs;
}

[[nodiscard]] static ChromosomeSet get_chromosomes(const std::vector<File>& clrs) {
  assert(clrs.size() > 1);
  auto chroms = clrs.front().chromosomes();

  std::for_each(clrs.begin() + 1, clrs.end(), [&](const File& clr) {
    if (chroms != clr.chromosomes()) {
      throw std::runtime_error(
          fmt::format(FMT_STRING("input coolers have different axes: found difference in the "
                                 "chromosome table of Coolers at the following URIs:\n"
                                 "- \"{}\""
                                 "- \"{}\""),
                      clrs.front().uri(), clr.uri()));
    }
  });

  return chroms;
}

template <typename N>
class PixelMerger {
  struct Node {
    Pixel<N> pixel{};
    std::size_t i{};

    bool operator<(const Node& other) const noexcept {
      assert(!!this->pixel);
      assert(!!other.pixel);
      return this->pixel < other.pixel;
    }
    bool operator==(const Node& other) const noexcept {
      return this->pixel.coords == other.pixel.coords;
    }
  };

  std::vector<Pixel<N>> _buffer{};
  std::priority_queue<Node, std::vector<Node>, std::less<>> _pqueue{};
  using PixelIt = decltype(std::declval<File>().begin<N>());

  std::vector<PixelIt> _heads{};
  std::vector<PixelIt> _tails{};

 public:
  PixelMerger() = delete;
  explicit PixelMerger(const std::vector<File>& input_coolers) {
    for (const auto& clr : input_coolers) {
      auto first = clr.begin<N>();
      auto last = clr.end<N>();
      if (first != last) {
        auto pixel = *first++;
        _heads.emplace_back(std::move(first));
        _tails.emplace_back(std::move(last));
        _pqueue.emplace(Node{std::move(pixel), _pqueue.size()});
      }
    }
  }

  void merge(File& clr, std::size_t capacity = 2'000'000) {
    this->_buffer.clear();
    this->_buffer.reserve(std::max(capacity, this->_buffer.capacity()));

    std::size_t pixels_processed{};
    while (true) {
      auto pixel = this->next();
      if (!pixel) {
        break;
      }
      this->_buffer.emplace_back(std::move(pixel));
      if (this->_buffer.size() == capacity) {
        clr.append_pixels(this->_buffer.begin(), this->_buffer.end());
        pixels_processed += this->_buffer.size();
        if (pixels_processed % std::max(capacity, std::size_t(1'000'000)) == 0) {
          fmt::print(stderr, FMT_STRING("Procesed {}M pixels...\n"), pixels_processed / 1'000'000);
        }
        this->_buffer.clear();
      }
    }

    if (!this->_buffer.empty()) {
      clr.append_pixels(this->_buffer.begin(), this->_buffer.end());
    }
  }

 private:
  void replace_top_node(std::size_t i) {
    assert(this->_pqueue.top().i == i);
    this->_pqueue.pop();
    if (auto& it = this->_heads[i]; it != this->_tails[i]) {
      // TODO: why is this a lot faster than Node{*(it++)}?
      this->_pqueue.emplace(Node{*it, i});
      ++it;
    }
  }

  [[nodiscard]] Pixel<N> next() {
    if (this->_pqueue.empty()) {
      return {};
    }

    auto current_node = this->_pqueue.top();
    this->replace_top_node(current_node.i);

    for (auto next_node = this->_pqueue.top(); !this->_pqueue.empty() && next_node == current_node;
         next_node = this->_pqueue.top()) {
      current_node.pixel.count += next_node.pixel.count;
      this->replace_top_node(next_node.i);
    }
    return current_node.pixel;
  }
};

template <typename N>
static void merge_coolers(const std::vector<File>& clrs, File&& out) {
  PixelMerger<N>{clrs}.merge(out);
}

void merge_subcmd(const MergeConfig& c) {
  const auto clrs = open_coolers(c.input_uris);
  assert(clrs.size() > 1);

  const auto bin_size = clrs.front().bin_size();
  const auto chroms = get_chromosomes(clrs);

  if (c.floating_point) {
    merge_coolers<double>(clrs,
                          File::create_new_cooler<double>(c.output_uri, chroms, bin_size, c.force));
  } else {
    merge_coolers<std::int64_t>(
        clrs, File::create_new_cooler<std::int64_t>(c.output_uri, chroms, bin_size, c.force));
  }
}

}  // namespace coolerpp::tools
