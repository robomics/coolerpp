// Copyright (C) 2023 Roberto Rossini <roberros@uio.no>
//
// SPDX-License-Identifier: MIT

#pragma once

#include <algorithm>
#include <queue>
#include <set>
#include <string_view>
#include <vector>

#include "coolerpp/coolerpp.hpp"

namespace coolerpp::utils {

namespace internal {
template <typename N>
inline bool PixelMerger<N>::Node::operator<(const Node& other) const noexcept {
  assert(!!this->pixel);
  assert(!!other.pixel);
  return this->pixel.coords < other.pixel.coords;
}

template <typename N>
inline bool PixelMerger<N>::Node::operator>(const Node& other) const noexcept {
  assert(!!this->pixel);
  assert(!!other.pixel);
  return this->pixel.coords > other.pixel.coords;
}

template <typename N>
inline bool PixelMerger<N>::Node::operator==(const Node& other) const noexcept {
  return this->pixel.coords == other.pixel.coords;
}

template <typename N>
inline bool PixelMerger<N>::Node::operator!=(const Node& other) const noexcept {
  return !(*this == other);
}

template <typename N>
inline PixelMerger<N>::PixelMerger(const std::vector<File>& input_coolers)
    : PixelMerger<N>(input_coolers.begin(), input_coolers.end()) {}

template <typename N>
template <typename FileIt>
inline PixelMerger<N>::PixelMerger(FileIt first_file, FileIt last_file) {
  std::for_each(first_file, last_file, [&](const auto& clr) {
    auto first = clr.template begin<N>();
    auto last = clr.template end<N>();
    if (first != last) {
      auto pixel = *first++;
      _heads.emplace_back(std::move(first));
      _tails.emplace_back(std::move(last));
      _pqueue.emplace(Node{std::move(pixel), _pqueue.size()});
    }
  });
}

template <typename N>
inline void PixelMerger<N>::merge(File& clr, MergeStrategy strat, std::size_t queue_capacity,
                                  bool quiet) {
  if (strat == MergeStrategy::IN_MEMORY) {
    return this->merge_in_memory(clr, quiet);
  }
  return this->merge_pqueue(clr, queue_capacity, quiet);
}

template <typename N>
inline void PixelMerger<N>::merge_pqueue(File& clr, std::size_t queue_capacity, bool quiet) {
  this->_buffer.clear();
  this->_buffer.reserve(std::max(queue_capacity, this->_buffer.capacity()));

  std::size_t pixels_processed{};
  while (true) {
    auto pixel = this->next();
    if (!pixel) {
      break;
    }
    this->_buffer.emplace_back(std::move(pixel));
    if (this->_buffer.size() == queue_capacity) {
      clr.append_pixels(this->_buffer.begin(), this->_buffer.end());
      pixels_processed += this->_buffer.size();
      if (!quiet && pixels_processed % std::max(queue_capacity, std::size_t(1'000'000)) == 0) {
        fmt::print(stderr, FMT_STRING("Procesed {}M pixels...\n"), pixels_processed / 1'000'000);
      }
      this->_buffer.clear();
    }
  }

  if (!this->_buffer.empty()) {
    clr.append_pixels(this->_buffer.begin(), this->_buffer.end());
  }
}

template <typename N>
inline void PixelMerger<N>::merge_in_memory(File& clr, bool quiet) {
  auto compare_pixel_coords = [](const Pixel<N>& p1, const Pixel<N>& p2) {
    return p1.coords < p2.coords;
  };

  using ComparePixelCoordsT = decltype(compare_pixel_coords);
  std::set<Pixel<N>, ComparePixelCoordsT> pixel_tank{compare_pixel_coords};

  std::size_t pixels_processed = 0;
  auto emplace_or_merge_pixel = [&](Pixel<N> p) {
    auto match = pixel_tank.find(p);
    if (match != pixel_tank.end()) {
      p.count += match->count;
      pixel_tank.erase(match);
    }
    pixel_tank.emplace(std::move(p));

    if (!quiet && ++pixels_processed % 1'000'000 == 0) {
      fmt::print(stderr, FMT_STRING("Procesed {}M pixels...\n"), pixels_processed / 1'000'000);
    }
  };

  while (!this->_pqueue.empty()) {
    emplace_or_merge_pixel(this->_pqueue.top().pixel);
    this->_pqueue.pop();
  }

  for (std::size_t i = 0; i < this->_heads.size(); ++i) {
    std::for_each(this->_heads[i], this->_tails[i], emplace_or_merge_pixel);
  }
  clr.append_pixels(std::make_move_iterator(pixel_tank.begin()),
                    std::make_move_iterator(pixel_tank.end()));
  clr.flush();
}

template <typename N>
inline void PixelMerger<N>::replace_top_node(std::size_t i) {
  assert(this->_pqueue.top().i == i);
  this->_pqueue.pop();
  if (auto& it = this->_heads[i]; it != this->_tails[i]) {
    this->_pqueue.emplace(Node{*it++, i});
  }
}

template <typename N>
inline Pixel<N> PixelMerger<N>::next() {
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

[[nodiscard]] inline std::uint32_t get_bin_size_checked(const std::vector<File>& coolers) {
  assert(coolers.size() > 1);
  const auto& clr1 = coolers.front();

  for (std::size_t i = 1; i < coolers.size(); ++i) {
    const auto& clr2 = coolers[i];
    if (clr1.bin_size() != clr2.bin_size()) {
      throw std::runtime_error(fmt::format(
          FMT_STRING(
              "cooler \"{}\" and \"{}\" have different resolutions ({}  and {} respectively)"),
          clr1.uri(), clr2.uri(), clr1.bin_size(), clr2.bin_size()));
    }
  }
  return clr1.bin_size();
}

[[nodiscard]] inline ChromosomeSet get_chromosomes_checked(const std::vector<File>& coolers) {
  assert(coolers.size() > 1);
  const auto& clr1 = coolers.front();

  for (std::size_t i = 1; i < coolers.size(); ++i) {
    const auto& clr2 = coolers[i];
    if (clr1.chromosomes() != clr2.chromosomes()) {
      throw std::runtime_error(
          fmt::format(FMT_STRING("cooler \"{}\" and \"{}\" use different reference genomes"),
                      clr1.uri(), clr2.uri()));
    }
  }
  return clr1.chromosomes();
}

[[nodiscard]] inline bool merging_requires_float_pixels(const std::vector<File>& coolers) {
  for (const auto& clr : coolers) {
    if (clr.has_float_pixels()) {
      return true;
    }
  }
  return false;
}

}  // namespace internal

template <typename Str>
inline void merge(Str first_file, Str last_file, std::string_view dest_uri,
                  bool overwrite_if_exists, MergeStrategy strat) {
  static_assert(std::is_constructible_v<std::string, decltype(*first_file)>);

  std::vector<File> clrs{};
  std::transform(first_file, last_file, std::back_inserter(clrs),
                 [&](const auto& uri) { return File::open_read_only_read_once(std::string{uri}); });

  if (clrs.size() < 2) {
    throw std::runtime_error("unable to merge less than 2 coolers");
  }

  const auto chroms = internal::get_chromosomes_checked(clrs);
  const auto bin_size = internal::get_bin_size_checked(clrs);
  const auto float_pixels = internal::merging_requires_float_pixels(clrs);

  auto dest =
      float_pixels
          ? File::create_new_cooler<double>(dest_uri, chroms, bin_size, overwrite_if_exists)
          : File::create_new_cooler<std::int32_t>(dest_uri, chroms, bin_size, overwrite_if_exists);
  try {
    if (float_pixels) {
      internal::PixelMerger<double>(clrs).merge(dest, strat);
    } else {
      internal::PixelMerger<std::int32_t>(clrs).merge(dest, strat);
    }
  } catch (const std::exception& e) {
    throw std::runtime_error(
        fmt::format(FMT_STRING("failed to merge {} cooler files: {}"), clrs.size(), e.what()));
  }
}
}  // namespace coolerpp::utils
