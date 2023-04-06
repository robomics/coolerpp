// Copyright (C) 2022 Roberto Rossini <roberros@uio.no>
//
// SPDX-License-Identifier: MIT

#pragma once

#include <fmt/format.h>

#include <algorithm>
#include <cassert>
#include <cstdint>
#include <iterator>
#include <stdexcept>
#include <string>
#include <string_view>
#include <utility>

#include "coolerpp/common.hpp"

namespace coolerpp {

inline Chromosome::Chromosome(std::string name_, std::uint32_t size_) noexcept
    : name(std::move(name_)), size(size_) {}

inline bool Chromosome::operator<(const Chromosome& other) const noexcept {
  return this->name < other.name;
}

inline bool Chromosome::operator>(const Chromosome& other) const noexcept {
  return this->name > other.name;
}

inline bool Chromosome::operator<=(const Chromosome& other) const noexcept {
  return this->name <= other.name;
}

inline bool Chromosome::operator>=(const Chromosome& other) const noexcept {
  return this->name >= other.name;
}

inline bool Chromosome::operator==(const Chromosome& other) const noexcept {
  // Comparing by ID only is not enough, as different sets of chromosomes may use the same ID for
  // different chromosomes
  return this->name == other.name && this->size == other.size;
}

inline bool Chromosome::operator!=(const Chromosome& other) const noexcept {
  return !(*this == other);
}

inline bool operator<(const Chromosome& a, std::string_view b_name) noexcept {
  return a.name < b_name;
}
inline bool operator>(const Chromosome& a, std::string_view b_name) noexcept {
  return a.name > b_name;
}
inline bool operator<=(const Chromosome& a, std::string_view b_name) noexcept {
  return a.name <= b_name;
}
inline bool operator>=(const Chromosome& a, std::string_view b_name) noexcept {
  return a.name >= b_name;
}
inline bool operator==(const Chromosome& a, std::string_view b_name) noexcept {
  return a.name == b_name;
}
inline bool operator!=(const Chromosome& a, std::string_view b_name) noexcept {
  return a.name != b_name;
}

inline bool operator<(std::string_view a_name, const Chromosome& b) noexcept {
  return a_name < b.name;
}
inline bool operator>(std::string_view a_name, const Chromosome& b) noexcept {
  return a_name > b.name;
}
inline bool operator<=(std::string_view a_name, const Chromosome& b) noexcept {
  return a_name <= b.name;
}
inline bool operator>=(std::string_view a_name, const Chromosome& b) noexcept {
  return a_name >= b.name;
}
inline bool operator==(std::string_view a_name, const Chromosome& b) noexcept {
  return a_name == b.name;
}
inline bool operator!=(std::string_view a_name, const Chromosome& b) noexcept {
  return a_name != b.name;
}

template <typename ChromosomeIt>
inline ChromosomeSet::ChromosomeSet(ChromosomeIt first_chrom, ChromosomeIt last_chrom)
    : _set(construct_set(first_chrom, last_chrom)) {}

template <typename ChromosomeNameIt, typename ChromosomeSizeIt>
inline ChromosomeSet::ChromosomeSet(ChromosomeNameIt first_chrom_name,
                                    ChromosomeNameIt last_chrom_name,
                                    ChromosomeSizeIt first_chrom_size)
    : _set(construct_set(first_chrom_name, last_chrom_name, first_chrom_size)) {}

inline ChromosomeSet::ChromosomeSet(std::initializer_list<Chromosome> chromosomes)
    : ChromosomeSet(chromosomes.begin(), chromosomes.end()) {}

inline auto ChromosomeSet::begin() const -> const_iterator { return this->cbegin(); }
inline auto ChromosomeSet::end() const -> const_iterator { return this->cend(); }
inline auto ChromosomeSet::cbegin() const -> const_iterator { return this->_set.cbegin(); }
inline auto ChromosomeSet::cend() const -> const_iterator { return this->_set.cend(); }

inline auto ChromosomeSet::rbegin() const -> const_reverse_iterator { return this->rcbegin(); }
inline auto ChromosomeSet::rend() const -> const_reverse_iterator { return this->rcend(); }
inline auto ChromosomeSet::rcbegin() const -> const_reverse_iterator {
  return this->_set.rcbegin();
}
inline auto ChromosomeSet::rcend() const -> const_reverse_iterator { return this->_set.rcend(); }

inline bool ChromosomeSet::empty() const noexcept { return this->size() == 0; }
inline std::size_t ChromosomeSet::size() const noexcept { return this->_set.size(); }

inline auto ChromosomeSet::find(std::uint32_t id) const -> const_iterator {
  if (static_cast<std::size_t>(id) > this->size()) {
    return this->end();
  }
  return this->_set.begin() + static_cast<std::ptrdiff_t>(id);
}

inline auto ChromosomeSet::find(std::string_view chrom_name) const -> const_iterator {
  return this->_set.find(chrom_name);
}

inline auto ChromosomeSet::find(const Chromosome& chrom) const -> const_iterator {
  return this->_set.find(chrom);
}

inline const Chromosome& ChromosomeSet::at(std::uint32_t id) const {
  this->validate_chrom_id(id);
  return *this->find(id);
}

inline const Chromosome& ChromosomeSet::at(std::string_view chrom_name) const {
  if (const auto match = this->find(chrom_name); match != this->end()) {
    return *match;
  }
  throw std::out_of_range(fmt::format(FMT_STRING("chromosome \"{}\" not found"), chrom_name));
}

inline const Chromosome& ChromosomeSet::operator[](std::uint32_t id) const noexcept {
  auto it = this->find(id);
  assert(it != this->end());
  return *it;
}
inline const Chromosome& ChromosomeSet::operator[](std::string_view chrom_name) const noexcept {
  auto it = this->find(chrom_name);
  assert(it != this->end());
  return *it;
}

inline bool ChromosomeSet::contains(std::uint32_t id) const {
  return this->find(id) != this->end();
}
inline bool ChromosomeSet::contains(const Chromosome& chrom) const {
  return this->find(chrom) != this->end();
}
inline bool ChromosomeSet::contains(std::string_view chrom_name) const {
  return this->find(chrom_name) != this->end();
}

inline std::uint32_t ChromosomeSet::get_id(const Chromosome& chrom) const {
  if (const auto match = this->find(chrom); match != this->end()) {
    return static_cast<std::uint32_t>(std::distance(this->begin(), match));
  }
  throw std::out_of_range(fmt::format(FMT_STRING("chromosome {} not found"), chrom));
}

inline std::uint32_t ChromosomeSet::get_id(std::string_view chrom_name) const {
  if (const auto match = this->find(chrom_name); match != this->end()) {
    return static_cast<std::uint32_t>(std::distance(this->begin(), match));
  }
  throw std::out_of_range(fmt::format(FMT_STRING("chromosome \"{}\" not found"), chrom_name));
}

inline bool ChromosomeSet::operator==(const ChromosomeSet& other) const {
  return this->_set == other._set;
}

inline bool ChromosomeSet::operator!=(const ChromosomeSet& other) const {
  return !(*this == other);
}

inline const Chromosome& ChromosomeSet::find_longest_chromosome() const {
  if (this->_set.empty()) {
    throw std::runtime_error("find_longest_chromosome() was called on an empty ChromosomeSet");
  }
  return *std::max_element(this->_set.begin(), this->_set.end(),
                           [&](const Chromosome& chrom1, const Chromosome& chrom2) {
                             return chrom1.size < chrom2.size;
                           });
}
inline const Chromosome& ChromosomeSet::find_chromosome_with_longest_name() const {
  if (this->_set.empty()) {
    throw std::runtime_error(
        "find_chromosome_with_longest_name() was called on an empty ChromosomeSet");
  }
  return *std::max_element(this->_set.begin(), this->_set.end(),
                           [&](const Chromosome& chrom1, const Chromosome& chrom2) {
                             return chrom1.name.size() < chrom2.name.size();
                           });
}

inline void ChromosomeSet::validate_chrom_id(std::uint32_t chrom_id) const {
  if (static_cast<std::size_t>(chrom_id) >= this->size()) {
    throw std::out_of_range(fmt::format(FMT_STRING("chromosome with id {} not found"), chrom_id));
  }
}

template <typename ChromosomeNameIt, typename ChromosomeSizeIt>
inline auto ChromosomeSet::construct_set(ChromosomeNameIt first_chrom_name,
                                         ChromosomeNameIt last_chrom_name,
                                         ChromosomeSizeIt first_chrom_size) -> SetT {
  const auto num_chroms =
      conditional_static_cast<std::size_t>(std::distance(first_chrom_name, last_chrom_name));
  SetT set(num_chroms);

  std::transform(
      first_chrom_name, last_chrom_name, std::inserter(set, set.begin()), [&](const auto& name) {
        auto chrom = Chromosome{std::string{name}, *(first_chrom_size++)};
        if (const auto& collision = set.find(chrom); collision != set.end()) {
          const auto id1 = set.size();
          const auto id2 = std::distance(set.begin(), collision);
          throw std::runtime_error(fmt::format(
              FMT_STRING("found duplicate chromosome: {}:{} (id={}) collides with {}:{} (id={})"),
              chrom.name, chrom.size, id1, collision->name, collision->size, id2));
        }
        return chrom;
      });

  return set;
}

template <typename ChromosomeIt>
inline auto ChromosomeSet::construct_set(ChromosomeIt first_chrom, ChromosomeIt last_chrom)
    -> SetT {
  const auto num_chroms =
      conditional_static_cast<std::size_t>(std::distance(first_chrom, last_chrom));
  SetT set(num_chroms);

  std::transform(first_chrom, last_chrom, std::inserter(set, set.begin()), [&](Chromosome chrom) {
    if (const auto& collision = set.find(chrom); collision != set.end()) {
      const auto id1 = set.size();
      const auto id2 = std::distance(set.begin(), collision);
      throw std::runtime_error(fmt::format(
          FMT_STRING("found duplicate chromosome: {}:{} (id={}) collides with {}:{} (id={})"),
          chrom.name, chrom.size, id1, collision->name, collision->size, id2));
    }
    return chrom;
  });

  return set;
}

namespace internal {

inline bool ChromEq::operator()(const Chromosome& a, const Chromosome& b) const noexcept {
  return a.name == b.name;
}
inline bool ChromEq::operator()(const Chromosome& a, std::string_view b_name) const noexcept {
  return a.name == b_name;
}
inline bool ChromEq::operator()(std::string_view a_name, const Chromosome& b) const noexcept {
  return a_name == b.name;
}

inline std::size_t ChromHasher::operator()(const Chromosome& c) const {
  return std::hash<Chromosome>{}(c);
}
inline std::size_t ChromHasher::operator()(std::string_view name) const {
  return std::hash<std::string_view>{}(name);
}

}  // namespace internal
}  // namespace coolerpp

inline std::size_t std::hash<coolerpp::Chromosome>::operator()(
    const coolerpp::Chromosome& k) const {
  return std::hash<std::string>{}(k.name);
}

constexpr auto fmt::formatter<coolerpp::Chromosome>::parse(format_parse_context& ctx)
    -> decltype(ctx.begin()) {
  if (ctx.begin() != ctx.end() && *ctx.begin() != '}') {
    throw fmt::format_error("invalid format");
  }
  return ctx.end();
}

template <typename FormatContext>
inline auto fmt::formatter<coolerpp::Chromosome>::format(const coolerpp::Chromosome& c,
                                                         FormatContext& ctx) const
    -> decltype(ctx.out()) {
  return fmt::format_to(ctx.out(), FMT_STRING("{}:{}"), c.name, c.size);
}
