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
#include <utility>

#include "coolerpp/common.hpp"

namespace coolerpp {

template <class ChromosomeIt>
inline ChromosomeSet::ChromosomeSet(ChromosomeIt first_chrom, ChromosomeIt last_chrom)
    : _set(construct_set(first_chrom, last_chrom)) {}

template <class ChromosomeNameIt, class ChromosomeSizeIt>
inline ChromosomeSet::ChromosomeSet(ChromosomeNameIt first_chrom_name,
                                    ChromosomeNameIt last_chrom_name,
                                    ChromosomeSizeIt first_chrom_size)
    : _set(construct_set(first_chrom_name, last_chrom_name, first_chrom_size)) {}

template <class ChromosomeNameIt, class ChromosomeSizeIt>
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

template <class ChromosomeIt>
auto ChromosomeSet::construct_set(ChromosomeIt first_chrom, ChromosomeIt last_chrom) -> SetT {
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

}  // namespace coolerpp

constexpr auto fmt::formatter<coolerpp::Chromosome>::parse(format_parse_context& ctx)
    -> decltype(ctx.begin()) {
  if (ctx.begin() != ctx.end() && *ctx.begin() != '}') {
    throw fmt::format_error("invalid format");
  }
  return ctx.end();
}

template <class FormatContext>
inline auto fmt::formatter<coolerpp::Chromosome>::format(const coolerpp::Chromosome& c,
                                                         FormatContext& ctx) const
    -> decltype(ctx.out()) {
  return fmt::format_to(ctx.out(), FMT_STRING("{}:{}"), c.name, c.size);
}
