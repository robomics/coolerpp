// Copyright (C) 2022 Roberto Rossini <roberros@uio.no>
//
// SPDX-License-Identifier: MIT

#pragma once

#include <array>  // for array
#include <atomic>
#include <functional>  // for hash
#include <memory>
#include <string>
#include <string_view>
#include <type_traits>  // for is_arithmetic_v

#include "coolerpp/chromosome.hpp"

namespace coolerpp {

class PixelCoordinates {
  const BinTableLazy *_bins{};

 public:
  const Chromosome *chrom1{};
  const Chromosome *chrom2{};

  std::uint32_t bin1_start{};
  std::uint32_t bin2_start{};

  PixelCoordinates() = delete;
  inline PixelCoordinates(const BinTableLazy &bins, std::uint32_t chrom1_id_,
                          std::uint32_t chrom2_id_, std::uint32_t bin1_start_,
                          std::uint32_t bin2_start_)
      : _bins(&bins),
        chrom1(&bins.chromosomes().at(chrom1_id_)),
        chrom2(&bins.chromosomes().at(chrom2_id_)),
        bin1_start(bin1_start_),
        bin2_start(bin2_start_) {}

  inline PixelCoordinates(const BinTableLazy &bins, std::uint32_t chrom_id,
                          std::uint32_t bin1_start_, std::uint32_t bin2_start_)
      : PixelCoordinates(bins, chrom_id, chrom_id, bin1_start_, bin2_start_) {}

  inline PixelCoordinates(const BinTableLazy &bins, std::string_view chrom1_name,
                          std::string_view chrom2_name, std::uint32_t bin1_start_,
                          std::uint32_t bin2_start_)
      : _bins(&bins),
        chrom1(&bins.chromosomes().at(chrom1_name)),
        chrom2(&bins.chromosomes().at(chrom2_name)),
        bin1_start(bin1_start_),
        bin2_start(bin2_start_) {}

  inline PixelCoordinates(const BinTableLazy &bins, std::string_view chrom_name,
                          std::uint32_t bin1_start_, std::uint32_t bin2_start_)
      : PixelCoordinates(bins, chrom_name, chrom_name, bin1_start_, bin2_start_) {}

  inline PixelCoordinates(const BinTableLazy &bins, std::uint64_t bin1_id, std::uint64_t bin2_id)
      : _bins(&bins) {
    auto coord1 = bins.bin_id_to_coords(bin1_id);
    auto coord2 = bins.bin_id_to_coords(bin2_id);

    chrom1 = &coord1.chrom;
    bin1_start = coord1.bin_start;

    chrom2 = &coord2.chrom;
    bin2_start = coord2.bin_start;
  }
  /*
        chrom1(&bins.chromosomes().at(chrom1_id_)),
        chrom2(&bins.chromosomes().at(chrom2_id_)),
        bin1_start(bin1_start_),
        bin2_start(bin2_start_) {}
*/
  [[nodiscard]] constexpr bool operator==(const PixelCoordinates &other) const {
    return this->_bins == other._bins && this->chrom1 == other.chrom1 &&
           this->chrom2 == other.chrom2 && this->bin1_start == other.bin1_start &&
           this->bin2_start == other.bin2_start;
  }

  [[nodiscard]] constexpr bool operator<(const PixelCoordinates &other) const {
    assert(this->_bins == other._bins);
    assert(this->chrom1 && this->chrom2);

    if (this->chrom1 != other.chrom1 && !!this->chrom1 && !!other.chrom2) {
      return *this->chrom1 < *other.chrom1;
    }
    if (this->chrom2 != other.chrom2 && !!this->chrom2 && !!other.chrom2) {
      return *this->chrom2 < *other.chrom2;
    }
    if (this->bin1_start != other.bin1_start) {
      return this->bin1_start < other.bin1_start;
    }
    return this->bin2_start < other.bin2_start;
  }

  [[nodiscard]] inline std::uint64_t bin1_id() const {
    assert(this->_bins);
    assert(this->chrom1);
    return this->_bins->coord_to_bin_id(*this->chrom1, this->bin1_start);
  }
  [[nodiscard]] inline std::uint64_t bin2_id() const {
    assert(this->_bins);
    assert(this->chrom2);
    return this->_bins->coord_to_bin_id(*this->chrom2, this->bin2_start);
  }
};

template <class N>
struct Pixel {
  static_assert(std::is_arithmetic_v<N>);
  using Coordinates = PixelCoordinates;

  Coordinates coords;
  N count;

  [[nodiscard]] constexpr bool operator==(const Pixel &other) const {
    return this->coords == other.coords && this->count == other.count;
  }
  [[nodiscard]] constexpr bool operator<(const Pixel &other) const {
    return this->coords < other.coords;
  }
};

}  // namespace coolerpp
