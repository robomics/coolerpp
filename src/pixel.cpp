// Copyright (C) 2022 Roberto Rossini <roberros@uio.no>
//
// SPDX-License-Identifier: MIT

#include "coolerpp/pixel.hpp"

#include <cassert>
#include <cstdint>
#include <string_view>

#include "coolerpp/bin_table.hpp"
#include "coolerpp/chromosome.hpp"

namespace coolerpp {

PixelCoordinates::PixelCoordinates(const BinTableLazy &bins, const Chromosome &chrom1_,
                                   const Chromosome &chrom2_, std::uint32_t bin1_start_,
                                   std::uint32_t bin2_start_)
    : _bins(&bins),
      chrom1(&chrom1_),
      chrom2(&chrom2_),
      bin1_start(bin1_start_),
      bin2_start(bin2_start_) {}

PixelCoordinates::PixelCoordinates(const BinTableLazy &bins, std::string_view chrom1_name,
                                   std::string_view chrom2_name, std::uint32_t bin1_start_,
                                   std::uint32_t bin2_start_)
    : PixelCoordinates(bins, bins.chromosomes().at(chrom1_name), bins.chromosomes().at(chrom2_name),
                       bin1_start_, bin2_start_) {}

PixelCoordinates::PixelCoordinates(const BinTableLazy &bins, std::uint32_t chrom1_id,
                                   std::uint32_t chrom2_id, std::uint32_t bin1_start_,
                                   std::uint32_t bin2_start_)
    : PixelCoordinates(bins, bins.chromosomes().at(chrom1_id), bins.chromosomes().at(chrom2_id),
                       bin1_start_, bin2_start_) {}

PixelCoordinates::PixelCoordinates(const BinTableLazy &bins, const Chromosome &chrom,
                                   std::uint32_t bin1_start_, std::uint32_t bin2_start_)
    : PixelCoordinates(bins, chrom, chrom, bin1_start_, bin2_start_) {}

PixelCoordinates::PixelCoordinates(const BinTableLazy &bins, std::uint32_t chrom_id,
                                   std::uint32_t bin1_start_, std::uint32_t bin2_start_)
    : PixelCoordinates(bins, chrom_id, chrom_id, bin1_start_, bin2_start_) {}

PixelCoordinates::PixelCoordinates(const BinTableLazy &bins, std::string_view chrom_name,
                                   std::uint32_t bin1_start_, std::uint32_t bin2_start_)
    : PixelCoordinates(bins, chrom_name, chrom_name, bin1_start_, bin2_start_) {}

PixelCoordinates::PixelCoordinates(const BinTableLazy &bins, std::uint64_t bin1_id,
                                   std::uint64_t bin2_id)
    : _bins(&bins) {
  auto coord1 = bins.bin_id_to_coords(bin1_id);
  auto coord2 = bins.bin_id_to_coords(bin2_id);

  chrom1 = &coord1.chrom;
  bin1_start = coord1.bin_start;

  chrom2 = &coord2.chrom;
  bin2_start = coord2.bin_start;
}

bool PixelCoordinates::operator==(const PixelCoordinates &other) const {
  if (this->_bins == other._bins) {
    // If pixels refer to the same bin table we can get away with comparing ptr to chromosomes
    // instead of the actual chromosomes
    return this->chrom1 == other.chrom1 && this->chrom2 == other.chrom2 &&
           this->bin1_start == other.bin1_start && this->bin2_start == other.bin2_start;
  }

  return *this->_bins == *other._bins && *this->chrom1 == *other.chrom1 &&
         *this->chrom2 == *other.chrom2 && this->bin1_start == other.bin1_start &&
         this->bin2_start == other.bin2_start;
}

bool PixelCoordinates::operator<(const PixelCoordinates &other) const {
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

std::uint64_t PixelCoordinates::bin1_id() const {
  assert(this->_bins);
  assert(this->chrom1);
  return this->_bins->coord_to_bin_id(*this->chrom1, this->bin1_start);
}

std::uint64_t PixelCoordinates::bin2_id() const {
  assert(this->_bins);
  assert(this->chrom2);
  return this->_bins->coord_to_bin_id(*this->chrom2, this->bin2_start);
}

std::uint32_t PixelCoordinates::bin_size() const noexcept {
  assert(this->_bins);
  return this->_bins->bin_size();
}

}  // namespace coolerpp
