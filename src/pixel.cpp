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

PixelCoordinates::PixelCoordinates(const BinTableLazy &bins, const Chromosome &chrom1,
                                   const Chromosome &chrom2, std::uint32_t bin1_start_,
                                   std::uint32_t bin2_start_)
    : PixelCoordinates(bins, bins.chromosomes().get_id(chrom1), bins.chromosomes().get_id(chrom2),
                       bin1_start_, bin2_start_) {}

PixelCoordinates::PixelCoordinates(const BinTableLazy &bins, std::string_view chrom1_name,
                                   std::string_view chrom2_name, std::uint32_t bin1_start_,
                                   std::uint32_t bin2_start_)
    : PixelCoordinates(bins, bins.chromosomes().get_id(chrom1_name),
                       bins.chromosomes().get_id(chrom2_name), bin1_start_, bin2_start_) {}

PixelCoordinates::PixelCoordinates(const BinTableLazy &bins, std::uint32_t chrom1_id_,
                                   std::uint32_t chrom2_id_, std::uint32_t bin1_start_,
                                   std::uint32_t bin2_start_)
    : _bins(&bins),
      chrom1_id(chrom1_id_),
      chrom2_id(chrom2_id_),
      bin1_start(bin1_start_),
      bin2_start(bin2_start_) {}

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

  chrom1_id = bins.chromosomes().get_id(coord1.chrom);
  chrom2_id = bins.chromosomes().get_id(coord2.chrom);

  bin1_start = coord1.bin_start;
  bin2_start = coord2.bin_start;
}

const Chromosome &PixelCoordinates::chrom1() const {
  assert(this->_bins);
  assert(!!*this);
  return this->_bins->chromosomes()[chrom1_id];
}

const Chromosome &PixelCoordinates::chrom2() const {
  assert(this->_bins);
  assert(!!*this);
  return this->_bins->chromosomes()[chrom2_id];
}

std::uint64_t PixelCoordinates::bin1_id() const {
  assert(this->_bins);
  assert(!!*this);
  return this->_bins->coord_to_bin_id(this->chrom1_id, this->bin1_start);
}

std::uint64_t PixelCoordinates::bin2_id() const {
  assert(this->_bins);
  assert(!!*this);
  return this->_bins->coord_to_bin_id(this->chrom2_id, this->bin2_start);
}

std::uint32_t PixelCoordinates::bin_size() const noexcept {
  assert(this->_bins);
  return this->_bins->bin_size();
}

}  // namespace coolerpp
