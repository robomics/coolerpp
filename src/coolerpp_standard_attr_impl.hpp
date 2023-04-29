// Copyright (C) 2022 Roberto Rossini <roberros@uio.no>
//
// SPDX-License-Identifier: MIT

#pragma once

#include <cassert>
#include <cstdint>
#include <type_traits>

namespace coolerpp {

template <class PixelT, class>
inline StandardAttributes StandardAttributes::init(std::uint32_t bin_size_) {
  StandardAttributes attrs{};
  attrs.bin_size = bin_size_;
  if constexpr (std::is_floating_point_v<PixelT>) {
    attrs.sum = 0.0;
    attrs.cis = 0.0;
  }
  if constexpr (std::is_integral_v<PixelT>) {
    attrs.sum = std::int64_t(0);
    attrs.cis = std::int64_t(0);
  }
  return attrs;
}

inline StandardAttributes StandardAttributes::init_empty() noexcept {
  StandardAttributes attrs{};

  attrs.bin_type.reset();
  attrs.creation_date.reset();
  attrs.format_url.reset();
  attrs.generated_by.reset();
  attrs.assembly.reset();
  attrs.nbins.reset();
  attrs.nchroms.reset();
  attrs.metadata.reset();
  attrs.storage_mode.reset();
  attrs.sum.reset();
  attrs.cis.reset();

  return attrs;
}

}  // namespace coolerpp
