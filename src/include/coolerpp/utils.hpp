// Copyright (C) 2022 Roberto Rossini <roberros@uio.no>
//
// SPDX-License-Identifier: MIT

#pragma once

#include <cstdint>
#include <string_view>
#include <utility>
#include <vector>

#include "coolerpp/coolerpp.hpp"

namespace coolerpp::utils {

[[nodiscard]] bool equal(std::string_view uri1, std::string_view uri2,
                         bool ignore_attributes = true);
[[nodiscard]] bool equal(const File& clr1, const File& clr2, bool ignore_attributes = true);


}  // namespace coolerpp::utils

#include "../../utils_equal_impl.hpp"
