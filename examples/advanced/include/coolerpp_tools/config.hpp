// Copyright (C) 2023 Roberto Rossini <roberros@uio.no>
//
// SPDX-License-Identifier: MIT

#pragma once

#include <cstdint>
#include <filesystem>
#include <string>
#include <variant>

#include "coolerpp/coolerpp.hpp"

namespace coolerpp::tools {

struct DumpConfig {
  std::string uri{};

  std::string range1{"all"};
  std::string range2{"all"};
  std::filesystem::path query_file{};

  std::string table{"pixels"};
  bool join{true};

  std::string balanced{};
  Weights::Type weight_type{Weights::Type::INFER};
};

struct LoadConfig {
  std::string uri{};

  std::filesystem::path path_to_chrom_sizes{};
  std::uint32_t bin_size{};
  std::string format{"bg2"};
  std::string assembly{"unknown"};
  bool count_as_float{false};
  bool assume_sorted{true};
  bool force{false};
};

struct MergeConfig {
  std::vector<std::string> input_uris{};
  std::string output_uri{};

  bool floating_point{false};
  bool force{false};
};

// clang-format off
using Config = std::variant<std::monostate,
                            DumpConfig,
                            LoadConfig,
                            MergeConfig>;
// clang-format on

}  // namespace coolerpp::tools
