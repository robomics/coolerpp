// Copyright (C) 2023 Roberto Rossini <roberros@uio.no>
//
// SPDX-License-Identifier: MIT

#include "coolerpp/utils.hpp"
#include "coolerpp_tools/config.hpp"
#include "coolerpp_tools/tools.hpp"
#include <filesystem>

namespace coolerpp::tools {

void merge_subcmd(const MergeConfig& c) {
  utils::merge(c.input_uris.begin(), c.input_uris.end(), c.output_uri, c.force, c.chunk_size, false);
}

}  // namespace coolerpp::tools
