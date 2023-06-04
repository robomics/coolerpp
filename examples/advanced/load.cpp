// Copyright (C) 2023 Roberto Rossini <roberros@uio.no>
//
// SPDX-License-Identifier: MIT

#include <fmt/compile.h>
#include <fmt/format.h>

#include <algorithm>
#include <cassert>
#include <fstream>
#include <string_view>

#include "coolerpp/coolerpp.hpp"
#include "coolerpp_tools/config.hpp"
#include "coolerpp_tools/tools.hpp"

namespace coolerpp::tools {

template <typename N>
[[nodiscard]] static Pixel<N> parse_pixel(const BinTable& bins, std::string_view line) {
  auto next_token = [&]() {
    assert(!line.empty());
    const auto pos = line.find('\t');
    auto tok = line.substr(0, pos);
    line.remove_prefix(pos + 1);
    return tok;
  };

  const auto chrom1 = next_token();
  const auto start1 = internal::parse_numeric_or_throw<std::uint32_t>(next_token());
  std::ignore = next_token();

  const auto chrom2 = next_token();
  const auto start2 = internal::parse_numeric_or_throw<std::uint32_t>(next_token());
  std::ignore = next_token();

  const auto count = internal::parse_numeric_or_throw<N>(next_token());
  return Pixel<N>{{bins.at(chrom1, start1), bins.at(chrom2, start2)}, count};
}

[[nodiscard]] static ChromosomeSet import_chromosomes(std::string_view chrom_sizes) {
  try {
    std::string line;
    std::vector<Chromosome> chroms{};

    std::ifstream f(std::string{chrom_sizes});
    while (std::getline(f, line)) {
      const auto delim_pos = line.find('\t');

      const auto chrom_name = line.substr(0, delim_pos);
      const auto chrom_size = std::stoull(line.substr(delim_pos + 1));

      chroms.emplace_back(0, chrom_name, chrom_size);
    }

    return {chroms.begin(), chroms.end()};
  } catch (const std::exception& e) {
    throw std::runtime_error(
        fmt::format(FMT_STRING("an error occurred while importing chromosomes from {}: {}"),
                    chrom_sizes, e.what()));
  }
}

template <typename N>
[[nodiscard]] static bool process_batch(const BinTable& bins, std::size_t batch_size,
                                        std::vector<Pixel<N>>& buffer) {
  buffer.clear();
  std::string line;
  try {
    while (std::getline(std::cin, line)) {
      if (buffer.size() == batch_size) {
        return true;
      }
      buffer.emplace_back(parse_pixel<N>(bins, line));
    }
    return false;
  } catch (const std::exception& e) {
    throw std::runtime_error(
        fmt::format(FMT_STRING("encountered error while processing the following line:\n"
                               "\"{}\"\n"
                               "Cause: {}"),
                    line, e.what()));
  }
}

template <typename N, std::size_t chunk_size = std::numeric_limits<std::size_t>::max()>
static std::string ingest_pixels(coolerpp::File&& clr, std::size_t batch_size = 1'000'000) {
  batch_size = std::min(batch_size, chunk_size);
  const auto& bins = clr.bins();

  std::vector<Pixel<N>> write_buffer(batch_size);

  std::size_t lines_processed = 0;
  while (process_batch(bins, batch_size, write_buffer)) {
    clr.append_pixels(write_buffer.begin(), write_buffer.end());
    lines_processed += write_buffer.size();

    if (lines_processed % std::max(batch_size, std::size_t(1'000'000)) == 0) {
      fmt::print(stderr, FMT_STRING("Read {}M pixels...\n"), lines_processed / 1'000'000);
    }
    if (lines_processed >= chunk_size) {
      break;
    }
  }

  if (!write_buffer.empty()) {
    clr.append_pixels(write_buffer.begin(), write_buffer.end());
    lines_processed += write_buffer.size();
  }
  if (lines_processed == 0) {
    return "";
  }
  return clr.uri();
}

void load_subcmd(const LoadConfig& c) {
  auto chroms = import_chromosomes(c.path_to_chrom_sizes);
  if (c.assume_sorted) {
    c.count_as_float
        ? ingest_pixels<double>(File::create_new_cooler<double>(c.uri, chroms, c.bin_size, c.force))
        : ingest_pixels<std::int32_t>(
              File::create_new_cooler<std::int32_t>(c.uri, chroms, c.bin_size, c.force));
    return;
  }

  std::vector<std::string> uris{};
  constexpr std::size_t chunk_size = 5'000'000;

  for (std::size_t i = 0; true; ++i) {
    const auto tmp_uri = fmt::format(FMT_STRING("/tmp/foo_{}.cool"), i);
    uris.emplace_back(
        c.count_as_float
            ? ingest_pixels<double, chunk_size>(
                  File::create_new_cooler<double>(tmp_uri, chroms, c.bin_size, c.force))
            : ingest_pixels<std::int32_t, chunk_size>(
                  File::create_new_cooler<std::int32_t>(tmp_uri, chroms, c.bin_size, c.force)));
    if (uris.back().empty()) {
      uris.pop_back();
      break;
    }
    fmt::print(stderr, FMT_STRING("Done writing to tmp file {}...\n"), tmp_uri);
  }

  merge_subcmd(MergeConfig{uris, c.uri, c.count_as_float, c.force});
  for (const auto& f : uris) {
    std::filesystem::remove(f);
  }
}

}  // namespace coolerpp::tools
