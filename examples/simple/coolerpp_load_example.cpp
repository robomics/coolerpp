// Copyright (C) 2022 Roberto Rossini <roberros@uio.no>
//
// SPDX-License-Identifier: MIT

#include <fmt/format.h>
#include <fmt/std.h>

#include <array>
#include <cassert>
#include <chrono>
#include <cstdint>
#include <fstream>
#include <string>
#include <string_view>

#include "coolerpp/coolerpp.hpp"

using namespace coolerpp;

struct BEDPE {
  std::string chrom1{};
  std::uint32_t start1{};
  std::uint32_t end1{};
  std::string chrom2{};
  std::uint32_t start2{};
  std::uint32_t end2{};

  std::int32_t count{};

  /// Construct a BEDPE object from a tab-separated line
  //
  // Example line:
  // chr1	5000	10000	chr1	87650000	87655000	1
  explicit BEDPE(std::string_view line, std::string_view delim = "\t") {
    auto next_token = [&]() {
      assert(!line.empty());
      const auto pos = line.find(delim);
      auto tok = line.substr(0, pos);
      line.remove_prefix(pos + 1);
      return std::string{tok};
    };

    chrom1 = next_token();
    start1 = static_cast<std::uint32_t>(std::stoul(next_token()));
    end1 = static_cast<std::uint32_t>(std::stoul(next_token()));

    chrom2 = next_token();
    start2 = static_cast<std::uint32_t>(std::stoul(next_token()));
    end2 = static_cast<std::uint32_t>(std::stoul(next_token()));

    count = static_cast<std::int32_t>(std::stol(next_token()));
  }

  [[nodiscard]] Pixel<std::int32_t> to_pixel(
      const std::shared_ptr<const BinTableLazy>& bins) const {
    return {PixelCoordinates{bins, this->chrom1, this->chrom2, this->start1, this->start2},
            this->count};
  }
};

/// Parse a 2-column .chrom.sizes file
//
// Example:
// chr1    248956422
// chr2    242193529
// ...
[[nodiscard]] static ChromosomeSet import_chromosomes(std::string_view chrom_sizes,
                                                      std::string_view delim = "\t") {
  try {
    std::string line;
    std::vector<Chromosome> chroms{};

    std::ifstream f(std::string{chrom_sizes});
    while (std::getline(f, line)) {
      const auto delim_pos = line.find(delim);

      const auto chrom_name = line.substr(0, delim_pos);
      const auto chrom_size = std::stoull(line.substr(delim_pos + 1));

      chroms.emplace_back(chrom_name, chrom_size);
    }

    return {chroms.begin(), chroms.end()};
  } catch (const std::exception& e) {
    throw std::runtime_error(
        fmt::format(FMT_STRING("an error occurred while importing chromosomes from {}: {}"),
                    chrom_sizes, e.what()));
  }
}

/// Process a chunk of pixels and return true as long as there's more data to be processed
[[nodiscard]] static bool process_chunk(const std::shared_ptr<const BinTableLazy>& bins,
                                        std::size_t batch_size,
                                        std::vector<Pixel<std::int32_t>>& buffer) {
  buffer.clear();
  std::string line;
  try {
    while (std::getline(std::cin, line)) {
      if (buffer.size() == batch_size) {
        return true;
      }
      buffer.emplace_back(BEDPE(line).to_pixel(bins));
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

static std::size_t ingest_pixels(std::string_view path_to_chrom_sizes,
                                 std::string_view path_to_output_cooler, std::uint32_t bin_size,
                                 std::size_t batch_size = 100'000) {
  auto cooler = File::create_new_cooler(path_to_output_cooler,
                                        import_chromosomes(path_to_chrom_sizes), bin_size);
  const auto& bins = cooler.bins_ptr();

  std::vector<Pixel<std::int32_t>> write_buffer;
  write_buffer.reserve(batch_size);

  std::size_t lines_processed = 0;
  while (process_chunk(bins, batch_size, write_buffer)) {
    lines_processed += write_buffer.size();
    cooler.append_pixels(write_buffer.begin(), write_buffer.end());

    if (lines_processed % 1'000'000 == 0) {
      fmt::print(stderr, FMT_STRING("Read {}M pixels...\n"), lines_processed / 1'000'000);
    }
  }

  if (!write_buffer.empty()) {
    cooler.append_pixels(write_buffer.begin(), write_buffer.end());
    lines_processed += write_buffer.size();
  }

  return lines_processed;
}

static void print_usage(std::string_view arg0) {
  fmt::print(
      stderr,
      FMT_STRING(
          "Usage:   {0} my_chroms.chrom.sizes bin_size path/to/output.cool < contacts.bedpe\n"
          "Example: {0} test/data/hg38.chrom.sizes 1000 /tmp/output.cool < contacts.bedpe\n"
          "Example: zcat contacts.bedpe.gz | {0} test/data/hg38.chrom.sizes 1000 "
          "/tmp/output.cool\n"),
      arg0);
}

int main(int argc, char** argv) {
  if (argc != 4) {
    std::string_view arg0 = argv[0];                   // NOLINT
    std::string_view arg1 = argc == 1 ? "" : argv[1];  // NOLINT
    print_usage(arg0);                                 // NOLINT
    if (arg1 == "--help" || arg1 == "-h") {
      return 0;
    }
    return 1;
  }

  try {
    const std::string_view path_to_chrom_sizes = argv[1];                   // NOLINT
    const auto bin_size = static_cast<std::uint32_t>(std::stoul(argv[2]));  // NOLINT
    const std::string_view path_to_output_cooler = argv[3];                 // NOLINT

    const auto t0 = std::chrono::steady_clock::now();
    ingest_pixels(path_to_chrom_sizes, path_to_output_cooler, bin_size);
    const auto t1 = std::chrono::steady_clock::now();

    const auto elapsed_time_ms =
        std::chrono::duration_cast<std::chrono::milliseconds>(t1 - t0).count();

    const auto nnz = File::open_read_only(path_to_output_cooler).attributes().nnz;
    assert(nnz.has_value());
    fmt::print(stderr, FMT_STRING("Written {} pixels in {}s!\n"), *nnz,
               static_cast<double>(elapsed_time_ms) / 1000.0);
  } catch (const std::exception& e) {
    fmt::print(stderr, FMT_STRING("The following error occurred while running coolerpp_load: {}\n"),
               e.what());
    return 1;
  }
}
