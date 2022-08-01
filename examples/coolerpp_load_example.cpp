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

/// Parse a line containing a BEDPE record
//
// Example line:
// 1	5000	10000	1	87650000	87655000	1
void parse_bedpe_record(std::string_view line, std::array<std::string, 7>& buff,
                        std::string_view delim = "\t") {
  assert(!line.empty());
  for (auto& tok : buff) {
    const auto pos = line.find(delim);
    tok = line.substr(0, pos);
    line.remove_prefix(pos + 1);
  }
}

[[nodiscard]] ChromosomeSet import_chromosomes(const std::filesystem::path& chrom_sizes,
                                               std::string_view delim = "\t") {
  try {
    std::string line;
    std::vector<Chromosome> chroms{};

    std::ifstream f(chrom_sizes);
    while (std::getline(f, line)) {
      const auto delim_pos = line.find(delim);

      std::string chrom_name = line.substr(0, delim_pos);
      std::string chrom_size_str = line.substr(delim_pos + 1);

      chroms.emplace_back(chrom_name, std::stoull(chrom_size_str));
    }

    return {chroms.begin(), chroms.end()};
  } catch (const std::exception& e) {
    throw std::runtime_error(
        fmt::format(FMT_STRING("an error occurred while importing chromosomes from {}: {}"),
                    chrom_sizes, e.what()));
  }
}

[[nodiscard]] Pixel<std::uint32_t> construct_pixel(const BinTableLazy& bins,
                                                   const std::array<std::string, 7>& bedpe_tokens) {
  const auto& chrom1_name = bedpe_tokens[0];
  const auto& chrom2_name = bedpe_tokens[3];
  const auto bin1_start = static_cast<std::uint32_t>(std::stoul(bedpe_tokens[1]));
  const auto bin2_start = static_cast<std::uint32_t>(std::stoul(bedpe_tokens[4]));

  PixelCoordinates coords{bins, chrom1_name, chrom2_name, bin1_start, bin2_start};
  const auto count = static_cast<std::uint32_t>(std::stoul(bedpe_tokens[6]));

  return {coords, count};
}

template <class ContactT>
void ingest_pixels(const std::string& path_to_chrom_sizes, const std::string& path_to_output_cooler,
                   std::uint32_t bin_size, std::size_t batch_size = 100'000) {
  auto chromosomes = import_chromosomes(path_to_chrom_sizes);
  auto cooler = File::create_new_cooler<ContactT>(path_to_output_cooler, chromosomes, bin_size);
  const auto& bins = cooler.bins();

  std::vector<Pixel<ContactT>> write_buffer;
  write_buffer.reserve(batch_size);

  std::size_t i = 0;
  std::string line;
  try {
    std::array<std::string, 7> bedpe_record{};
    while (std::getline(std::cin, line)) {
      if (write_buffer.size() == batch_size) {
        cooler.append_pixels(write_buffer.begin(), write_buffer.end());
        write_buffer.clear();
      }
      parse_bedpe_record(line, bedpe_record);
      write_buffer.emplace_back(construct_pixel(bins, bedpe_record));

      if (++i % 1'000'000 == 0) {
        fmt::print(stderr, FMT_STRING("Read {}M pixels from stdin...\n"), i / 1'000'000);
      }
    }
    if (!write_buffer.empty()) {
      cooler.append_pixels(write_buffer.begin(), write_buffer.end());
    }
  } catch (const std::exception& e) {
    throw std::runtime_error(
        fmt::format(FMT_STRING("encountered error while processing line {} of file {}:\n"
                               "{}\n"
                               "line that caused the error: \"{}\""),
                    i, path_to_chrom_sizes, e.what(), line));
  }
}

int main(int argc, char** argv) {
  try {
    std::ios_base::sync_with_stdio(false);

    if (argc != 4) {
      fmt::print(
          stderr,
          FMT_STRING(
              "Usage:   {0} my_chroms.chrom.sizes bin_size path/to/output.cool < contacts.bedpe\n"
              "Example: {0} test/data/hg38.chrom.sizes 1000 /tmp/output.cool < contacts.bedpe\n"
              "Example: cat contacts.bedpe | {0} test/data/hg38.chrom.sizes 1000 "
              "/tmp/output.cool\n"),
          argv[0]);  // NOLINT(cppcoreguidelines-pro-bounds-pointer-arithmetic)
      return 1;
    }

    // clang-format off
    const std::string path_to_chrom_sizes = argv[1];                        // NOLINT(cppcoreguidelines-pro-bounds-pointer-arithmetic)
    const auto bin_size = static_cast<std::uint32_t>(std::stoul(argv[2]));  // NOLINT(cppcoreguidelines-pro-bounds-pointer-arithmetic)
    const std::string path_to_output_cooler = argv[3];                      // NOLINT(cppcoreguidelines-pro-bounds-pointer-arithmetic)
    // clang-format on

    using ContactT = std::uint32_t;

    const auto t0 = std::chrono::steady_clock::now();
    ingest_pixels<ContactT>(path_to_chrom_sizes, path_to_output_cooler, bin_size);
    const auto t1 = std::chrono::steady_clock::now();

    const auto elapsed_time_ms =
        std::chrono::duration_cast<std::chrono::milliseconds>(t1 - t0).count();

    auto cooler = File::open_read_only(path_to_output_cooler);
    fmt::print(stderr, FMT_STRING("Written {} pixels in {}s!\n"), *cooler.attributes().nnz,
               static_cast<double>(elapsed_time_ms) / 1000.0);
  } catch (const std::exception& e) {
    fmt::print(stderr, FMT_STRING("The following error occurred while running coolerpp_load: {}\n"),
               e.what());
    return 1;
  }
}
