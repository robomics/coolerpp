// Copyright (C) 2022 Roberto Rossini <roberros@uio.no>
//
// SPDX-License-Identifier: MIT

#pragma once

#include <fmt/chrono.h>

#include <array>
#include <chrono>
#include <cstdint>
#include <highfive/H5DataSet.hpp>
#include <highfive/H5DataSpace.hpp>
#include <highfive/H5File.hpp>
#include <highfive/H5Group.hpp>
#include <optional>
#include <string>
#include <string_view>
#include <variant>
#include <vector>

#include "coolerpp/bin_table.hpp"
#include "coolerpp/chromosome.hpp"
#include "coolerpp/dataset.hpp"
#include "coolerpp/group.hpp"
#include "coolerpp/index.hpp"
#include "coolerpp/internal/numeric_variant.hpp"
#include "coolerpp/pixel.hpp"

namespace coolerpp {

enum class IO_MODE : std::uint_fast16_t {
  ReadOnly = HighFive::File::ReadOnly,
  ReadWrite = HighFive::File::ReadWrite,
  Truncate = HighFive::File::Truncate,
  Excl = HighFive::File::Excl,
  Debug = HighFive::File::Debug,
  Create = HighFive::File::Create,
  Overwrite = Truncate,
  OpenOrCreate = ReadWrite | Create
};

struct StandardAttributes {
  // Mandatory attributes
  std::uint32_t bin_size{0};
  std::optional<std::string> bin_type{"fixed"};  // Mandatory in v3
  std::string format{COOL_MAGIC};
  std::uint8_t format_version{3};
  std::optional<std::string> storage_mode{"symmetric-upper"};  // Mandatory in v3

  // Reserved attributes
  std::optional<std::string> creation_date{fmt::format(
      FMT_STRING("{:%FT%T}"), fmt::gmtime(std::time(nullptr)))};  // e.g. 2022-07-26T20:35:19
  std::optional<std::string> generated_by{"coolerpp-v0.0.1"};     // TODO fixme
  std::optional<std::string> assembly{"unknown"};
  std::optional<std::string> metadata{"{}"};

  // Optional but common
  std::optional<std::string> format_url{"https://github.com/open2c/cooler"};
  std::optional<std::uint64_t> nbins{0};
  std::optional<std::uint32_t> nchroms{0};
  std::optional<std::uint64_t> nnz{0};
  using SumVar = std::variant<double, std::int64_t, std::uint64_t>;
  SumVar sum{std::uint64_t(0)};

  template <class PixelT = std::uint32_t, class = std::enable_if_t<std::is_arithmetic_v<PixelT>>>
  [[nodiscard]] static StandardAttributes init(std::uint32_t bin_size_);
  [[nodiscard]] static StandardAttributes init_empty() noexcept;
};

template <class InputIt>
void init_mcool(std::string_view file_path, InputIt first_resolution, InputIt last_resolution,
                bool force_overwrite = false);
void init_mcool(std::string_view file_path, bool force_overwrite = false);

template <class ChromSizeInputIt, class CellIDInputIt>
void init_scool(std::string_view file_path, ChromSizeInputIt first_chrom,
                ChromSizeInputIt last_chrom, CellIDInputIt first_cell_id,
                CellIDInputIt last_cell_id, bool force_overwrite = false);
template <class InputIt>
void init_scool(std::string_view file_path, InputIt first_chrom, InputIt last_chrom,
                bool force_overwrite = false);

[[nodiscard]] Chromosome find_longest_chromosome(const ChromosomeSet &chroms);
[[nodiscard]] Chromosome find_chromosome_with_longest_name(const ChromosomeSet &chroms);

class File {
 public:
  using BinTable = BinTableLazy;
  template <class N>
  class iterator;

 private:
  IO_MODE _mode{IO_MODE::ReadOnly};
  HighFive::File _fp;
  RootGroup _root_group;
  GroupMap _groups;
  DatasetMap _datasets;
  StandardAttributes _attrs;
  internal::NumericVariant _pixel_variant;
  BinTable _bins;
  Index _index{};

  // Constructors are private. Cooler files are opened using factory methods
  explicit File(std::string_view uri, bool validate = true);

  template <class PixelT>
  explicit File(std::string_view uri, ChromosomeSet chroms, PixelT pixel,
                StandardAttributes attributes);

 public:
  File() = delete;
  File(const File &other) = delete;
  File(File &&other) noexcept = delete;

  // Simple constructor. Open file in read-only mode. Automatically detects pixel count type
  [[nodiscard]] static File open_read_only(std::string_view uri, bool validate = true);
  template <class PixelT = std::uint32_t>
  [[nodiscard]] static File create_new_cooler(
      std::string_view uri, const ChromosomeSet &chroms, std::uint32_t bin_size,
      bool overwrite_if_exists = false,
      StandardAttributes attributes = StandardAttributes::init<PixelT>(0));

  ~File() noexcept;

  File &operator=(const File &other) = delete;
  File &operator=(File &&other) noexcept = delete;

  // template <class PixelT, class InputIt>
  // [[nodiscard]] static  File create_new_mcool(std::string_view file_path,
  //                                                   InputIt first_resolution,
  //                                                   InputIt last_resolution,
  //                                                   bool force_overwrite = false);

  // template <class ChromSizeInputIt, class CellIDInputIt>
  // [[nodiscard]] static  File create_scool(
  //     std::string_view file_path, ChromSizeInputIt first_chrom, ChromSizeInputIt last_chrom,
  //     CellIDInputIt first_cell_id, CellIDInputIt last_cell_id, bool force_overwrite = false);
  // template <class InputIt>
  // [[nodiscard]] static  File create_scool(std::string_view file_path, InputIt first_chrom,
  //                                               InputIt last_chrom, bool force_overwrite =
  //                                               false);

  [[nodiscard]] std::string uri() const;
  [[nodiscard]] std::string hdf5_path() const;
  [[nodiscard]] std::string path() const;

  [[nodiscard]] constexpr std::uint32_t bin_size() const noexcept;
  [[nodiscard]] constexpr auto chromosomes() const noexcept -> const ChromosomeSet &;
  [[nodiscard]] constexpr auto bins() const noexcept -> const BinTable &;

  [[nodiscard]] auto attributes() const noexcept -> const StandardAttributes &;
  [[nodiscard]] auto group(std::string_view group_name) -> Group &;
  [[nodiscard]] auto dataset(std::string_view dataset_name) -> Dataset &;
  [[nodiscard]] auto group(std::string_view group_name) const -> const Group &;
  [[nodiscard]] auto dataset(std::string_view dataset_name) const -> const Dataset &;

  [[nodiscard]] const internal::NumericVariant &pixel_variant() const noexcept;
  template <class T>
  [[nodiscard]] bool has_pixel_of_type() const noexcept;

  template <class PixelIt, class = std::enable_if_t<is_iterator_v<PixelIt>>>
  void append_pixels(PixelIt first_pixel, PixelIt last_pixel, bool validate = false,
                     std::size_t chunk_size = 64 * 1024);

  template <class N>
  [[nodiscard]] auto begin() const -> iterator<N>;
  template <class N>
  [[nodiscard]] auto end() const -> iterator<N>;

  template <class N>
  [[nodiscard]] auto cbegin() const -> iterator<N>;
  template <class N>
  [[nodiscard]] auto cend() const -> iterator<N>;

  void flush();

 private:
  [[nodiscard]] static HighFive::File open_file(std::string_view uri, IO_MODE mode, bool validate);

  [[nodiscard]] static auto open_or_create_root_group(HighFive::File &f, std::string_view uri)
      -> RootGroup;

  // Open/read groups, datasets and attributes
  [[nodiscard]] static auto open_root_group(const HighFive::File &f, std::string_view uri)
      -> RootGroup;
  [[nodiscard]] static auto open_groups(const RootGroup &root_grp) -> GroupMap;
  [[nodiscard]] static auto open_datasets(const RootGroup &root_grp) -> DatasetMap;
  [[nodiscard]] static auto read_standard_attributes(const RootGroup &root_grp,
                                                     bool initialize_missing = false)
      -> StandardAttributes;

  // Create/write groups, datasets and attributes
  [[nodiscard]] static auto create_root_group(HighFive::File &f, std::string_view uri) -> RootGroup;
  [[nodiscard]] static auto create_groups(RootGroup &root_grp) -> GroupMap;
  template <class PixelT>
  [[nodiscard]] static auto create_datasets(RootGroup &root_grp, const ChromosomeSet &chroms)
      -> DatasetMap;
  static void write_standard_attributes(RootGroup &root_grp, const StandardAttributes &attributes);

  [[nodiscard]] static auto import_chroms(const Dataset &chrom_names, const Dataset &chrom_sizes,
                                          bool missing_ok) -> ChromosomeSet;

  [[nodiscard]] static Index import_indexes(const Dataset &chrom_offset_dset,
                                            const Dataset &bin_offset_dset,
                                            const ChromosomeSet &chroms,
                                            const BinTableLazy &bintable,
                                            std::uint64_t expected_nnz, bool missing_ok);

  void validate_bins() const;

  template <class PixelIt>
  void validate_pixels_before_append(PixelIt first_pixel, PixelIt last_pixel) const;

  [[nodiscard]] static internal::NumericVariant detect_pixel_type(
      const RootGroup &root_grp, std::string_view path = "pixels/count");
  void write_attributes();
  void write_chromosomes();

  template <class ChromIt, class UnaryOperation = identity,
            class = std::enable_if_t<is_unary_operation_on_iterator<UnaryOperation, ChromIt>>>
  static void write_chromosomes(Dataset &name_dset, Dataset &size_dset, ChromIt first_chrom,
                                ChromIt last_chrom, UnaryOperation op = identity());

  void write_bin_table();
  static void write_bin_table(Dataset &chrom_dset, Dataset &start_dset, Dataset &end_dset,
                              const BinTable &bin_table);
  template <class PixelIt>
  void update_indexes(PixelIt first_pixel, PixelIt last_pixel);

  void write_indexes();
  static void write_indexes(Dataset &chrom_offset_dset, Dataset &bin_offset_dset, const Index &idx);

  void finalize();

 public:
  template <class N>
  class iterator {
    static_assert(std::is_arithmetic_v<N>);
    friend File;

    const BinTableLazy *_bins{};
    Dataset::iterator<std::uint32_t, DEFAULT_HDF5_CHUNK_SIZE> _bin1_id_it{};
    Dataset::iterator<std::uint32_t, DEFAULT_HDF5_CHUNK_SIZE> _bin2_id_it{};
    Dataset::iterator<N, DEFAULT_HDF5_CHUNK_SIZE> _count_it{};
    Dataset::iterator<N, DEFAULT_HDF5_CHUNK_SIZE> _count_last{};

    explicit iterator(const File &f) noexcept;

    [[nodiscard]] static auto make_end_iterator(const File &f) -> iterator;

   public:
    using difference_type = std::ptrdiff_t;
    using value_type = Pixel<N>;
    using pointer = value_type *;
    using reference = value_type &;
    using iterator_category = std::forward_iterator_tag;

    iterator() = default;

    [[nodiscard]] bool operator==(const iterator &other) const noexcept;
    [[nodiscard]] bool operator!=(const iterator &other) const noexcept;
    [[nodiscard]] auto operator*() const -> value_type;
    auto operator++() -> iterator &;
    auto operator++(int) -> iterator;
  };
};

}  // namespace coolerpp

#include "../../coolerpp_impl.hpp"
