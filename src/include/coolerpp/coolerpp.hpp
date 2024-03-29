// Copyright (C) 2022 Roberto Rossini <roberros@uio.no>
//
// SPDX-License-Identifier: MIT

#pragma once

#include <fmt/chrono.h>

#include <cstdint>
// clang-format off
#include "coolerpp/internal/suppress_warnings.hpp"
// clang-format on
DISABLE_WARNING_PUSH
DISABLE_WARNING_NULL_DEREF
#include <highfive/H5DataSet.hpp>
#include <highfive/H5DataSpace.hpp>
#include <highfive/H5File.hpp>
#include <highfive/H5Group.hpp>
DISABLE_WARNING_POP
#include <memory>
#include <optional>
#include <string>
#include <string_view>
#include <variant>
#include <vector>

#include "coolerpp/balancing.hpp"
#include "coolerpp/bin_table.hpp"
#include "coolerpp/chromosome.hpp"
#include "coolerpp/dataset.hpp"
#include "coolerpp/group.hpp"
#include "coolerpp/index.hpp"
#include "coolerpp/internal/numeric_variant.hpp"
#include "coolerpp/pixel.hpp"
#include "coolerpp/pixel_selector.hpp"

namespace coolerpp {

using DefaultPixelT = std::int32_t;

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
  std::optional<std::string> generated_by{COOLERPP_VERSION_STRING};
  std::optional<std::string> assembly{"unknown"};
  std::optional<std::string> metadata{"{}"};

  // Optional but common
  std::optional<std::string> format_url{"https://github.com/open2c/cooler"};
  std::optional<std::int64_t> nbins{0};
  std::optional<std::int32_t> nchroms{0};
  std::optional<std::int64_t> nnz{0};
  using SumVar = std::variant<double, std::int64_t>;
  std::optional<SumVar> sum{std::int64_t(0)};
  std::optional<SumVar> cis{std::int64_t(0)};

  template <typename PixelT = DefaultPixelT,
            typename = std::enable_if_t<std::is_arithmetic_v<PixelT>>>
  [[nodiscard]] static StandardAttributes init(std::uint32_t bin_size_);
  [[nodiscard]] static StandardAttributes init_empty() noexcept;
  [[nodiscard]] bool operator==(const StandardAttributes& other) const noexcept;
  [[nodiscard]] bool operator!=(const StandardAttributes& other) const noexcept;

 private:
  // Use the init factory methods to construct a StandardAttribute object
  StandardAttributes() = default;
};

template <typename InputIt>
void init_mcool(std::string_view file_path, InputIt first_resolution, InputIt last_resolution,
                bool force_overwrite = false);
void init_mcool(std::string_view file_path, bool force_overwrite = false);

// template <typename ChromSizeInputIt, typename CellIDInputIt>
// void init_scool(std::string_view file_path, ChromSizeInputIt first_chrom,
//                 ChromSizeInputIt last_chrom, CellIDInputIt first_cell_id,
//                 CellIDInputIt last_cell_id, bool force_overwrite = false);
// template <typename InputIt>
// void init_scool(std::string_view file_path, InputIt first_chrom, InputIt last_chrom,
//                 bool force_overwrite = false);

class File {
 public:
  enum class QUERY_TYPE { BED, UCSC };

 private:
  unsigned int _mode{HighFive::File::ReadOnly};
  std::unique_ptr<HighFive::File> _fp{};
  RootGroup _root_group{};
  GroupMap _groups{};
  DatasetMap _datasets{};
  mutable WeightMap _weights{};
  StandardAttributes _attrs{StandardAttributes::init(0)};
  internal::NumericVariant _pixel_variant{};
  std::shared_ptr<const BinTable> _bins{};
  std::shared_ptr<Index> _index{};
  bool _finalize{false};

  // Constructors are private. Cooler files are opened using factory methods
  explicit File(std::string_view uri, unsigned mode = HighFive::File::ReadOnly,
                std::size_t cache_size_bytes = DEFAULT_HDF5_CACHE_SIZE,
                double w0 = DEFAULT_HDF5_CACHE_W0, bool validate = true);

  template <typename PixelT>
  explicit File(std::string_view uri, ChromosomeSet chroms, PixelT pixel,
                StandardAttributes attributes,
                std::size_t cache_size_bytes = DEFAULT_HDF5_CACHE_SIZE,
                double w0 = DEFAULT_HDF5_CACHE_W0);

 public:
  File() = default;
  File(const File &other) = delete;
  File(File &&other) noexcept(noexcept_move_ctor()) = default;  // NOLINT

  // Simple constructor. Open file in read-only mode. Automatically detects pixel count type
  [[nodiscard]] static File open_read_only(std::string_view uri,
                                           std::size_t cache_size_bytes = DEFAULT_HDF5_CACHE_SIZE,
                                           bool validate = true);
  [[nodiscard]] static File open_read_only_random_access(
      std::string_view uri, std::size_t cache_size_bytes = DEFAULT_HDF5_CACHE_SIZE,
      bool validate = true);
  [[nodiscard]] static File open_read_only_read_once(
      std::string_view uri, std::size_t cache_size_bytes = DEFAULT_HDF5_CACHE_SIZE,
      bool validate = true);
  template <typename PixelT = DefaultPixelT>
  [[nodiscard]] static File create_new_cooler(
      std::string_view uri, const ChromosomeSet &chroms, std::uint32_t bin_size,
      bool overwrite_if_exists = false,
      StandardAttributes attributes = StandardAttributes::init<PixelT>(0),
      std::size_t cache_size_bytes = DEFAULT_HDF5_CACHE_SIZE * 4);

  ~File() noexcept;

  File &operator=(const File &other) = delete;
  File &operator=(File &&other) noexcept(noexcept_move_assigment_op()) = default;  // NOLINT

  [[nodiscard]] explicit operator bool() const noexcept;

  void open(std::string_view uri, bool validate = true);
  template <typename PixelT = DefaultPixelT>
  void create(std::string_view uri, const ChromosomeSet &chroms, std::uint32_t bin_size,
              bool overwrite_if_exists = false,
              StandardAttributes attributes = StandardAttributes::init<PixelT>(0));
  void close();

  // template <typename PixelT, typename InputIt>
  // [[nodiscard]] static  File create_new_mcool(std::string_view file_path,
  //                                                   InputIt first_resolution,
  //                                                   InputIt last_resolution,
  //                                                   bool force_overwrite = false);

  // template <typename ChromSizeInputIt, typename CellIDInputIt>
  // [[nodiscard]] static  File create_scool(
  //     std::string_view file_path, ChromSizeInputIt first_chrom, ChromSizeInputIt last_chrom,
  //     CellIDInputIt first_cell_id, CellIDInputIt last_cell_id, bool force_overwrite = false);
  // template <typename InputIt>
  // [[nodiscard]] static  File create_scool(std::string_view file_path, InputIt first_chrom,
  //                                               InputIt last_chrom, bool force_overwrite =
  //                                               false);

  [[nodiscard]] std::string uri() const;
  [[nodiscard]] std::string hdf5_path() const;
  [[nodiscard]] std::string path() const;

  [[nodiscard]] std::uint32_t bin_size() const noexcept;
  [[nodiscard]] auto chromosomes() const noexcept -> const ChromosomeSet &;
  [[nodiscard]] auto bins() const noexcept -> const BinTable &;
  [[nodiscard]] auto bins_ptr() const noexcept -> std::shared_ptr<const BinTable>;

  [[nodiscard]] auto attributes() const noexcept -> const StandardAttributes &;
  [[nodiscard]] auto group(std::string_view group_name) -> Group &;
  [[nodiscard]] auto dataset(std::string_view dataset_name) -> Dataset &;
  [[nodiscard]] auto group(std::string_view group_name) const -> const Group &;
  [[nodiscard]] auto dataset(std::string_view dataset_name) const -> const Dataset &;

  [[nodiscard]] const internal::NumericVariant &pixel_variant() const noexcept;
  template <typename T>
  [[nodiscard]] bool has_pixel_of_type() const noexcept;

  [[nodiscard]] bool has_signed_pixels() const noexcept;
  [[nodiscard]] bool has_unsigned_pixels() const noexcept;
  [[nodiscard]] bool has_integral_pixels() const noexcept;
  [[nodiscard]] bool has_float_pixels() const noexcept;

  template <typename PixelIt, typename = std::enable_if_t<is_iterable_v<PixelIt>>>
  void append_pixels(PixelIt first_pixel, PixelIt last_pixel, bool validate = false);

  template <typename N, std::size_t CHUNK_SIZE = DEFAULT_HDF5_DATASET_ITERATOR_BUFFER_SIZE>
  [[nodiscard]] typename PixelSelector<N, CHUNK_SIZE>::iterator begin() const;
  template <typename N, std::size_t CHUNK_SIZE = DEFAULT_HDF5_DATASET_ITERATOR_BUFFER_SIZE>
  [[nodiscard]] typename PixelSelector<N, CHUNK_SIZE>::iterator end() const;

  template <typename N, std::size_t CHUNK_SIZE = DEFAULT_HDF5_DATASET_ITERATOR_BUFFER_SIZE>
  [[nodiscard]] typename PixelSelector<N, CHUNK_SIZE>::iterator cbegin() const;
  template <typename N, std::size_t CHUNK_SIZE = DEFAULT_HDF5_DATASET_ITERATOR_BUFFER_SIZE>
  [[nodiscard]] typename PixelSelector<N, CHUNK_SIZE>::iterator cend() const;

  template <typename N, std::size_t CHUNK_SIZE = DEFAULT_HDF5_DATASET_ITERATOR_BUFFER_SIZE>
  [[nodiscard]] PixelSelector<N, CHUNK_SIZE> fetch(std::string_view query,
                                                   QUERY_TYPE query_type = QUERY_TYPE::UCSC) const;
  template <typename N, std::size_t CHUNK_SIZE = DEFAULT_HDF5_DATASET_ITERATOR_BUFFER_SIZE>
  [[nodiscard]] PixelSelector<N, CHUNK_SIZE> fetch(std::string_view chrom_name, std::uint32_t start,
                                                   std::uint32_t end) const;

  template <typename N, std::size_t CHUNK_SIZE = DEFAULT_HDF5_DATASET_ITERATOR_BUFFER_SIZE>
  [[nodiscard]] PixelSelector<N, CHUNK_SIZE> fetch(std::string_view range1, std::string_view range2,
                                                   QUERY_TYPE query_type = QUERY_TYPE::UCSC) const;
  template <typename N, std::size_t CHUNK_SIZE = DEFAULT_HDF5_DATASET_ITERATOR_BUFFER_SIZE>
  [[nodiscard]] PixelSelector<N, CHUNK_SIZE> fetch(std::string_view chrom1_name,
                                                   std::uint32_t start1, std::uint32_t end1,
                                                   std::string_view chrom2_name,
                                                   std::uint32_t start2, std::uint32_t end2) const;

  bool has_weights(std::string_view name) const;
  std::shared_ptr<const Weights> read_weights(std::string_view name) const;
  std::shared_ptr<const Weights> read_weights(std::string_view name, Weights::Type type) const;

  bool purge_weights(std::string_view name = "");

  void flush();

  template <typename It>
  static void write_weights(std::string_view uri, std::string_view name, It first_weight,
                            It last_weight, bool overwrite_if_exists = false,
                            bool divisive = false);
  template <typename It>
  void write_weights(std::string_view name, It first_weight, It last_weight,
                     bool overwrite_if_exists = false, bool divisive = false);

 private:
  [[nodiscard]] auto index() const noexcept -> const Index &;
  [[nodiscard]] auto index() noexcept -> Index &;

  [[nodiscard]] static HighFive::File open_file(std::string_view uri, unsigned int mode,
                                                bool validate);

  [[nodiscard]] static auto open_or_create_root_group(HighFive::File &f, std::string_view uri)
      -> RootGroup;

  // Open/read groups, datasets and attributes
  [[nodiscard]] static auto open_root_group(const HighFive::File &f, std::string_view uri)
      -> RootGroup;
  [[nodiscard]] static auto open_groups(const RootGroup &root_grp) -> GroupMap;
  [[nodiscard]] static auto open_datasets(const RootGroup &root_grp, std::size_t cache_size_bytes,
                                          double w0) -> DatasetMap;
  [[nodiscard]] static auto read_standard_attributes(const RootGroup &root_grp,
                                                     bool initialize_missing = false)
      -> StandardAttributes;

  // Create/write groups, datasets and attributes
  [[nodiscard]] static auto create_root_group(HighFive::File &f, std::string_view uri,
                                              bool write_sentinel_attr = true) -> RootGroup;
  [[nodiscard]] static auto create_groups(RootGroup &root_grp) -> GroupMap;
  template <typename PixelT>
  [[nodiscard]] static auto create_datasets(RootGroup &root_grp, const ChromosomeSet &chroms,
                                            std::size_t cache_size_bytes, double w0) -> DatasetMap;
  static void write_standard_attributes(RootGroup &root_grp, const StandardAttributes &attributes,
                                        bool skip_sentinel_attr = true);

  [[nodiscard]] static auto import_chroms(const Dataset &chrom_names, const Dataset &chrom_sizes,
                                          bool missing_ok) -> ChromosomeSet;

  [[nodiscard]] static Index import_indexes(const Dataset &chrom_offset_dset,
                                            const Dataset &bin_offset_dset,
                                            const ChromosomeSet &chroms,
                                            std::shared_ptr<const BinTable> bin_table,
                                            std::uint64_t expected_nnz, bool missing_ok);

  void validate_bins() const;

  template <typename PixelIt>
  void validate_pixels_before_append(PixelIt first_pixel, PixelIt last_pixel) const;

  [[nodiscard]] static internal::NumericVariant detect_pixel_type(
      const RootGroup &root_grp, std::string_view path = "pixels/count");
  void write_attributes(bool skip_sentinel_attr = true);
  void write_chromosomes();

  template <typename ChromIt, typename UnaryOperation = identity,
            typename = std::enable_if_t<is_unary_operation_on_iterator<UnaryOperation, ChromIt>>>
  static void write_chromosomes(Dataset &name_dset, Dataset &size_dset, ChromIt first_chrom,
                                ChromIt last_chrom, UnaryOperation op = identity());

  void write_bin_table();
  static void write_bin_table(Dataset &chrom_dset, Dataset &start_dset, Dataset &end_dset,
                              const BinTable &bin_table);
  template <typename PixelIt>
  void update_indexes(PixelIt first_pixel, PixelIt last_pixel);

  void write_indexes();
  static void write_indexes(Dataset &chrom_offset_dset, Dataset &bin_offset_dset, const Index &idx);

  void finalize();

  static void write_sentinel_attr(HighFive::Group grp);
  [[nodiscard]] static bool check_sentinel_attr(const HighFive::Group &grp);
  void write_sentinel_attr();
  [[nodiscard]] bool check_sentinel_attr();

  [[nodiscard]] Bin get_last_bin_written() const;

  template <typename N, bool cis = false>
  void update_pixel_sum(N partial_sum);

  template <typename PixelT>
  void validate_pixel_type() const noexcept;

  // IMPORTANT: the private fetch() methods interpret queries as open-open
  template <typename N, std::size_t CHUNK_SIZE>
  [[nodiscard]] PixelSelector<N, CHUNK_SIZE> fetch(PixelCoordinates coord) const;
  template <typename N, std::size_t CHUNK_SIZE>
  [[nodiscard]] PixelSelector<N, CHUNK_SIZE> fetch(PixelCoordinates coord1,
                                                   PixelCoordinates coord2) const;
};

}  // namespace coolerpp

#include "../../coolerpp_accessors_impl.hpp"
#include "../../coolerpp_impl.hpp"
#include "../../coolerpp_read_impl.hpp"
#include "../../coolerpp_standard_attr_impl.hpp"
#include "../../coolerpp_validation_impl.hpp"
#include "../../coolerpp_write_impl.hpp"
