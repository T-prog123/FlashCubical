#pragma once

#include "zero_lookup/runtime/zero_table.hpp"

#include <array>
#include <cstdint>
#include <filesystem>
#include <vector>

namespace zero_lookup::detail {

constexpr std::uint32_t lookup_schema_version = 1;
constexpr std::uint32_t lookup_table_version = 1;
constexpr std::size_t table2d_entry_count = 256;
constexpr std::size_t table3d_es_count = 64u * 4096u;

struct LoadedTable2D {
    std::array<Entry2D, table2d_entry_count> entries{};
};

struct LoadedTable3D {
    std::vector<std::uint16_t> es;
    std::vector<std::uint16_t> full;
    std::vector<Entry3D> pairs;
};

LoadedTable2D read_table2d_file(const std::filesystem::path& path);
LoadedTable3D read_table3d_file(const std::filesystem::path& path);

} // namespace zero_lookup::detail
