#include "zero_lookup/runtime/zero_table.hpp"
#include "zero_lookup/runtime/zero_table_file.hpp"

#include <filesystem>
#include <memory>
#include <mutex>
#include <stdexcept>

namespace zero_lookup {
namespace {

std::mutex& load_mutex() {
    static std::mutex m;
    return m;
}

std::unique_ptr<detail::LoadedTable2D>& loaded_table2d() {
    static std::unique_ptr<detail::LoadedTable2D> table;
    return table;
}

std::unique_ptr<detail::LoadedTable3D>& loaded_table3d() {
    static std::unique_ptr<detail::LoadedTable3D> table;
    return table;
}

const detail::LoadedTable2D& require_table2d() {
    const std::unique_ptr<detail::LoadedTable2D>& table = loaded_table2d();
    if (!table) {
        throw std::runtime_error(
            "flash_cubical: 2D zero lookup table is not loaded");
    }
    return *table;
}

const detail::LoadedTable3D& require_table3d() {
    const std::unique_ptr<detail::LoadedTable3D>& table = loaded_table3d();
    if (!table) {
        throw std::runtime_error(
            "flash_cubical: 3D zero lookup table is not loaded");
    }
    return *table;
}

} // namespace

const Entry2D& lookup2d(std::uint8_t mask) {
    return require_table2d().entries[mask];
}

const Entry3D& lookup3d(std::uint32_t mask26) {
    const detail::LoadedTable3D& table = require_table3d();
    const std::uint32_t a = mask26 & 0x3fu;
    const std::uint32_t f = (mask26 >> 6) & 0xfffu;
    const std::uint32_t b = (mask26 >> 18) & 0xffu;
    const std::uint16_t es = table.es[a * 4096u + f];
    const std::uint16_t geom =
        table.full[static_cast<std::size_t>(es) * 256u + b];
    return table.pairs[geom];
}

void warm_zero_tables() {
    warm_zero_table2d();
    warm_zero_table3d();
}

void warm_zero_table2d() {
    (void)require_table2d();
}

void warm_zero_table3d() {
    (void)require_table3d();
}

Table3DStats table3d_stats() {
    const detail::LoadedTable3D& table = require_table3d();
    return {table.full.size() / 256u,
            table.pairs.size(),
            table.es.size() * sizeof(std::uint16_t),
            table.full.size() * sizeof(std::uint16_t),
            table.pairs.size() * sizeof(Entry3D)};
}

std::size_t entry3d_count() {
    return require_table3d().pairs.size();
}

const Entry3D& entry3d_at(std::size_t idx) {
    const detail::LoadedTable3D& table = require_table3d();
    if (idx >= table.pairs.size()) {
        throw std::out_of_range(
            "flash_cubical: 3D lookup entry index out of range");
    }
    return table.pairs[idx];
}

void load_zero_table2d(const std::filesystem::path& path) {
    auto loaded = std::make_unique<detail::LoadedTable2D>(
        detail::read_table2d_file(path));
    std::lock_guard<std::mutex> lock(load_mutex());
    loaded_table2d() = std::move(loaded);
}

void load_zero_table3d(const std::filesystem::path& path) {
    auto loaded = std::make_unique<detail::LoadedTable3D>(
        detail::read_table3d_file(path));
    std::lock_guard<std::mutex> lock(load_mutex());
    loaded_table3d() = std::move(loaded);
}

void load_zero_tables(const std::filesystem::path& path2d,
                      const std::filesystem::path& path3d) {
    auto loaded2d = std::make_unique<detail::LoadedTable2D>(
        detail::read_table2d_file(path2d));
    auto loaded3d = std::make_unique<detail::LoadedTable3D>(
        detail::read_table3d_file(path3d));
    std::lock_guard<std::mutex> lock(load_mutex());
    loaded_table2d() = std::move(loaded2d);
    loaded_table3d() = std::move(loaded3d);
}

void load_zero_tables_from_directory(const std::filesystem::path& dir) {
    load_zero_tables(dir / "zero_lookup_2d.bin",
                     dir / "zero_lookup_3d.bin");
}

bool zero_table2d_loaded() {
    return static_cast<bool>(loaded_table2d());
}

bool zero_table3d_loaded() {
    return static_cast<bool>(loaded_table3d());
}

} // namespace zero_lookup
