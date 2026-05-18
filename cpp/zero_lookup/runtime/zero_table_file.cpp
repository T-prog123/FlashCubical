#include "zero_lookup/runtime/zero_table_file.hpp"

#include <array>
#include <cstddef>
#include <fstream>
#include <stdexcept>
#include <string>

namespace zero_lookup::detail {
namespace {

static_assert(sizeof(std::uint8_t) == 1, "uint8_t must be one byte");
static_assert(sizeof(std::uint16_t) == 2, "uint16_t must be two bytes");
static_assert(sizeof(std::uint32_t) == 4, "uint32_t must be four bytes");
static_assert(sizeof(std::uint64_t) == 8, "uint64_t must be eight bytes");

constexpr std::array<char, 8> magic2d{{'F', 'C', 'Z', 'L', 'U', 'T', '2', 'D'}};
constexpr std::array<char, 8> magic3d{{'F', 'C', 'Z', 'L', 'U', 'T', '3', 'D'}};
constexpr std::uint64_t fnv_offset = 14695981039346656037ull;
constexpr std::uint64_t fnv_prime = 1099511628211ull;

std::runtime_error file_error(const std::filesystem::path& path,
                              const std::string& message) {
    return std::runtime_error("flash_cubical: lookup table file '" +
                              path.string() + "': " + message);
}

class Reader {
public:
    Reader(const std::filesystem::path& path,
           const std::vector<std::uint8_t>& data)
        : path_(path), data_(data) {}

    std::uint8_t u8() {
        require(1);
        return data_[pos_++];
    }

    std::uint16_t u16() {
        require(2);
        const std::uint16_t v =
            static_cast<std::uint16_t>(data_[pos_]) |
            static_cast<std::uint16_t>(data_[pos_ + 1] << 8u);
        pos_ += 2;
        return v;
    }

    std::uint32_t u32() {
        require(4);
        std::uint32_t v = 0;
        for (unsigned i = 0; i < 4; ++i) {
            v |= static_cast<std::uint32_t>(data_[pos_ + i]) << (8u * i);
        }
        pos_ += 4;
        return v;
    }

    std::uint64_t u64() {
        require(8);
        std::uint64_t v = 0;
        for (unsigned i = 0; i < 8; ++i) {
            v |= static_cast<std::uint64_t>(data_[pos_ + i]) << (8u * i);
        }
        pos_ += 8;
        return v;
    }

    void done() const {
        if (pos_ != data_.size()) {
            throw file_error(path_, "unexpected trailing payload bytes");
        }
    }

private:
    const std::filesystem::path& path_;
    const std::vector<std::uint8_t>& data_;
    std::size_t pos_ = 0;

    void require(std::size_t n) const {
        if (n > data_.size() - pos_) {
            throw file_error(path_, "truncated payload");
        }
    }
};

std::uint64_t checksum(const std::vector<std::uint8_t>& payload) {
    std::uint64_t h = fnv_offset;
    for (std::uint8_t b : payload) {
        h ^= b;
        h *= fnv_prime;
    }
    return h;
}

template <std::size_t N>
void read_u8_array(Reader& r, std::array<std::uint8_t, N>& a) {
    for (std::uint8_t& v : a) {
        v = r.u8();
    }
}

template <std::size_t N>
void read_u16_array(Reader& r, std::array<std::uint16_t, N>& a) {
    for (std::uint16_t& v : a) {
        v = r.u16();
    }
}

Direct2DPair read_direct2d_pair(Reader& r) {
    Direct2DPair p;
    p.edge = r.u8();
    p.square = r.u8();
    return p;
}

Direct3DH1Pair read_direct3d_h1_pair(Reader& r) {
    Direct3DH1Pair p;
    p.edge = r.u8();
    p.square = r.u8();
    return p;
}

Direct3DH2Pair read_direct3d_h2_pair(Reader& r) {
    Direct3DH2Pair p;
    p.square = r.u8();
    p.cube = r.u8();
    return p;
}

Entry2D read_entry2d(Reader& r) {
    Entry2D e;
    e.edge_mask = r.u8();
    e.square_mask = r.u8();
    e.edge_h0_negative_mask = r.u8();
    e.edge_h1_positive_mask = r.u8();
    e.square_h1_negative_mask = r.u8();
    e.survivor_edge_mask = r.u8();
    e.survivor_square_mask = r.u8();
    e.direct_h0_edge = r.u8();
    e.direct_h1_dual_count = r.u8();
    e.edge_h1_positive_count = r.u8();
    e.square_h1_negative_count = r.u8();
    e.survivor_edge_count = r.u8();
    e.survivor_square_count = r.u8();
    for (Direct2DPair& p : e.direct_h1_dual_pairs) {
        p = read_direct2d_pair(r);
    }
    read_u8_array(r, e.edge_second_rank);
    read_u8_array(r, e.square_second_rank);
    read_u8_array(r, e.square_birth_local_edge);
    read_u8_array(r, e.edge_order_asc);
    read_u8_array(r, e.edge_order_desc);
    read_u8_array(r, e.square_order_asc);
    read_u8_array(r, e.square_order_desc);
    e.residual_edge_order_asc_count = r.u8();
    e.residual_edge_order_desc_count = r.u8();
    e.residual_square_order_asc_count = r.u8();
    read_u8_array(r, e.residual_edge_order_asc);
    read_u8_array(r, e.residual_edge_order_desc);
    read_u8_array(r, e.residual_square_order_asc);
    return e;
}

Entry3D read_entry3d(Reader& r) {
    Entry3D e;
    e.edge_mask = r.u8();
    e.square_mask = r.u16();
    e.cube_mask = r.u8();
    e.edge_h0_negative_mask = r.u8();
    e.edge_h1_positive_mask = r.u8();
    e.square_h1_negative_mask = r.u16();
    e.square_h2_negative_mask = r.u16();
    e.square_h2_skip_mask = r.u16();
    e.survivor_edge_mask = r.u8();
    e.survivor_square_mask = r.u16();
    e.survivor_cube_mask = r.u8();
    e.h2_residual_square_mask = r.u16();
    e.direct_h0_edge = r.u8();
    e.zero_h1_count = r.u8();
    e.direct_h2_count = r.u8();
    for (Direct3DH1Pair& p : e.zero_h1_pairs) {
        p = read_direct3d_h1_pair(r);
    }
    read_u8_array(r, e.h1_square_semireduced_boundary);
    read_u16_array(r, e.h1_square_substitution_square_mask);
    read_u8_array(r, e.zero_h1_rewrite_masks);
    for (Direct3DH2Pair& p : e.direct_h2_pairs) {
        p = read_direct3d_h2_pair(r);
    }
    e.edge_count = r.u8();
    e.square_count = r.u8();
    e.cube_count = r.u8();
    read_u8_array(r, e.edge_second_rank);
    read_u8_array(r, e.square_second_rank);
    read_u8_array(r, e.cube_second_rank);
    read_u8_array(r, e.square_birth_local_edge);
    read_u8_array(r, e.cube_birth_local_square);
    read_u8_array(r, e.edge_order_asc);
    read_u8_array(r, e.edge_order_desc);
    read_u8_array(r, e.square_order_asc);
    read_u8_array(r, e.square_order_desc);
    read_u8_array(r, e.cube_order_asc);
    read_u8_array(r, e.cube_order_desc);
    e.residual_edge_asc_count = r.u8();
    e.residual_edge_desc_count = r.u8();
    e.residual_square_asc_count = r.u8();
    e.residual_square_desc_count = r.u8();
    e.h2_residual_square_asc_count = r.u8();
    e.residual_cube_asc_count = r.u8();
    read_u8_array(r, e.residual_edge_order_asc);
    read_u8_array(r, e.residual_edge_order_desc);
    read_u8_array(r, e.residual_square_order_asc);
    read_u8_array(r, e.residual_square_order_desc);
    read_u8_array(r, e.h2_residual_square_order_asc);
    read_u8_array(r, e.residual_cube_order_asc);
    return e;
}

std::vector<std::uint8_t> read_payload(const std::filesystem::path& path,
                                       const std::array<char, 8>& magic) {
    std::ifstream in(path, std::ios::binary);
    if (!in) {
        throw file_error(path, "could not open file");
    }

    std::array<char, 8> got_magic{};
    in.read(got_magic.data(), static_cast<std::streamsize>(got_magic.size()));
    if (!in || got_magic != magic) {
        throw file_error(path, "bad magic");
    }

    std::array<std::uint8_t, 24> header{};
    in.read(reinterpret_cast<char*>(header.data()),
            static_cast<std::streamsize>(header.size()));
    if (!in) {
        throw file_error(path, "truncated header");
    }

    const std::vector<std::uint8_t> header_payload(header.begin(),
                                                   header.end());
    Reader hr(path, header_payload);
    const std::uint32_t schema = hr.u32();
    const std::uint32_t table_version = hr.u32();
    const std::uint64_t payload_size = hr.u64();
    const std::uint64_t expected_checksum = hr.u64();
    if (schema != lookup_schema_version) {
        throw file_error(path, "schema version mismatch");
    }
    if (table_version != lookup_table_version) {
        throw file_error(path, "table version mismatch");
    }
    if (payload_size > static_cast<std::uint64_t>(1) << 34) {
        throw file_error(path, "payload size is unreasonable");
    }

    std::vector<std::uint8_t> payload(static_cast<std::size_t>(payload_size));
    if (!payload.empty()) {
        in.read(reinterpret_cast<char*>(payload.data()),
                static_cast<std::streamsize>(payload.size()));
    }
    if (!in) {
        throw file_error(path, "truncated payload");
    }
    char extra = 0;
    if (in.read(&extra, 1)) {
        throw file_error(path, "unexpected trailing file bytes");
    }
    if (checksum(payload) != expected_checksum) {
        throw file_error(path, "checksum mismatch");
    }
    return payload;
}

} // namespace

LoadedTable2D read_table2d_file(const std::filesystem::path& path) {
    const std::vector<std::uint8_t> payload = read_payload(path, magic2d);
    Reader r(path, payload);
    const std::uint32_t count = r.u32();
    if (count != table2d_entry_count) {
        throw file_error(path, "unexpected 2D entry count");
    }
    LoadedTable2D out;
    for (Entry2D& e : out.entries) {
        e = read_entry2d(r);
    }
    r.done();
    return out;
}

LoadedTable3D read_table3d_file(const std::filesystem::path& path) {
    const std::vector<std::uint8_t> payload = read_payload(path, magic3d);
    Reader r(path, payload);
    const std::uint32_t es_count = r.u32();
    const std::uint32_t full_count = r.u32();
    const std::uint32_t pair_count = r.u32();
    if (es_count != table3d_es_count) {
        throw file_error(path, "unexpected 3D ES table size");
    }
    if (pair_count == 0 || pair_count > 0xffffu) {
        throw file_error(path, "unexpected 3D pair table size");
    }
    if (full_count == 0 || full_count % 256u != 0) {
        throw file_error(path, "unexpected 3D FULL table size");
    }

    LoadedTable3D out;
    out.es.resize(es_count);
    out.full.resize(full_count);
    out.pairs.resize(pair_count);
    for (std::uint16_t& v : out.es) {
        v = r.u16();
        if (static_cast<std::uint32_t>(v) * 256u >= full_count) {
            throw file_error(path, "ES index out of range");
        }
    }
    for (std::uint16_t& v : out.full) {
        v = r.u16();
        if (v >= pair_count) {
            throw file_error(path, "FULL index out of range");
        }
    }
    for (Entry3D& e : out.pairs) {
        e = read_entry3d(r);
    }
    r.done();
    return out;
}

} // namespace zero_lookup::detail
