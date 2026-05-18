#pragma once

#include <array>
#include <cstddef>
#include <cstdint>
#include <vector>

namespace zero_lookup {

enum Edge2DBit : std::uint8_t {
    edge2_n = 0,
    edge2_e = 1,
    edge2_s = 2,
    edge2_w = 3
};

enum Square2DBit : std::uint8_t {
    square2_ne = 0,
    square2_se = 1,
    square2_sw = 2,
    square2_nw = 3
};

struct Direct2DPair {
    std::uint8_t edge = 0xffu;
    std::uint8_t square = 0xffu;
};

struct Entry2D {
    std::uint8_t edge_mask = 0;
    std::uint8_t square_mask = 0;
    std::uint8_t edge_h0_negative_mask = 0;
    std::uint8_t edge_h1_positive_mask = 0;
    std::uint8_t square_h1_negative_mask = 0;
    std::uint8_t survivor_edge_mask = 0;
    std::uint8_t survivor_square_mask = 0;
    std::uint8_t direct_h0_edge = 0xffu;
    std::uint8_t direct_h1_dual_count = 0;
    std::uint8_t edge_h1_positive_count = 0;
    std::uint8_t square_h1_negative_count = 0;
    std::uint8_t survivor_edge_count = 0;
    std::uint8_t survivor_square_count = 0;
    std::array<Direct2DPair, 4> direct_h1_dual_pairs{};
    std::array<std::uint8_t, 4> edge_second_rank{};
    std::array<std::uint8_t, 4> square_second_rank{};
    std::array<std::uint8_t, 4> square_birth_local_edge{};
    std::array<std::uint8_t, 4> edge_order_asc{};
    std::array<std::uint8_t, 4> edge_order_desc{};
    std::array<std::uint8_t, 4> square_order_asc{};
    std::array<std::uint8_t, 4> square_order_desc{};
    std::uint8_t residual_edge_order_asc_count = 0;
    std::uint8_t residual_edge_order_desc_count = 0;
    std::uint8_t residual_square_order_asc_count = 0;
    std::array<std::uint8_t, 4> residual_edge_order_asc{};
    std::array<std::uint8_t, 4> residual_edge_order_desc{};
    std::array<std::uint8_t, 4> residual_square_order_asc{};
};

enum Edge3DBit : std::uint8_t {
    edge3_x_neg = 0,
    edge3_x_pos = 1,
    edge3_y_neg = 2,
    edge3_y_pos = 3,
    edge3_z_neg = 4,
    edge3_z_pos = 5
};

struct Direct3DH1Pair {
    std::uint8_t edge = 0xffu;
    std::uint8_t square = 0xffu;
};

struct Direct3DH2Pair {
    std::uint8_t square = 0xffu;
    std::uint8_t cube = 0xffu;
};

struct Entry3D {
    std::uint8_t edge_mask = 0;
    std::uint16_t square_mask = 0;
    std::uint8_t cube_mask = 0;

    std::uint8_t edge_h0_negative_mask = 0;
    std::uint8_t edge_h1_positive_mask = 0;
    std::uint16_t square_h1_negative_mask = 0;
    std::uint16_t square_h2_negative_mask = 0;
    std::uint16_t square_h2_skip_mask = 0;
    std::uint8_t survivor_edge_mask = 0;
    std::uint16_t survivor_square_mask = 0;
    std::uint8_t survivor_cube_mask = 0;
    std::uint16_t h2_residual_square_mask = 0;

    std::uint8_t direct_h0_edge = 0xffu;
    std::uint8_t zero_h1_count = 0;
    std::uint8_t direct_h2_count = 0;
    std::array<Direct3DH1Pair, 12> zero_h1_pairs{};
    std::array<std::uint8_t, 12> h1_square_semireduced_boundary{};
    std::array<std::uint16_t, 12> h1_square_substitution_square_mask{};
    std::array<std::uint8_t, 12> zero_h1_rewrite_masks{};
    std::array<Direct3DH2Pair, 8> direct_h2_pairs{};

    std::uint8_t edge_count = 0;
    std::uint8_t square_count = 0;
    std::uint8_t cube_count = 0;
    std::array<std::uint8_t, 6>  edge_second_rank{};
    std::array<std::uint8_t, 12> square_second_rank{};
    std::array<std::uint8_t, 8>  cube_second_rank{};
    std::array<std::uint8_t, 12> square_birth_local_edge{};
    std::array<std::uint8_t, 8>  cube_birth_local_square{};
    std::array<std::uint8_t, 6>  edge_order_asc{};
    std::array<std::uint8_t, 6>  edge_order_desc{};
    std::array<std::uint8_t, 12> square_order_asc{};
    std::array<std::uint8_t, 12> square_order_desc{};
    std::array<std::uint8_t, 8>  cube_order_asc{};
    std::array<std::uint8_t, 8>  cube_order_desc{};
    std::uint8_t residual_edge_asc_count = 0;
    std::uint8_t residual_edge_desc_count = 0;
    std::uint8_t residual_square_asc_count = 0;
    std::uint8_t residual_square_desc_count = 0;
    std::uint8_t h2_residual_square_asc_count = 0;
    std::uint8_t residual_cube_asc_count = 0;
    std::array<std::uint8_t, 6>  residual_edge_order_asc{};
    std::array<std::uint8_t, 6>  residual_edge_order_desc{};
    std::array<std::uint8_t, 12> residual_square_order_asc{};
    std::array<std::uint8_t, 12> residual_square_order_desc{};
    std::array<std::uint8_t, 12> h2_residual_square_order_asc{};
    std::array<std::uint8_t, 8>  residual_cube_order_asc{};
};

const Entry2D& lookup2d(std::uint8_t mask);
const Entry3D& lookup3d(std::uint32_t mask26);
std::size_t    entry3d_count();
const Entry3D& entry3d_at(std::size_t idx);

inline std::uint8_t direct_h0_edge2d(std::uint8_t edge_mask) {
    if ((edge_mask & (1u << edge2_n)) != 0) {
        return edge2_n;
    }
    if ((edge_mask & (1u << edge2_w)) != 0) {
        return edge2_w;
    }
    if ((edge_mask & (1u << edge2_e)) != 0) {
        return edge2_e;
    }
    if ((edge_mask & (1u << edge2_s)) != 0) {
        return edge2_s;
    }
    return 0xffu;
}

inline std::uint8_t direct_h0_edge3d(std::uint8_t edge_mask) {
    if ((edge_mask & (1u << edge3_z_neg)) != 0) {
        return edge3_z_neg;
    }
    if ((edge_mask & (1u << edge3_y_neg)) != 0) {
        return edge3_y_neg;
    }
    if ((edge_mask & (1u << edge3_x_neg)) != 0) {
        return edge3_x_neg;
    }
    if ((edge_mask & (1u << edge3_x_pos)) != 0) {
        return edge3_x_pos;
    }
    if ((edge_mask & (1u << edge3_y_pos)) != 0) {
        return edge3_y_pos;
    }
    if ((edge_mask & (1u << edge3_z_pos)) != 0) {
        return edge3_z_pos;
    }
    return 0xffu;
}

inline std::uint8_t local_cube_mask3d(std::uint32_t mask26) {
    std::uint8_t cubes = 0;
    for (std::uint8_t cube = 0; cube < 8; ++cube) {
        const bool xp = (cube & 1u) != 0;
        const bool yp = (cube & 2u) != 0;
        const bool zp = (cube & 4u) != 0;
        const std::uint32_t ex = xp ? edge3_x_pos : edge3_x_neg;
        const std::uint32_t ey = yp ? edge3_y_pos : edge3_y_neg;
        const std::uint32_t ez = zp ? edge3_z_pos : edge3_z_neg;
        const std::uint32_t xy = 6u + (xp ? 1u : 0u) + 2u * (yp ? 1u : 0u);
        const std::uint32_t xz = 10u + (xp ? 1u : 0u) + 2u * (zp ? 1u : 0u);
        const std::uint32_t yz = 14u + (yp ? 1u : 0u) + 2u * (zp ? 1u : 0u);
        const std::uint32_t body =
            18u + (xp ? 1u : 0u) + 2u * (yp ? 1u : 0u) +
            4u * (zp ? 1u : 0u);
        const std::uint32_t required =
            (1u << ex) | (1u << ey) | (1u << ez) |
            (1u << xy) | (1u << xz) | (1u << yz) | (1u << body);
        if ((mask26 & required) == required) {
            cubes = static_cast<std::uint8_t>(cubes | (1u << cube));
        }
    }
    return cubes;
}

struct Table3DStats {
    std::size_t edge_square_states = 0;
    std::size_t geometry_states = 0;
    std::size_t es_bytes = 0;
    std::size_t full_bytes = 0;
    std::size_t pair_bytes = 0;

    std::size_t total_bytes() const {
        return es_bytes + full_bytes + pair_bytes;
    }
};

struct RunStats {
    double mask_build_ms = 0.0;
    double table_apply_ms = 0.0;

    std::size_t direct_h0_pairs_found = 0;
    std::size_t direct_h0_pairs_applied = 0;
    std::size_t direct_h1_dual_pairs_found = 0;
    std::size_t direct_h1_dual_pairs_applied = 0;
    std::size_t direct_h2_dual_pairs_found = 0;
    std::size_t direct_h2_dual_pairs_applied = 0;

    std::size_t lookup_h1_positive_edges = 0;
    std::size_t lookup_h1_positive_edges_applied = 0;
    std::size_t lookup_h1_negative_squares = 0;
    std::size_t lookup_h2_negative_squares = 0;
    std::size_t lookup_survivor_edges = 0;
    std::size_t lookup_survivor_squares = 0;
    std::size_t lookup_survivor_cubes = 0;
};

void warm_zero_table2d();
void warm_zero_table3d();
void warm_zero_tables();
Table3DStats table3d_stats();

}
