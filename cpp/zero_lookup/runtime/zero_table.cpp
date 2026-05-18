#include "zero_lookup/runtime/zero_table.hpp"

#include <algorithm>
#include <cstddef>
#include <stdexcept>
#include <unordered_map>

namespace zero_lookup {
namespace {

constexpr std::uint8_t no_cell = 0xffu;

std::uint8_t bit(std::uint8_t i) {
    return static_cast<std::uint8_t>(1u << i);
}

std::uint8_t popcount4(std::uint8_t x) {
    x &= 0x0fu;
    return static_cast<std::uint8_t>((x & 1u) + ((x >> 1) & 1u) +
                                     ((x >> 2) & 1u) + ((x >> 3) & 1u));
}

std::uint16_t bit16(std::uint8_t i) {
    return static_cast<std::uint16_t>(1u << i);
}

struct EdgeColumn {
    std::uint8_t bits = 0;
};

struct SquareColumn {
    std::uint16_t bits = 0;
};

std::uint8_t square_boundary_edges_2d(std::uint8_t s) {
    switch (s) {
    case square2_ne:
        return bit(edge2_n) | bit(edge2_e);
    case square2_se:
        return bit(edge2_s) | bit(edge2_e);
    case square2_sw:
        return bit(edge2_s) | bit(edge2_w);
    case square2_nw:
        return bit(edge2_n) | bit(edge2_w);
    default:
        return 0;
    }
}

constexpr std::array<std::uint8_t, 4> rule_c_edge_order{{
    edge2_n,
    edge2_e,
    edge2_s,
    edge2_w
}};

std::uint8_t rule_c_square_label(std::uint8_t s) {
    switch (s) {
    case square2_nw:
        return 0;
    case square2_ne:
        return 1;
    case square2_se:
        return 2;
    case square2_sw:
        return 3;
    default:
        return no_cell;
    }
}

std::uint8_t square_from_rule_c_label(std::uint8_t label) {
    switch (label) {
    case 0:
        return square2_nw;
    case 1:
        return square2_ne;
    case 2:
        return square2_se;
    case 3:
        return square2_sw;
    default:
        return no_cell;
    }
}

std::uint8_t smallest_rule_c_square(std::uint8_t mask) {
    for (std::uint8_t label = 0; label < 4; ++label) {
        const std::uint8_t s = square_from_rule_c_label(label);
        if ((mask & bit(s)) != 0) {
            return s;
        }
    }
    return no_cell;
}

std::uint8_t square_adjacency_2d(std::uint8_t s, std::uint8_t square_mask) {
    std::uint8_t adjacent = 0;
    const std::uint8_t boundary = square_boundary_edges_2d(s);
    for (std::uint8_t t = 0; t < 4; ++t) {
        if (t == s || (square_mask & bit(t)) == 0) {
            continue;
        }
        if ((boundary & square_boundary_edges_2d(t)) != 0) {
            adjacent = static_cast<std::uint8_t>(adjacent | bit(t));
        }
    }
    return adjacent;
}

std::uint8_t square_degree_2d(std::uint8_t s, std::uint8_t square_mask) {
    return popcount4(square_adjacency_2d(s, square_mask));
}

std::uint8_t square_component_2d(std::uint8_t start,
                                 std::uint8_t square_mask) {
    std::uint8_t component = 0;
    std::uint8_t frontier = bit(start);
    while (frontier != 0) {
        const std::uint8_t s = smallest_rule_c_square(frontier);
        frontier = static_cast<std::uint8_t>(frontier & ~bit(s));
        if ((component & bit(s)) != 0) {
            continue;
        }
        component = static_cast<std::uint8_t>(component | bit(s));
        frontier = static_cast<std::uint8_t>(
            frontier | (square_adjacency_2d(s, square_mask) & ~component));
    }
    return component;
}

std::uint8_t rule_c_component_start(std::uint8_t component) {
    std::uint8_t degree_one = 0;
    for (std::uint8_t s = 0; s < 4; ++s) {
        if ((component & bit(s)) != 0 &&
            square_degree_2d(s, component) == 1u) {
            degree_one = static_cast<std::uint8_t>(degree_one | bit(s));
        }
    }
    return smallest_rule_c_square(degree_one != 0 ? degree_one : component);
}

std::array<std::uint8_t, 4> rule_c_square_visit_order(
    std::uint8_t square_mask,
    std::uint8_t& count) {
    std::array<std::uint8_t, 4> order{};
    order.fill(no_cell);
    count = 0;

    std::uint8_t remaining = square_mask;
    while (remaining != 0) {
        const std::uint8_t seed = smallest_rule_c_square(remaining);
        const std::uint8_t component =
            square_component_2d(seed, remaining);
        std::uint8_t visited = 0;
        std::uint8_t current = rule_c_component_start(component);

        while (visited != component) {
            if ((visited & bit(current)) == 0) {
                order[count++] = current;
                visited = static_cast<std::uint8_t>(visited | bit(current));
                remaining = static_cast<std::uint8_t>(remaining & ~bit(current));
            }

            const std::uint8_t unvisited_adjacent =
                static_cast<std::uint8_t>(
                    square_adjacency_2d(current, component) & ~visited);
            if (unvisited_adjacent != 0) {
                current = smallest_rule_c_square(unvisited_adjacent);
                continue;
            }

            const std::uint8_t unvisited =
                static_cast<std::uint8_t>(component & ~visited);
            if (unvisited == 0) {
                break;
            }

            std::uint8_t touching = 0;
            for (std::uint8_t s = 0; s < 4; ++s) {
                if ((unvisited & bit(s)) != 0 &&
                    (square_adjacency_2d(s, component) & visited) != 0) {
                    touching = static_cast<std::uint8_t>(touching | bit(s));
                }
            }
            if (touching == 0) {
                current = smallest_rule_c_square(unvisited);
                continue;
            }

            std::uint8_t degree_one = 0;
            for (std::uint8_t s = 0; s < 4; ++s) {
                if ((touching & bit(s)) != 0 &&
                    square_degree_2d(s, component) == 1u) {
                    degree_one = static_cast<std::uint8_t>(degree_one | bit(s));
                }
            }
            current = smallest_rule_c_square(
                degree_one != 0 ? degree_one : touching);
        }
    }

    return order;
}

std::uint8_t youngest_edge_2d(std::uint8_t mask,
                              const std::array<std::uint8_t, 4>& rank) {
    std::uint8_t best = no_cell;
    for (std::uint8_t e = 0; e < 4; ++e) {
        if ((mask & bit(e)) == 0) {
            continue;
        }
        if (best == no_cell || rank[best] < rank[e]) {
            best = e;
        }
    }
    return best;
}

void fill_rule_c_metadata_2d(Entry2D& out) {
    out.edge_second_rank.fill(no_cell);
    out.square_second_rank.fill(no_cell);
    out.square_birth_local_edge.fill(no_cell);
    out.edge_order_asc.fill(no_cell);
    out.edge_order_desc.fill(no_cell);
    out.square_order_asc.fill(no_cell);
    out.square_order_desc.fill(no_cell);
    out.residual_edge_order_asc.fill(no_cell);
    out.residual_edge_order_desc.fill(no_cell);
    out.residual_square_order_asc.fill(no_cell);
    out.residual_edge_order_asc_count = 0;
    out.residual_edge_order_desc_count = 0;
    out.residual_square_order_asc_count = 0;

    std::array<bool, 4> edge_added{{false, false, false, false}};
    std::uint8_t edge_count = 0;
    auto add_edge = [&](std::uint8_t e) {
        if (e == no_cell || edge_added[e]) {
            return;
        }
        edge_added[e] = true;
        out.edge_second_rank[e] = edge_count;
        out.edge_order_asc[edge_count++] = e;
    };

    std::uint8_t square_count = 0;
    auto add_square = [&](std::uint8_t s) {
        out.square_second_rank[s] = square_count;
        out.square_order_asc[square_count++] = s;
    };

    std::uint8_t target_square_count = 0;
    const std::array<std::uint8_t, 4> target_squares =
        rule_c_square_visit_order(out.square_mask, target_square_count);
    for (std::uint8_t i = 0; i < target_square_count; ++i) {
        const std::uint8_t s = target_squares[i];
        const std::uint8_t boundary =
            static_cast<std::uint8_t>(
                square_boundary_edges_2d(s) & out.edge_mask);
        for (std::uint8_t e : rule_c_edge_order) {
            if ((boundary & bit(e)) != 0) {
                add_edge(e);
            }
        }
        add_square(s);
    }

    for (std::uint8_t e : rule_c_edge_order) {
        if ((out.edge_mask & bit(e)) != 0) {
            add_edge(e);
        }
    }

    for (std::uint8_t i = 0; i < edge_count; ++i) {
        out.edge_order_desc[i] = out.edge_order_asc[edge_count - 1u - i];
    }
    for (std::uint8_t i = 0; i < square_count; ++i) {
        out.square_order_desc[i] =
            out.square_order_asc[square_count - 1u - i];
    }

    for (std::uint8_t s = 0; s < 4; ++s) {
        if ((out.square_mask & bit(s)) == 0) {
            continue;
        }
        const std::uint8_t boundary =
            static_cast<std::uint8_t>(
                square_boundary_edges_2d(s) & out.edge_mask);
        std::uint8_t best = no_cell;
        for (std::uint8_t e = 0; e < 4; ++e) {
            if ((boundary & bit(e)) == 0) {
                continue;
            }
            if (best == no_cell ||
                out.edge_second_rank[best] < out.edge_second_rank[e]) {
                best = e;
            }
        }
        out.square_birth_local_edge[s] = best;
    }
}

Entry2D build_entry2d(std::uint8_t raw) {
    Entry2D out{};
    out.edge_mask = static_cast<std::uint8_t>(raw & 0x0fu);
    if ((raw & (bit(0) | bit(1) | bit(4))) == (bit(0) | bit(1) | bit(4))) {
        out.square_mask |= bit(square2_ne);
    }
    if ((raw & (bit(2) | bit(1) | bit(5))) == (bit(2) | bit(1) | bit(5))) {
        out.square_mask |= bit(square2_se);
    }
    if ((raw & (bit(2) | bit(3) | bit(6))) == (bit(2) | bit(3) | bit(6))) {
        out.square_mask |= bit(square2_sw);
    }
    if ((raw & (bit(0) | bit(3) | bit(7))) == (bit(0) | bit(3) | bit(7))) {
        out.square_mask |= bit(square2_nw);
    }
    fill_rule_c_metadata_2d(out);

    const std::uint8_t h0_edge = out.edge_order_asc[0];
    out.direct_h0_edge = h0_edge;
    if (h0_edge != no_cell) {
        out.edge_h0_negative_mask = bit(h0_edge);
    }

    const std::uint8_t positive_edges =
        static_cast<std::uint8_t>(out.edge_mask & ~out.edge_h0_negative_mask);
    std::array<EdgeColumn, 4> rewrite{};
    std::uint8_t has_rewrite = 0;

    for (std::uint8_t s : out.square_order_asc) {
        if (s == no_cell) {
            continue;
        }
        const std::uint8_t original =
            static_cast<std::uint8_t>(square_boundary_edges_2d(s) & positive_edges);
        std::uint8_t column = original;
        while (column != 0) {
            const std::uint8_t pivot =
                youngest_edge_2d(column, out.edge_second_rank);
            if ((has_rewrite & bit(pivot)) == 0) {
                break;
            }
            column = static_cast<std::uint8_t>(
                (column & ~bit(pivot)) ^ rewrite[pivot].bits);
        }
        if (column == 0) {
            continue;
        }

        const std::uint8_t pivot =
            youngest_edge_2d(column, out.edge_second_rank);
        out.edge_h1_positive_mask |= bit(pivot);
        out.square_h1_negative_mask |= bit(s);
        if ((original & bit(pivot)) != 0 && out.direct_h1_dual_count < 4) {
            out.direct_h1_dual_pairs[out.direct_h1_dual_count++] = {pivot, s};
        }
        has_rewrite |= bit(pivot);
        rewrite[pivot].bits = static_cast<std::uint8_t>(column & ~bit(pivot));
    }

    out.survivor_edge_mask =
        static_cast<std::uint8_t>(out.edge_mask &
                                  ~out.edge_h0_negative_mask &
                                  ~out.edge_h1_positive_mask);
    out.survivor_square_mask =
        static_cast<std::uint8_t>(out.square_mask & ~out.square_h1_negative_mask);
    out.edge_h1_positive_count = popcount4(out.edge_h1_positive_mask);
    out.square_h1_negative_count = popcount4(out.square_h1_negative_mask);
    out.survivor_edge_count = popcount4(out.survivor_edge_mask);
    out.survivor_square_count = popcount4(out.survivor_square_mask);
    for (std::uint8_t e : out.edge_order_asc) {
        if (e != no_cell && (out.survivor_edge_mask & bit(e)) != 0) {
            out.residual_edge_order_asc[out.residual_edge_order_asc_count++] = e;
        }
    }
    for (std::uint8_t e : out.edge_order_desc) {
        if (e != no_cell && (out.survivor_edge_mask & bit(e)) != 0) {
            out.residual_edge_order_desc[out.residual_edge_order_desc_count++] = e;
        }
    }
    for (std::uint8_t s : out.square_order_asc) {
        if (s != no_cell && (out.survivor_square_mask & bit(s)) != 0) {
            out.residual_square_order_asc[out.residual_square_order_asc_count++] = s;
        }
    }
    return out;
}

std::array<Entry2D, 256> build_table2d() {
    std::array<Entry2D, 256> table{};
    for (std::size_t i = 0; i < table.size(); ++i) {
        table[i] = build_entry2d(static_cast<std::uint8_t>(i));
    }
    return table;
}

const std::array<Entry2D, 256>& table2d() {
    static const std::array<Entry2D, 256> table = build_table2d();
    return table;
}

int local_vertex_id_3d(int dx, int dy, int dz) {
    const int lx = dx + 1;
    const int ly = dy + 1;
    const int lz = dz + 1;
    return (lz * 3 + ly) * 3 + lx;
}

int edge_source_vertex_3d(std::uint8_t e) {
    switch (e) {
    case edge3_x_neg:
        return local_vertex_id_3d(-1, 0, 0);
    case edge3_x_pos:
        return local_vertex_id_3d(0, 0, 0);
    case edge3_y_neg:
        return local_vertex_id_3d(0, -1, 0);
    case edge3_y_pos:
        return local_vertex_id_3d(0, 0, 0);
    case edge3_z_neg:
        return local_vertex_id_3d(0, 0, -1);
    case edge3_z_pos:
        return local_vertex_id_3d(0, 0, 0);
    default:
        return 0;
    }
}

int edge_dir_3d(std::uint8_t e) {
    switch (e) {
    case edge3_x_neg:
    case edge3_x_pos:
        return 0;
    case edge3_y_neg:
    case edge3_y_pos:
        return 1;
    default:
        return 2;
    }
}

int edge_id_3d(std::uint8_t e) {
    return 3 * edge_source_vertex_3d(e) + edge_dir_3d(e);
}

int sign_bit(bool positive) {
    return positive ? 1 : 0;
}

std::uint8_t edge_x(bool positive) {
    return positive ? edge3_x_pos : edge3_x_neg;
}

std::uint8_t edge_y(bool positive) {
    return positive ? edge3_y_pos : edge3_y_neg;
}

std::uint8_t edge_z(bool positive) {
    return positive ? edge3_z_pos : edge3_z_neg;
}

std::uint8_t square_xy(bool xp, bool yp) {
    return static_cast<std::uint8_t>(sign_bit(xp) + 2 * sign_bit(yp));
}

std::uint8_t square_xz(bool xp, bool zp) {
    return static_cast<std::uint8_t>(4 + sign_bit(xp) + 2 * sign_bit(zp));
}

std::uint8_t square_yz(bool yp, bool zp) {
    return static_cast<std::uint8_t>(8 + sign_bit(yp) + 2 * sign_bit(zp));
}

int square_source_vertex_3d(std::uint8_t s) {
    if (s < 4) {
        const bool xp = (s & 1u) != 0;
        const bool yp = (s & 2u) != 0;
        return local_vertex_id_3d(xp ? 0 : -1, yp ? 0 : -1, 0);
    }
    if (s < 8) {
        const std::uint8_t t = static_cast<std::uint8_t>(s - 4);
        const bool xp = (t & 1u) != 0;
        const bool zp = (t & 2u) != 0;
        return local_vertex_id_3d(xp ? 0 : -1, 0, zp ? 0 : -1);
    }
    const std::uint8_t t = static_cast<std::uint8_t>(s - 8);
    const bool yp = (t & 1u) != 0;
    const bool zp = (t & 2u) != 0;
    return local_vertex_id_3d(0, yp ? 0 : -1, zp ? 0 : -1);
}

int square_orient_3d(std::uint8_t s) {
    if (s < 4) {
        return 0;
    }
    if (s < 8) {
        return 1;
    }
    return 2;
}

int square_id_3d(std::uint8_t s) {
    return 3 * square_source_vertex_3d(s) + square_orient_3d(s);
}

std::uint8_t square_boundary_edges_3d(std::uint8_t s) {
    if (s < 4) {
        const bool xp = (s & 1u) != 0;
        const bool yp = (s & 2u) != 0;
        return static_cast<std::uint8_t>(bit(edge_x(xp)) | bit(edge_y(yp)));
    }
    if (s < 8) {
        const std::uint8_t t = static_cast<std::uint8_t>(s - 4);
        const bool xp = (t & 1u) != 0;
        const bool zp = (t & 2u) != 0;
        return static_cast<std::uint8_t>(bit(edge_x(xp)) | bit(edge_z(zp)));
    }
    const std::uint8_t t = static_cast<std::uint8_t>(s - 8);
    const bool yp = (t & 1u) != 0;
    const bool zp = (t & 2u) != 0;
    return static_cast<std::uint8_t>(bit(edge_y(yp)) | bit(edge_z(zp)));
}

int square_birth_edge_id_3d(std::uint8_t s) {
    const std::uint8_t boundary = square_boundary_edges_3d(s);
    int best = -1;
    for (std::uint8_t e = 0; e < 6; ++e) {
        if ((boundary & bit(e)) != 0) {
            best = std::max(best, edge_id_3d(e));
        }
    }
    return best;
}

bool square_less_3d(std::uint8_t a, std::uint8_t b) {
    const int ea = square_birth_edge_id_3d(a);
    const int eb = square_birth_edge_id_3d(b);
    if (ea != eb) {
        return ea < eb;
    }
    return square_id_3d(a) < square_id_3d(b);
}

int cube_id_3d(std::uint8_t c) {
    return c;
}

std::uint16_t cube_boundary_squares_3d(std::uint8_t c) {
    const bool xp = (c & 1u) != 0;
    const bool yp = (c & 2u) != 0;
    const bool zp = (c & 4u) != 0;
    return static_cast<std::uint16_t>(
        bit16(square_xy(xp, yp)) |
        bit16(square_xz(xp, zp)) |
        bit16(square_yz(yp, zp)));
}

std::uint8_t cube_birth_square_3d(std::uint8_t c) {
    const std::uint16_t boundary = cube_boundary_squares_3d(c);
    std::uint8_t best = no_cell;
    for (std::uint8_t s = 0; s < 12; ++s) {
        if ((boundary & bit16(s)) == 0) {
            continue;
        }
        if (best == no_cell || square_less_3d(best, s)) {
            best = s;
        }
    }
    return best;
}

bool cube_less_3d(std::uint8_t a, std::uint8_t b) {
    const std::uint8_t sa = cube_birth_square_3d(a);
    const std::uint8_t sb = cube_birth_square_3d(b);
    if (square_less_3d(sa, sb)) {
        return true;
    }
    if (square_less_3d(sb, sa)) {
        return false;
    }
    return cube_id_3d(a) < cube_id_3d(b);
}

std::uint8_t youngest_edge_3d(std::uint8_t mask) {
    std::uint8_t best = no_cell;
    for (std::uint8_t e = 0; e < 6; ++e) {
        if ((mask & bit(e)) == 0) {
            continue;
        }
        if (best == no_cell || edge_id_3d(best) < edge_id_3d(e)) {
            best = e;
        }
    }
    return best;
}

std::uint8_t youngest_square_3d(std::uint16_t mask) {
    std::uint8_t best = no_cell;
    for (std::uint8_t s = 0; s < 12; ++s) {
        if ((mask & bit16(s)) == 0) {
            continue;
        }
        if (best == no_cell || square_less_3d(best, s)) {
            best = s;
        }
    }
    return best;
}

std::uint16_t square_mask_from_af(std::uint8_t edges, std::uint16_t faces) {
    std::uint16_t squares = 0;
    for (std::uint8_t s = 0; s < 12; ++s) {
        const std::uint8_t boundary = square_boundary_edges_3d(s);
        if ((faces & bit16(s)) != 0 &&
            (edges & boundary) == boundary) {
            squares |= bit16(s);
        }
    }
    return squares;
}

std::uint8_t cube_mask_from_esb(std::uint8_t edges,
                                std::uint16_t squares,
                                std::uint8_t bodies) {
    (void)edges;
    std::uint8_t cubes = 0;
    for (std::uint8_t c = 0; c < 8; ++c) {
        const bool xp = (c & 1u) != 0;
        const bool yp = (c & 2u) != 0;
        const bool zp = (c & 4u) != 0;
        const std::uint8_t required_edges =
            static_cast<std::uint8_t>(bit(edge_x(xp)) |
                                      bit(edge_y(yp)) |
                                      bit(edge_z(zp)));
        const std::uint16_t required_squares = cube_boundary_squares_3d(c);
        if ((bodies & bit(c)) != 0 &&
            (edges & required_edges) == required_edges &&
            (squares & required_squares) == required_squares) {
            cubes |= bit(c);
        }
    }
    return cubes;
}

void fill_ordering_metadata_3d(Entry3D& out) {
    out.edge_second_rank.fill(no_cell);
    out.square_second_rank.fill(no_cell);
    out.cube_second_rank.fill(no_cell);
    out.square_birth_local_edge.fill(no_cell);
    out.cube_birth_local_square.fill(no_cell);
    out.edge_order_asc.fill(no_cell);
    out.edge_order_desc.fill(no_cell);
    out.square_order_asc.fill(no_cell);
    out.square_order_desc.fill(no_cell);
    out.cube_order_asc.fill(no_cell);
    out.cube_order_desc.fill(no_cell);
    out.residual_edge_order_asc.fill(no_cell);
    out.residual_edge_order_desc.fill(no_cell);
    out.residual_square_order_asc.fill(no_cell);
    out.residual_square_order_desc.fill(no_cell);
    out.h2_residual_square_order_asc.fill(no_cell);
    out.residual_cube_order_asc.fill(no_cell);

    std::array<std::uint8_t, 6> edge_sorted{};
    out.edge_count = 0;
    for (std::uint8_t e = 0; e < 6; ++e) {
        if ((out.edge_mask & bit(e)) != 0) {
            edge_sorted[out.edge_count++] = e;
        }
    }
    std::sort(edge_sorted.begin(), edge_sorted.begin() + out.edge_count,
              [](std::uint8_t a, std::uint8_t b) {
                  return edge_id_3d(a) < edge_id_3d(b);
              });
    for (std::uint8_t i = 0; i < out.edge_count; ++i) {
        const std::uint8_t e = edge_sorted[i];
        out.edge_second_rank[e] = i;
        out.edge_order_asc[i] = e;
        out.edge_order_desc[out.edge_count - 1u - i] = e;
    }

    std::array<std::uint8_t, 12> square_sorted{};
    out.square_count = 0;
    for (std::uint8_t s = 0; s < 12; ++s) {
        if ((out.square_mask & bit16(s)) != 0) {
            square_sorted[out.square_count++] = s;
        }
    }
    std::sort(square_sorted.begin(), square_sorted.begin() + out.square_count,
              [](std::uint8_t a, std::uint8_t b) {
                  return square_less_3d(a, b);
              });
    for (std::uint8_t i = 0; i < out.square_count; ++i) {
        const std::uint8_t s = square_sorted[i];
        out.square_second_rank[s] = i;
        out.square_order_asc[i] = s;
        out.square_order_desc[out.square_count - 1u - i] = s;
    }

    std::array<std::uint8_t, 8> cube_sorted{};
    out.cube_count = 0;
    for (std::uint8_t c = 0; c < 8; ++c) {
        if ((out.cube_mask & bit(c)) != 0) {
            cube_sorted[out.cube_count++] = c;
        }
    }
    std::sort(cube_sorted.begin(), cube_sorted.begin() + out.cube_count,
              [](std::uint8_t a, std::uint8_t b) {
                  return cube_less_3d(a, b);
              });
    for (std::uint8_t i = 0; i < out.cube_count; ++i) {
        const std::uint8_t c = cube_sorted[i];
        out.cube_second_rank[c] = i;
        out.cube_order_asc[i] = c;
        out.cube_order_desc[out.cube_count - 1u - i] = c;
    }

    for (std::uint8_t s = 0; s < 12; ++s) {
        if ((out.square_mask & bit16(s)) == 0) {
            continue;
        }
        const std::uint8_t boundary = square_boundary_edges_3d(s);
        std::uint8_t best = no_cell;
        for (std::uint8_t e = 0; e < 6; ++e) {
            if ((boundary & bit(e)) == 0) {
                continue;
            }
            if (best == no_cell || edge_id_3d(best) < edge_id_3d(e)) {
                best = e;
            }
        }
        out.square_birth_local_edge[s] = best;
    }

    for (std::uint8_t c = 0; c < 8; ++c) {
        if ((out.cube_mask & bit(c)) == 0) {
            continue;
        }
        const std::uint16_t boundary = cube_boundary_squares_3d(c);
        std::uint8_t best = no_cell;
        for (std::uint8_t s = 0; s < 12; ++s) {
            if ((boundary & bit16(s)) == 0) {
                continue;
            }
            if (best == no_cell || square_less_3d(best, s)) {
                best = s;
            }
        }
        out.cube_birth_local_square[c] = best;
    }

    for (std::uint8_t i = 0; i < out.edge_count; ++i) {
        const std::uint8_t e = out.edge_order_asc[i];
        if ((out.survivor_edge_mask & bit(e)) != 0) {
            out.residual_edge_order_asc[out.residual_edge_asc_count++] = e;
        }
    }
    for (std::uint8_t i = 0; i < out.edge_count; ++i) {
        const std::uint8_t e = out.edge_order_desc[i];
        if ((out.survivor_edge_mask & bit(e)) != 0) {
            out.residual_edge_order_desc[out.residual_edge_desc_count++] = e;
        }
    }

    for (std::uint8_t i = 0; i < out.square_count; ++i) {
        const std::uint8_t s = out.square_order_asc[i];
        if ((out.survivor_square_mask & bit16(s)) != 0) {
            out.residual_square_order_asc[out.residual_square_asc_count++] = s;
        }
    }
    for (std::uint8_t i = 0; i < out.square_count; ++i) {
        const std::uint8_t s = out.square_order_desc[i];
        if ((out.survivor_square_mask & bit16(s)) != 0) {
            out.residual_square_order_desc[out.residual_square_desc_count++] = s;
        }
    }

    for (std::uint8_t i = 0; i < out.square_count; ++i) {
        const std::uint8_t s = out.square_order_asc[i];
        if ((out.h2_residual_square_mask & bit16(s)) != 0) {
            out.h2_residual_square_order_asc[out.h2_residual_square_asc_count++] = s;
        }
    }

    for (std::uint8_t i = 0; i < out.cube_count; ++i) {
        const std::uint8_t c = out.cube_order_asc[i];
        if ((out.survivor_cube_mask & bit(c)) != 0) {
            out.residual_cube_order_asc[out.residual_cube_asc_count++] = c;
        }
    }
}

Entry3D build_entry3d(std::uint8_t edge_mask,
                      std::uint16_t square_mask,
                      std::uint8_t cube_mask) {
    Entry3D out{};
    out.edge_mask = edge_mask;
    out.square_mask = square_mask;
    out.cube_mask = cube_mask;

    std::uint8_t h0_edge = no_cell;
    for (std::uint8_t e = 0; e < 6; ++e) {
        if ((edge_mask & bit(e)) == 0) {
            continue;
        }
        if (h0_edge == no_cell || edge_id_3d(e) < edge_id_3d(h0_edge)) {
            h0_edge = e;
        }
    }
    out.direct_h0_edge = h0_edge;
    if (h0_edge != no_cell) {
        out.edge_h0_negative_mask = bit(h0_edge);
    }

    const std::uint8_t positive_edges =
        static_cast<std::uint8_t>(edge_mask & ~out.edge_h0_negative_mask);
    std::array<EdgeColumn, 6> edge_rewrite{};
    std::array<std::uint16_t, 6> edge_rewrite_substitutions{};
    std::uint8_t has_edge_rewrite = 0;

    std::array<std::uint8_t, 12> square_order{};
    for (std::uint8_t i = 0; i < 12; ++i) {
        square_order[i] = i;
    }
    std::sort(square_order.begin(), square_order.end(),
              [](std::uint8_t a, std::uint8_t b) {
                  return square_less_3d(a, b);
              });

    for (std::uint8_t s : square_order) {
        if ((square_mask & bit16(s)) == 0) {
            continue;
        }
        const std::uint8_t original =
            static_cast<std::uint8_t>(square_boundary_edges_3d(s) & positive_edges);
        std::uint8_t column = original;
        std::uint16_t substituted_squares = 0;
        while (column != 0) {
            const std::uint8_t pivot = youngest_edge_3d(column);
            if ((has_edge_rewrite & bit(pivot)) == 0) {
                break;
            }
            column = static_cast<std::uint8_t>(
                (column & ~bit(pivot)) ^ edge_rewrite[pivot].bits);
            substituted_squares = static_cast<std::uint16_t>(
                substituted_squares ^ edge_rewrite_substitutions[pivot]);
        }
        out.h1_square_substitution_square_mask[s] = substituted_squares;
        if (column == 0) {
            continue;
        }
        out.h1_square_semireduced_boundary[s] = column;
        const std::uint8_t pivot = youngest_edge_3d(column);
        out.edge_h1_positive_mask |= bit(pivot);
        out.square_h1_negative_mask |= bit16(s);
        out.square_h2_skip_mask |= bit16(s);
        if (out.zero_h1_count < 12) {
            const std::uint8_t pair_index = out.zero_h1_count++;
            out.zero_h1_pairs[pair_index] = {pivot, s};
            out.zero_h1_rewrite_masks[pair_index] =
                static_cast<std::uint8_t>(column & ~bit(pivot));
        }
        has_edge_rewrite |= bit(pivot);
        edge_rewrite[pivot].bits =
            static_cast<std::uint8_t>(column & ~bit(pivot));
        edge_rewrite_substitutions[pivot] =
            static_cast<std::uint16_t>(substituted_squares ^ bit16(s));
    }

    const std::uint16_t positive_squares =
        static_cast<std::uint16_t>(square_mask & ~out.square_h1_negative_mask);
    std::array<SquareColumn, 12> square_rewrite{};
    std::uint16_t has_square_rewrite = 0;
    std::uint8_t h2_negative_cube_mask = 0;
    std::uint16_t direct_h2_square_mask = 0;

    std::array<std::uint8_t, 8> cube_order{};
    for (std::uint8_t i = 0; i < 8; ++i) {
        cube_order[i] = i;
    }
    std::sort(cube_order.begin(), cube_order.end(),
              [](std::uint8_t a, std::uint8_t b) {
                  return cube_less_3d(a, b);
              });

    for (std::uint8_t c : cube_order) {
        if ((cube_mask & bit(c)) == 0) {
            continue;
        }
        const std::uint16_t original =
            static_cast<std::uint16_t>(cube_boundary_squares_3d(c) & positive_squares);
        std::uint16_t column = original;
        while (column != 0) {
            const std::uint8_t pivot = youngest_square_3d(column);
            if ((has_square_rewrite & bit16(pivot)) == 0) {
                break;
            }
            column = static_cast<std::uint16_t>(
                (column & ~bit16(pivot)) ^ square_rewrite[pivot].bits);
        }
        if (column == 0) {
            continue;
        }
        const std::uint8_t pivot = youngest_square_3d(column);
        out.square_h2_negative_mask |= bit16(pivot);
        h2_negative_cube_mask |= bit(c);
        if ((original & bit16(pivot)) != 0 && out.direct_h2_count < 8) {
            out.direct_h2_pairs[out.direct_h2_count++] = {pivot, c};
            direct_h2_square_mask =
                static_cast<std::uint16_t>(direct_h2_square_mask | bit16(pivot));
        }
        has_square_rewrite |= bit16(pivot);
        square_rewrite[pivot].bits =
            static_cast<std::uint16_t>(column & ~bit16(pivot));
    }

    out.survivor_edge_mask =
        static_cast<std::uint8_t>(edge_mask &
                                  ~out.edge_h0_negative_mask &
                                  ~out.edge_h1_positive_mask);
    out.survivor_square_mask =
        static_cast<std::uint16_t>(square_mask &
                                   ~out.square_h1_negative_mask &
                                   ~out.square_h2_negative_mask);
    out.survivor_cube_mask =
        static_cast<std::uint8_t>(cube_mask & ~h2_negative_cube_mask);
    out.h2_residual_square_mask =
        static_cast<std::uint16_t>(square_mask &
                                   ~out.square_h1_negative_mask &
                                   ~direct_h2_square_mask);
    fill_ordering_metadata_3d(out);
    return out;
}

class ZeroTable3DData {
public:
    ZeroTable3DData() {
        build();
    }

    const Entry3D& lookup(std::uint32_t mask26) const {
        const std::uint32_t a = mask26 & 0x3fu;
        const std::uint32_t f = (mask26 >> 6) & 0xfffu;
        const std::uint32_t b = (mask26 >> 18) & 0xffu;
        const std::uint16_t es = es_[a * 4096u + f];
        const std::uint16_t geom = full_[static_cast<std::size_t>(es) * 256u + b];
        return pairs_[geom];
    }

    Table3DStats stats() const {
        return {es_states_.size(),
                pairs_.size(),
                es_.size() * sizeof(std::uint16_t),
                full_.size() * sizeof(std::uint16_t),
                pairs_.size() * sizeof(Entry3D)};
    }

    std::size_t pairs_count() const { return pairs_.size(); }
    const Entry3D& pairs_at(std::size_t idx) const { return pairs_[idx]; }

private:
    struct ESState {
        std::uint8_t edges = 0;
        std::uint16_t squares = 0;
    };

    std::vector<std::uint16_t> es_;
    std::vector<std::uint16_t> full_;
    std::vector<Entry3D> pairs_;
    std::vector<ESState> es_states_;

    void build() {
        es_.resize(64u * 4096u);
        std::unordered_map<std::uint32_t, std::uint16_t> es_ids;
        es_ids.reserve(8192);
        for (std::uint32_t a = 0; a < 64u; ++a) {
            for (std::uint32_t f = 0; f < 4096u; ++f) {
                const std::uint8_t edges = static_cast<std::uint8_t>(a);
                const std::uint16_t squares =
                    square_mask_from_af(edges, static_cast<std::uint16_t>(f));
                const std::uint32_t key =
                    static_cast<std::uint32_t>(edges) |
                    (static_cast<std::uint32_t>(squares) << 6);
                auto it = es_ids.find(key);
                if (it == es_ids.end()) {
                    if (es_states_.size() > 0xffffu) {
                        throw std::runtime_error("too many 3D edge-square states");
                    }
                    const std::uint16_t id =
                        static_cast<std::uint16_t>(es_states_.size());
                    it = es_ids.emplace(key, id).first;
                    es_states_.push_back({edges, squares});
                }
                es_[a * 4096u + f] = it->second;
            }
        }

        full_.resize(es_states_.size() * 256u);
        std::unordered_map<std::uint32_t, std::uint16_t> geom_ids;
        geom_ids.reserve(20000);
        for (std::size_t es_id = 0; es_id < es_states_.size(); ++es_id) {
            const ESState state = es_states_[es_id];
            for (std::uint32_t b = 0; b < 256u; ++b) {
                const std::uint8_t cubes =
                    cube_mask_from_esb(state.edges, state.squares,
                                       static_cast<std::uint8_t>(b));
                const std::uint32_t key =
                    static_cast<std::uint32_t>(state.edges) |
                    (static_cast<std::uint32_t>(state.squares) << 6) |
                    (static_cast<std::uint32_t>(cubes) << 18);
                auto it = geom_ids.find(key);
                if (it == geom_ids.end()) {
                    if (pairs_.size() > 0xffffu) {
                        throw std::runtime_error("too many 3D zero table states");
                    }
                    const std::uint16_t id =
                        static_cast<std::uint16_t>(pairs_.size());
                    it = geom_ids.emplace(key, id).first;
                    pairs_.push_back(build_entry3d(state.edges, state.squares, cubes));
                }
                full_[es_id * 256u + b] = it->second;
            }
        }
    }
};

const ZeroTable3DData& table3d() {
    static const ZeroTable3DData table;
    return table;
}

}

const Entry2D& lookup2d(std::uint8_t mask) {
    return table2d()[mask];
}

const Entry3D& lookup3d(std::uint32_t mask26) {
    return table3d().lookup(mask26);
}

void warm_zero_tables() {
    warm_zero_table2d();
    warm_zero_table3d();
}

void warm_zero_table2d() {
    (void)table2d();
}

void warm_zero_table3d() {
    (void)table3d();
}

Table3DStats table3d_stats() {
    return table3d().stats();
}

std::size_t entry3d_count() {
    return table3d().pairs_count();
}

const Entry3D& entry3d_at(std::size_t idx) {
    return table3d().pairs_at(idx);
}

}
