#pragma once

#include <cstdint>

namespace smart_core {

inline constexpr std::uint8_t edge_state_h0_negative = 1u << 0;
inline constexpr std::uint8_t edge_state_h1_positive = 1u << 1;
inline constexpr std::uint8_t edge_state_h0_apparent = 1u << 2;
inline constexpr std::uint8_t edge_state_h1_apparent = 1u << 3;

inline constexpr std::uint8_t square_state_dual_h0_negative = 1u << 0;
inline constexpr std::uint8_t square_state_dual_h1_positive = 1u << 1;
inline constexpr std::uint8_t square_state_h2_apparent = 1u << 2;
inline constexpr std::uint8_t square_state_dual_h1_apparent = 1u << 3;

} // namespace smart_core
