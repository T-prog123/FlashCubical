# Zero Lookup Runtime

This directory contains only the runtime side of the local zero-persistence
lookup tables used by the Python package.

Runtime is load-only. Normal package install, import, and computation load these
precomputed artifacts:

```text
src/flash_cubical/data/zero_lookup_2d.bin
src/flash_cubical/data/zero_lookup_3d.bin
```

If either file is missing, corrupted, or has an incompatible schema/table
version, the package raises an error. Runtime does not rebuild the tables.

## Runtime API

- `zero_lookup::load_from_files(path2d, path3d)`
- `zero_lookup::lookup2d(mask8)`
- `zero_lookup::lookup3d(mask26)`

Each vertex gets a local older-neighbour mask. The mask selects an `Entry2D` or
`Entry3D` record containing precomputed local lower-star information.

## File Format

The lookup files use an explicit little-endian binary format:

- fixed-width integer fields
- field-by-field serialization
- schema and table version metadata
- payload checksum

The format does not dump raw C++ structs or vectors, so it does not depend on
compiler-specific padding, alignment, or STL layout.

## Table Contents

The 2D file stores the dense 256-entry local lookup table.

The 3D file stores the split lookup layout:

- `ES`: axial plus face-diagonal state table
- `FULL`: edge/square state plus body-diagonal table
- `PAIRS`: compact zero-classification records

Important 3D fields include direct H0/H2 pairs, zero H1 rewrite data, local
cell ranks, local birth-cell metadata, and residual local orders.
