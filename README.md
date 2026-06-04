# FlashCubical

Fast V-filtration cubical persistent homology for 2D and 3D scalar grids induced by images, impleneting the ideas proposed in https://arxiv.org/abs/2606.04801.

Computes H0 (connected components), H1 (loops / tunnels), and H2 (voids) using a lower-star filtration with a precomputed lookup table for local zero-persistence pairs.



Requirements

* Python 3.9+
* A C++ compiler with C++17 support (MSVC, GCC, or Clang)

## Installation

```
git clone git@github.com:T-prog123/FlashCubical.git
cd FlashCubical
pip install -e .
```

That is all. pip compiles the C++ extension automatically.

## Usage

```python
import flash_cubical as fc
import numpy as np

image = np.random.rand(256, 256)      # 2D scalar field
ph = fc.compute(image)

ph.values                             # float64 array, shape (n, 3)
                                      # columns: birth, death, dimension

ph.plot()                             # persistence diagram
```

For a 3D volume:

```python
volume = np.random.rand(64, 64, 64)
ph = fc.compute(volume)               # computes H0, H1, H2
ph = fc.compute(volume, h1=False)     # skip H1, only H0 and H2
```

## API

### `fc.compute(x, min_persistence=0.0, h1=True)`

|Parameter|Description|
|-|-|
|`x`|2D or 3D array-like, converted to float64|
|`min_persistence`|Minimal distance between birth and death values. Must be non-negative.|
|`h1`|Compute H1. Must be `True` for 2D inputs. Can be `False` for 3D to skip H1|

Returns a `Persistence` object.

### `ph.values`

NumPy array of shape `(n, 3)`, dtype float64.

|Column|Meaning|
|-|-|
|0|Birth|
|1|Death (essential H0 class has death = +inf)|
|2|Homology dimension (0, 1, or 2)|

### `ph.plot(ax=None)`

Draws a persistence diagram coloured by dimension. Returns a matplotlib `Axes`.

## Running the tests

```
pip install pytest
pytest tests/
```

