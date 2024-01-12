# Braid Visualiser

Simple python library for the creation and rendering of Artin braids.

## Installation

This package is available on PyPi and can be installed to your project with `pip`:

```
pip install braidvisualiser
```

## Example

- Rendering a simple 3-strand braid:

```
import braidvisualiser as bv

bv.Braid(3, 1, -2, 2)

bv.draw()
```

See the `draw` method's configurable parameters for rendering options.