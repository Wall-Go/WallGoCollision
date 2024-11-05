<img src="docs/source/figures/wallgo.svg" alt="WallGoLogo" width="200"/>

# WallGoCollision

Scientific library for computing Boltzmann collision integrals, written in C++ with Python bindings provided. **WallGoCollision** is part of the [**WallGo**](https://github.com/Wall-Go/WallGo) toolset for the computation of the bubble wall velocity and bubble wall width in first-order cosmological phase transitions.

Documentation for **WallGoCollision**: https://wallgocollision.readthedocs.io.

## Requirements

Python 3.10 or newer is required for the Python bindings. If building from source, you need a C++17 compliant compiler and CMake 3.18 or newer. See the [documentation](https://wallgocollision.readthedocs.io/en/latest/install.html) for more details on source builds.

## Installation

Stable releases of the Python bindings can be found on the [releases](https://github.com/Wall-Go/WallGoCollision/releases) page. These are also available on PyPI and can be installed with pip:

    pip install WallGoCollision

For source builds, see the [documentation](https://wallgocollision.readthedocs.io/en/latest/install.html).

## Quickstart

The [quickstart page in documentation](https://wallgocollision.readthedocs.io/en/latest/quickstart.html) should get you started. Example implementations of concrete physics models are available on the [**WallGo** repository](https://github.com/Wall-Go/WallGo/tree/main/Models).

## Feedback and further questions

For feedback and frequently asked questions, please see the
[WallGoCollision homepage](https://wallgocollision.readthedocs.io)
and the [WallGo homepage](https://wallgocollision.readthedocs.io).


## License

Copyright (c) 2024 Andreas Ekstedt, Oliver Gould, Joonas Hirvonen,
Benoit Laurent, Lauri Niemi, Philipp Schicho, and Jorinde van de Vis.

This program is free software: you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License.

This program is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with this program.  If not, see <https://www.gnu.org/licenses/>.
