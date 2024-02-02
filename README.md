# GRMONTY: A Relativistic Monte Carlo Code

GRMONTY is a Monte Carlo code designed for calculating the emergent spectrum from a model using a Monte Carlo technique. It is particularly suited for studying plasmas near black holes described by Kerr-Schild coordinates, radiating through thermal synchrotron and inverse Compton scattering. The code is based on the work presented in Dolence et al. 2009 Astrophysical Journal Supplement.

## Getting Started

To use GRMONTY, you can download the original version from the [Astrophysical Code Library](http://rainman.astro.illinois.edu/codelib/). The code is parallelized using [OpenMP](https://en.wikipedia.org/wiki/OpenMP) and is configured to use input files from the HARM code available on the same site.

## H-AMR Branch

This branch of GRMONTY has been modified to work with data from the H-AMR code ([Liska et al. 2019](https://arxiv.org/abs/1912.10192)). To use H-AMR data, a notebook is provided to convert H-AMR's dump files into a binary file with the components ordered appropriately for GRMONTY.

### H-AMR Data

We utilize a notebook to convert H-AMR's dump files into a binary file with the correct component order for GRMONTY. An example file, `HAMR_GRMONTY_DUMP323.bin`, is provided for a 2D simulation with dimensions $(256 \times 256)$ in $r - \theta$ dimensions. We do not provide the notebook for this conversion yet.

### Changes in the Code

Several modifications have been made to handle H-AMR data. Changes include adjustments to the conversion functions for spatial coordinates, differences in the correlation formula, and modifications in handling $x_1, x_2,$ and $x_3$ coordinates. Therefore, we use different functions to calculate:

* Gcov $(g_{\mu\nu})$ and Gcon $(g^{\mu \nu})$ components are calculated in functions ```gcov_func_hamr``` and ```gcon_func_hamr```
* $x_1, x_2$ and $x_3$ coordinates based on the cell indexes in function ```coord_hamr```
All the declarations are in harm_model.h file.

The pointers ```p``` and ```geom```  used to be declared as nested arrays but we changed them by flattening the 3-dimensional spatial indexes into a 1-dimensional array with a 3D index, like this:
```p[NPRIM][i][j][k] -> p[NPRIM][SPATIAL_INDEX3D(i,j,k)]``` and also ```geom[i][j] -> geom[SPATIAL_INDEX2D(i,j)]```.

## Run the Code
To run the code to read H-AMR data, switch on the corresponding flags in `decs.h`: 

``` #define HAMR (0)```
   
If you want to run 3D data, please also activate the switch:

   ``` #define HAMR3D (1)```
   
Note that both switches must be activated for H-AMR 3D data.

Compile the code by using

```
make
```

Set number of threads for bash

```
export OMP_NUM_THREADS=8
```

```
./grmonty 50000 ./data/HAMR_GRMONTY_DUMP323.bin 3.2e10
```


# LICENSE 

`grmonty` is free software: you can redistribute it and/or modify it under the terms of the GNU General Public License as published by the Free Software Foundation, either version 3 of the License, or (at your option) any later version.

This program is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the `LICENSE` file or the [GNU General Public License](http://www.gnu.org/licenses/) for more details.
