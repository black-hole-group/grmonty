GRMONTY: A relativistic Monte Carlo code
==========================================

Based on [Dolence et al. 2009 ApJ](http://adsabs.harvard.edu/abs/2009ApJS..184..387D). Originally downloaded from [Astrophysical Code Library](http://rainman.astro.illinois.edu/codelib/) @ [UI](http://illinois.edu).

GRMONTY is parallelized using [OpenMP](https://en.wikipedia.org/wiki/OpenMP). This version is configured to use input from [`harm2d`](http://rainman.astro.illinois.edu/codelib/codes/ham2d/src/).

# H-AMR branch
In this branch, we made the code work with data from the state-of-the-art code H-AMR ([Liska et al. 2019](https://arxiv.org/abs/1912.10192)).

## H-AMR data
We do not use H-AMR data directly, we use a notebook to convert H-AMR's dump files into a binary file with the components in the right order for grmonty to handle. In this example, I've but the file HAMR_GRMONTY_DUMP323.bin as an example. This is a 2D simulation with dimensions $(256 \times 256)$ in $r - \theta$ dimensions.

## Changes in the code
H-AMR deals with the $x_2$ dimension differently than HARM's usual correlation formula: $\theta = \pi x_2 + \frac{1}{2}(1-h)sin(2\pi x_2)$, which maps $0 \leq \theta \leq \pi$ to $0 \leq x_2 \leq 1$, because of this, we use a different function to convert $(x_1, x_2, x_3)$ to $(r, \theta, \phi)$. We also used different functions to calculate:

* Gcov $(g_{\mu\nu})$ and Gcon $(g^{\mu \nu})$ components are calculated in functions ```gcov_func_hamr``` and ```gcon_func_hamr```
* $x_1, x_2$ and $x_3$ coordinates based on the cell indexes in function ```coord_hamr```
* All the declarations are in harm_model.h file.

```p``` and ```geom```  were declared as nested arrays but we changed them by flattening the 3 dimensional spatial indexes into a 1 dimensional array with a 3D index, like this:
```p[NPRIM][i][j][k] -> p[NPRIM][SPATIAL_INDEX3D(i,j,k)]``` and also ```geom[i][j] -> geom[SPATIAL_INDEX2D(i,j)]```.

## Run the code
To run the code using H-AMR mode, you need to activate the switch within decs.h:72:

   ``` #define HAMR (0)```
   
change the number 0 to 1. If you want to run 3D data, please also activate the switch:
   ``` #define HAMR3D (1)```
Note that both switches must be activated for H-AMR 3D data.

Then you need to set the number of threads for the calculation using the right command listed below:

set number of threads for `csh` and 8 threads:

    setenv OMP_NUM_THREADS 8

if using `bash`:

    export OMP_NUM_THREADS=8

run the code on the supplied harm output file:

    ./grmonty 50000 ./data/HAMR_GRMONTY_DUMP323.bin 4.e19 


This will output the spectrum to `grmonty_hamr.spec` inside the output folder.

# LICENSE 

`grmonty` is free software: you can redistribute it and/or modify it under the terms of the GNU General Public License as published by the Free Software Foundation, either version 3 of the License, or (at your option) any later version.

This program is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the `LICENSE` file or the [GNU General Public License](http://www.gnu.org/licenses/) for more details.
