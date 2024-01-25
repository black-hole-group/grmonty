GRMONTY: A relativistic Monte Carlo code
==========================================

Based on [Dolence et al. 2009 ApJ](http://adsabs.harvard.edu/abs/2009ApJS..184..387D). Originally downloaded from [Astrophysical Code Library](http://rainman.astro.illinois.edu/codelib/) @ [UI](http://illinois.edu).

GRMONTY is parallelized using [OpenMP](https://en.wikipedia.org/wiki/OpenMP). This version is configured to use input from [`harm2d`](http://rainman.astro.illinois.edu/codelib/codes/ham2d/src/).

# H-AMR branch
In this branch, we tried to make the code work with data from the state-of-the-art code H-AMR ([Liska et al. 2019](https://arxiv.org/abs/1912.10192)).

## H-AMR data
We do not use H-AMR data directly, we use a notebook to convert H-AMR's dump files into a binary file with the components in the right order for grmonty to handle. In this example, I've but the file HAMR_GRMONTY_DUMP323.bin as an example. This is a 2D simulation with dimensions $(256 \times 256)$ in $r - \theta$ dimensions.

## Changes in the code
H-AMR deals with the $x_2$ dimension differently than HARM's usual correlation formula: $\theta = \pi x_2 + \frac{1}{2}(1-h)sin(2\pi x_2)$, which maps $0 \leq \theta \leq \pi$ to $0 \leq x_2 \leq 1$, because of this, we use a different function to convert $(x_1, x_2, x_3)$ to $(r, \theta, \phi)$. We also used different functions to calculate:

* Connection components $(\Gamma^\alpha_{\beta \gamma})$ are calculated in function ```conn_func```
* Gcov $(g_{\mu\nu})$ and Gcon $(g^{\mu \nu})$ components are calculated in functions ```gcov_func_hamr``` and ```gcon_func_hamr```

set number of threads for `csh` and 8 threads:

    setenv OMP_NUM_THREADS 8

if using `bash`:

    export OMP_NUM_THREADS=8

run the code on the supplied harm output file:

    ./grmonty 5000000 dump019 4.e19 


This will output the spectrum to `grmonty.spec`  which should be identical to `grmonty_spec_verify`.

# Plotting

Use python and the [`nmmn`](https://github.com/rsnemmen/nmmn) module:

```python
import nmmn.sed
s=nmmn.sed.SED()
s.grmonty('grmonty.spec')
plot(s.lognu, s.ll)
```

Old-fashioned way: Use the [SM](http://www.astro.princeton.edu/~rhl/sm/) script `plspec.m` to plot up broad-band spectral energy distribution.

# Calculate spectra from other sources

Replace `harm_model.c` with your own source model.  Begin by modifying `harm_model.c`. You must supply

```
init_model 
make_super_photon
bias_func
get_fluid_params
report_spectrum
stop_criterion
record_criterion

gcon_func 
gcov_func 
get_connection
```

in the model file.

# TODO

- [ ] make it work with [HARMPI](https://github.com/atchekho/harmpi)
- [ ] GPU support: OpenCL
- [ ] parallelize with MPI
- [ ] add bremsstrahlung
- [ ] nonthermal electron distribution
- [ ] dynamic metrics as input
- [x] add LICENSE

# References

- Code and methods: [Dolence et al. 2009 ApJ](http://adsabs.harvard.edu/abs/2009ApJS..184..387D)
- An early application: Sgr A\* SED model ([Moscibrodzka et al. 2009 ApJ](http://iopscience.iop.org/article/10.1088/0004-637X/706/1/497/meta)). A more recent Sgr A\* model by the same authors: [Moscibrodzka et al. 2014 A&A](http://www.aanda.org/articles/aa/abs/2014/10/aa24358-14/aa24358-14.html)
- More recent applications: M87 jet/RIAF SED ( [Moscibrodzka et al. 2016 A&A](http://www.aanda.org/articles/aa/abs/2016/02/aa26630-15/aa26630-15.html)), jet and RIAF SEDs for stellar mass black holes ([O'Riordan, Pe'er & McKinney 2016 ApJ](http://iopscience.iop.org/article/10.3847/0004-637X/819/2/95/meta))

# LICENSE 

`grmonty` is free software: you can redistribute it and/or modify it under the terms of the GNU General Public License as published by the Free Software Foundation, either version 3 of the License, or (at your option) any later version.

This program is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the `LICENSE` file or the [GNU General Public License](http://www.gnu.org/licenses/) for more details.
