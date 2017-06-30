GRMONTY: A relativistic Monte Carlo code
==========================================

Based on [Dolence et al. 2009 ApJ](http://adsabs.harvard.edu/abs/2009ApJS..184..387D). Originally downloaded from [Astrophysical Code Library](http://rainman.astro.illinois.edu/codelib/) @ [UI](http://illinois.edu).

GRMONTY is parallelized using [OpenMP](https://en.wikipedia.org/wiki/OpenMP). This version is configured to use input from [`harm2d`](http://rainman.astro.illinois.edu/codelib/codes/ham2d/src/).

# quick start

unpack the tarball:

    tar -xzvf grmonty.tgz

make (requires openmp enabled gcc):

    make

set number of threads for `csh` and 8 threads:

    setenv OMP_NUM_THREADS 8

if using `bash`:

    export OMP_NUM_THREADS=8

run the code on the supplied harm output file:

    ./grmonty 5000000 dump019 4.e19 >&! grmonty.err

Arguments are:

- estimate of photon number (actual number is probabilistic due to scattering)
- harm dump file for model
- mass unit (few x 10^19 is appropriate for Sgr A*)

This will output spectrum to `grmonty.spec`  which should be identical to `grmonty_spec_verify`.

![GIF created with `ttygif`](./tty.gif)

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

<<<<<<< HEAD
- [ ] make it work with HARMPI/HARM3D
- [ ] make it work with RAISHIN
=======
- [ ] make it work with [HARMPI](https://github.com/atchekho/harmpi)
- [ ] GPU support: OpenCL
- [ ] parallelize with MPI
- [ ] add bremsstrahlung
- [ ] nonthermal electron distribution
- [ ] dynamic metrics as input
>>>>>>> progressbar
- [x] add LICENSE

# References

- Code and methods: [Dolence et al. 2009 ApJ](http://adsabs.harvard.edu/abs/2009ApJS..184..387D)
- An early application: Sgr A\* SED model ([Moscibrodzka et al. 2009 ApJ](http://iopscience.iop.org/article/10.1088/0004-637X/706/1/497/meta)). A more recent Sgr A\* model by the same authors: [Moscibrodzka et al. 2014 A&A](http://www.aanda.org/articles/aa/abs/2014/10/aa24358-14/aa24358-14.html)
- More recent applications: M87 jet/RIAF SED ( [Moscibrodzka et al. 2016 A&A](http://www.aanda.org/articles/aa/abs/2016/02/aa26630-15/aa26630-15.html)), jet and RIAF SEDs for stellar mass black holes ([O'Riordan, Pe'er & McKinney 2016 ApJ](http://iopscience.iop.org/article/10.3847/0004-637X/819/2/95/meta))

# LICENSE 

`grmonty` is free software: you can redistribute it and/or modify it under the terms of the GNU General Public License as published by the Free Software Foundation, either version 3 of the License, or (at your option) any later version.

This program is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the `LICENSE` file or the [GNU General Public License](http://www.gnu.org/licenses/) for more details.
