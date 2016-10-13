## *lautils-mat*: Luigi Acerbi's utility functions for MATLAB

This is a collection of MATLAB utility functions that I wrote (unless otherwise specified), and that I share among a number of projects. Functions in the following list have been consistently used or at least somewhat tested. Undocumented functions in the repository might be untested/deprecated, so use only at your own risk.

Functions marked with (+) are pretty useful functions that I use regularly.
Functions marked with (++) have been game-changers for several projects.

Functions are loosely divided in folders by category.

### Math (`./math`)

Numerical and math-related functions. 

- (++) `bsxfunandsum.m`:  Element-by-element matrix operation followed by sum (prevents memory overflow with huge matrices)
- `derivcheck.m`: Check analytical vs numerical differentiation for a function
- `fgrad.m`: Estimation of function gradient via finite differences
- `fhess.m`: Estimation of function Hessian via finite differences
- `intpower.m`: Array integer power (faster than `.^` with large matrices)
- (+) `lininterp1.m`: Linear 1D interpolation on a regular grid (like `interp1.m`, but faster)
- `logsumexp.m`: Compute `log(sum(exp(X)))` while avoiding numerical underflow
- `nansumall.m`: Sum of all elements, ignoring NaN values
- (+) `qcumtrapz.m`: Quick cumulative trapezoidal numerical integration
- (++) `qtrapz.m`: Quick trapezoidal numerical integration (3-4 times faster than `trapz.m` on large arrays)

### Markov Chain Monte Carlo sampling (`./mcmc`)

MCMC samplers and related files.

- `maxsample.m`: Stochastic exploration of high-valued regions of a function (untested)
- `mhsamplecon.m`: Metropolis-Hastings with reflective constraints (untested)

### Optimization (`./optim`)

Optimization algorithms and related functions.

- `cmaes_wrapper.m`: A wrapper with a simplified interface for the CMA-ES algorithm (requires [`cmaes.m`](https://www.lri.fr/~hansen/cmaes_inmatlab.html))
- `coordtransf.m`: Coordinate transform
- `fminfold.m`: Multi-fold optimization (work in progress); uses `fminfold_wrap.m`
- `fmingrid.m`: Multidimensional minimization via grid search
- (+) `fminmulti.m`: Function minimization via recursive multi-start method (usable but needs improvement)
- (+) `funlogger.m`: Record log of calls to a specific function
- `lhs.m`: Latin hypercube sample
- (+) `qargmax1`: Quick and dirty numerical argmax via 1-D quadratic interpolation
- `transvars.m`: Another coordinate transform function

### Plotting (`./plot`)

Plotting and graphics functions.

- `graph.m`: Prepare a boxed, annotated graph
- `multigraph.m`: Prepare a graph with multiple panels (usable, but needs rehauling)
- `scatterbars.m`: Scatter/bubble plot with error bars (unused?)

### PsychToolbox Add-ons (`./psych`)

A bunch of functions for [PsychToolbox](http://psychtoolbox.org/) that I wrote during my PhD. Nothing particularly useful.

### Statistics (`./stats`)

Probability distributions and other statistics-related functions.

- `binbuild.m`: Build binned data
- `binolike.m`: Negative log-likelihood for the binomial distribution
- (++) `bsxfun_normcdf.m`, `bsxfun_normlogpdf.m`, `bsxfun_normpdf.m`: Vectorized normal cdf, log pdf and pdf
- (+) `fisher2kappa.m`: Fisher information to Von Mises concentration parameter
- `gaussks1.m`: 1-D Gaussian kernel smooth function approximation
- `gpinterp1.m`: 1-D interpolation using Gaussian process (requires [GPML toolbox](http://www.gaussianprocess.org/gpml/code/matlab/doc/))
- `moments.m`: Central moments of discrete distribution
- (++) `qrandvm.m`: Quick Von Mises distributed pseudorandom numbers (faster than `circ_vmrnd.m` from the [Circular Statistics Toolbox](http://bethgelab.org/software/circstat/))
- `randntrim.m`: Normally distributed trimmed pseudorandom numbers (dumb rejection method)
- `randx.m`: Uniformly distributed non-overlapping random points
- `samplemeanpdf.m`: Sampling distribution of the mean (see also `samplemeanpdf_test.m`)
- `sqmean.m`, `sqsum.m`: Return root of squared mean/sum of data
- `stderr.m`: Standard error of the mean (ignores `NaN`s)

### Standard and input/output functions (`./stdio`)

Functions that deal with format conversion, argument parsing, structure manipulation and related miscellanea.
