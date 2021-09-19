# Monte-Carlo Simulations and Analysis of Stochastic Differential Equations

------------------------------------------------------------------------

GitHub      | MacOS      | Linux      | Windows    | Code Coverge | J. Stat. Softw 
------------|------------|------------|------------|--------------|-------------
[![Project Status](https://www.repostatus.org/badges/latest/active.svg?color=green)](https://github.com/acguidoum/Sim.DiffProc) | [![R](https://github.com/acguidoum/Sim.DiffProc/workflows/R-CMD-check/badge.svg?color=green)](https://github.com/acguidoum/Sim.DiffProc/actions) | [![Travis Build Status](https://img.shields.io/travis/com/acguidoum/Sim.DiffProc.svg?logo=travis)](https://travis-ci.org/acguidoum/Sim.DiffProc) | [![AppVeyor Build status](https://ci.appveyor.com/api/projects/status/16a70vyf8rk7nn1i?svg=true)](https://ci.appveyor.com/project/acguidoum/sim-diffproc-xal8n) | [![codecov](https://codecov.io/gh/acguidoum/Sim.DiffProc/branch/master/graph/badge.svg?color=brightgreen)](https://codecov.io/gh/acguidoum/Sim.DiffProc) | [![JSS](https://img.shields.io/badge/JSS-10.18637%2Fjss.v096.i02-brightgreen)](https://dx.doi.org/10.18637/jss.v096.i02)


------------------------------------------------------------------------

License    | Depends | Version   | Published | Ago
-----------|---------|-----------|-----------|----
[![CRAN](https://img.shields.io/cran/l/Sim.DiffProc?color=blue)](https://cran.r-project.org/web/licenses/GPL-2) | [![minimal R version](https://img.shields.io/badge/R%3E%3D-3.0.0-blue.svg)](https://cran.r-project.org/) | ![](https://www.r-pkg.org/badges/version/Sim.DiffProc?color=blue) | ![](https://www.r-pkg.org/badges/last-release/Sim.DiffProc?color=blue) | ![](https://www.r-pkg.org/badges/ago/Sim.DiffProc?color=blue)

------------------------------------------------------------------------

Total    | Monthly |  Weekly |  Daily
---------|---------|---------|-----
[![](https://cranlogs.r-pkg.org/badges/grand-total/Sim.DiffProc?color=yellowgreen)](https://cran.r-project.org/package=Sim.DiffProc) | [![](https://cranlogs.r-pkg.org/badges/Sim.DiffProc?color=yellowgreen)](https://cran.r-project.org/package=Sim.DiffProc) | [![](https://cranlogs.r-pkg.org/badges/last-week/Sim.DiffProc?color=yellowgreen)](https://cran.r-project.org/package=Sim.DiffProc) | [![](https://cranlogs.r-pkg.org/badges/last-day/Sim.DiffProc?color=yellowgreen)](https://cran.r-project.org/package=Sim.DiffProc)

------------------------------------------------------------------------


Package Overview
---------------------

The package [Sim.DiffProc](https://doi.org/10.18637/jss.v096.i02) is an object created in R for symbolic and numerical computations on scalar and multivariate systems of stochastic differential equations. It provides users with a wide range of tools to simulate, estimate, analyze, and visualize the dynamics of these systems in both forms Ito and Stratonovich. The project was officially launched in September 2010 and is under active development by the authors. The current feature set of the package can be split in more main categories: Computing the stochastic integrals of Ito or Stratonovich type. Simulation sde's and bridge sde's of Ito or Stratonovich type (1,2 and 3-dim), with different methods. Approximate transition density and random number generators for SDE's. Density approximation for First-passage-time (f.p.t) in SDE's (1,2 and 3-dim). Statistical analysis with Parallel Monte-Carlo and moment equations methods of SDE's (1,2 and 3-dim). Estimate drift and diffusion parameters using pseudo-maximum likelihood estimators of 1-dim SDE's. Displaying an object inheriting from a class of SDE's.

The package includes the following categories (where `k=1,2,3`):

1. [`snssdekd()` & `dsdekd()` & `rsdekd()`- Monte-Carlo Simulation and Analysis of Stochastic Differential Equations](https://CRAN.R-project.org/package=Sim.DiffProc/vignettes/snssde.html).
2. [`bridgesdekd()` & `dsdekd()` & `rsdekd()` - Constructs and Analysis of Bridges Stochastic Differential Equations](https://CRAN.R-project.org/package=Sim.DiffProc/vignettes/bridgesde.html).
3. [`fptsdekd()` & `dfptsdekd()` - Monte-Carlo Simulation and Kernel Density Estimation of First passage time](https://CRAN.R-project.org/package=Sim.DiffProc/vignettes/fptsde.html).
4. [`MCM.sde()` & `MEM.sde()` - Parallel Monte-Carlo and Moment Equations for SDEs](https://CRAN.R-project.org/package=Sim.DiffProc/vignettes/mcmsde.html).
5. [`TEX.sde()` - Converting Sim.DiffProc Objects to LaTeX](https://CRAN.R-project.org/package=Sim.DiffProc/vignettes/sdetotex.html).
6. [`fitsde()` - Parametric Estimation of 1-D Stochastic Differential Equation](https://CRAN.R-project.org/package=Sim.DiffProc/vignettes/fitsde.html).

Obtaining and installation
-----------------------

As `Sim.DiffProc` is an `R` package, it requires `R version 3.0.0` or higher to be installed, distributed as open source software under the GPL-2/GPL-3 license. The package is available from CRAN at URL https://CRAN.R-project.org/package=Sim.DiffProc, or from GitHub at URL https://github.com/acguidoum/Sim.DiffProc. To download, install and load the current release, just type the code below in your current `R` session:

```r
install.packages("Sim.DiffProc")
## Or 
install.packages("devtools")
devtools::install_github("acguidoum/Sim.DiffProc")
library("Sim.DiffProc")

-- by @runifsubset, based on A.C. ''Guidoum \and K. Boukhetala''' [https://cran.r-project.org/web/packages/Sim.DiffProc/vignettes/snssde.html#fn1]

\begin{equation}\label{eq:}
\begin{cases}
\begin{split}
\frac{d}{dt} m(t) &= 0.5 \, \left( \theta^2 \, m(t) \right) \\
\frac{d}{dt} S(t) &= \theta^2 \, \left( 2 \, S(t) + m(t)^2 \right)
\end{split}
\end{cases}
\end{equation}
