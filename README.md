
# Relativistic viscous hydrodynamics code for the simulation of QGP in heavy-ion collisions for 1 dimenson [![Build Status](https://travis-ci.com/fleer/HIC_Solver_1D.svg?token=sGDhmF9VD4p3i1s4BjQK&branch=master)](https://travis-ci.com/fleer/HIC_Solver_1D)

This program simulates the space-time evolution of the Quark-Gluon-Plasma in one spacial dimension, using numerical relativistic viscous hydrodynamics. 
The algorithm is based on [1], which is a modified version of the *Two-Shock Riemann Solver* presented in [2].
It utilizes the second-order Godunov method and is able to solve different kind of Riemann problems.

A Riemann problem is an initial-value (\\(t=0\\)) problem of the kind
\\[
\mathbf{A}(x,0) = \left\lbrace \begin{array}{l r}
\mathbf{L} & \textrm{für} \; x <0 \\ & \\ \mathbf{R} & \textrm{für} \; x >0
\end{array}
\right. ,
\\]
where \\(\mathbf{L}\\) and \\(\mathbf{R}\\) characterize the left and right initial states that are separeted by a singularity.






## Build

### Building on Linux

Install the **GNU Scientific Library (GSL)**, eg.
     
```
 wget ftp://ftp.gnu.org/gnu/gsl/gsl-2.4.tar.gz;
 tar -xzvf gsl-2.4.tar.gz;
 cd gsl-2.4 && ./configure && make && sudo make install
```
Clone the repository and build the code via make

```
git clone https://github.com/fleer/HIC_Solver_1D.git
cd HIC_Solver_1D && make
```

### Building on OS X

Install [Homebrew](https://brew.sh/index_de.html) via:

```
/usr/bin/ruby -e "$(curl -fsSL https://raw.githubusercontent.com/Homebrew/install/master/install)"
```

Then install **gcc** (tested for gcc-5, gcc-6, gcc-7) with openmp support and the **GNU Scientific Library (GSL)**.

```
brew install gcc --without-multilib
brew install gsl
```

Clone the repository and build the code via make

```
git clone https://github.com/fleer/HIC_Solver_1D.git
cd HIC_Solver_1D && make CXX=g++-7
```


## Simulate

**At the current point, you have to create the ditectory 'DATA' where the results are written manually!**

The parameters for the simulation are stored in the `input.ini`. The program is able to solve simple Riemann Problems as the shocktube problem, but also more complicated problems.

The initial states \\(\mathbf{L}_1, \mathbf{L}_2, \mathbf{R}_1, \mathbf{R}_2\\) are characterized through their initial baryon density \\(n_{\textrm{B}}\\), velocity \\(v\\) and pressure \\(p\\).

The code includes two different equations of state:

1. Approximated EoS of QCD, introduced in [3]
2. Eos of a free Gluon Gas [1][4]

Configuration files for the different types of problems are in the `Examples` folder, together with [gnuplot](http://www.gnuplot.info/) scripts for plotting *contourplots* and generating *gif-animations*.

### Shocktube Problem

![Example](./images/shock_tube_problem.gif)


### Landau Model

![Example](./images/landau.gif)


#### Landau Model with small Pertubation
---
[1] Akamatsu, Yukinao, et al. "A new scheme of causal viscous hydrodynamics for relativistic heavy-ion collisions: A Riemann solver for quark–gluon plasma." Journal of Computational Physics 256 (2014): 34-54.

[2] A. Mignone, T. Plewa and G. Bodo, arXiv.org (2005), astro-ph/0505200v1.

[3] Z. Fodor et al., Journal of High Energy Physics 2010 (2010) 77.

[4] E. Molnar, H. Niemi and D.H. Rischke, The European Physical Journal C 65 (2009) 615.
