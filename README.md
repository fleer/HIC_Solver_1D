
# Relativistic viscous hydrodynamics code for the simulation of QGP in heavy-ion collisions for 1 dimenson [![Build Status](https://travis-ci.com/fleer/HIC_Solver_1D.svg?token=sGDhmF9VD4p3i1s4BjQK&branch=master)](https://travis-ci.com/fleer/HIC_Solver_1D)

This algorithm, based on [1], simulates the space-time evolution of Quark-Gluon Plasma, created during heavy-ion collisions, using relativistic viscous hydrodynamic.
It utilizes the second-order Godunov method and is able to solve different test problems.

[1] Akamatsu, Yukinao, et al. "A new scheme of causal viscous hydrodynamics for relativistic heavy-ion collisions: A Riemann solver for quark–gluon plasma." Journal of Computational Physics 256 (2014): 34-54.



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


## Simulate

**At the current point, you have to create the ditectory 'DATA' where the results are written manually!**

### Shocktube Problem


### Landau Model

![Example](./images/animate.gif)