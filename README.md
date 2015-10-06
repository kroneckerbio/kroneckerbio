KroneckerBio
============

KroneckerBio is a systems biology modeling toolbox for Matlab. It provides an easy-to-use programming interface for building, simulating, and analyzing ODE models of biological systems. The toolbox incorporates numerous methods developed in the Tidor lab at MIT.

For mass action models, simulations can be run using sparse matrix operations, a fact that KroneckerBio exploits. Because Matlab has fast sparse matrix algorithms, simulating mass action models and running analyses that are dependent on simulations in KroneckerBio is very fast.

Included in KroneckerBio is a rich set of methods for quantifying the uncertainty in the parameters and topology of a model and doing optimal experimental design to predict which experiments would be best for reducing the remaining uncertainty.

Installation
------------

KroneckerBio is written entirely in Matlab. No compilation is necessary. Simply download the entire source of the [most recent stable release](https://github.com/kroneckerbio/kroneckerbio/archive/stable.zip). Once the toolbox is downloaded, run the `InitKronecker.m` Matlab script. This script simply adds several KroneckerBio folders to the Matlab path so that the functions are available.

At the moment, KroneckerBio has a brief tutorial that covers the basics of model building and import, simulation, and fitting. The `T0x` scripts in the `Testing` folder cover a broad range of common uses of toolbox. The help sections of each Matlab function should be mostly up-to-date.

In addition to the tutorial, email the lead developer [David Hagen](`david@drhagen.com`) for help in using KroneckerBio. Bugs can also be reported to him or in the GitHub issues tab.
