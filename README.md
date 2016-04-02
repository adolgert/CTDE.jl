CTDE
===========

Continuous-time, stochastic discrete event simulation.
An exact stochastic simulation algorithm for time-dependent hazard rates.

* **Stochastic**&mdash;Uses a random number generator to create a sequence of events whose probability distribution is well-defined but whose values change each time it is run.
* **Continuous time**&mdash;There isn't a clock-tick to the simulation, so events happen at the rate with which they are prescribed to happen. No two events happen at the same time. No events happen immediately after another event.
* **Gillespie-type algorithm**&mdash;Rates for transitions are specified by hazard rates, also known as propensities.
* **Non-constant hazard rates**&mdash;Most Gillespie-type simulations specify each reaction with an exponential distribution, but these simulations can use Weibull, Gamma, or other distributions, including estimators of distributions measured in the field. They still fire correctly, according to Gillespie's algorithm.

[Documentation](https://ctdejl.readthedocs.org/)

Drew Dolgert, ajd27@cornell.edu
Contact me with any questions.

[![Build Status](https://travis-ci.org/adolgert/CTDE.jl.svg?branch=master)](https://travis-ci.org/adolgert/CTDE.jl)
