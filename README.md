CTDE
===========

Continuous time, discrete event system in Julia

Drew Dolgert, ajd27@cornell.edu

This is an implementation, in the Julia language, of a
generalized semi-Markov Petri net (GSPN). It is suitable for general
discrete-event modeling in continuous time with non-Exponential
distributions. This is similar to the Semi-Markov library, written
in C++, but this implementation is actually a more formal GSPN, because
it obeys stoichiometry exactly.

The code is alpha. Features include

* GSPN implementation
* Gillespie's first reaction method
* Gibson and Bruck's next reaction method (via Anderson)

Contact me with any questions.

[![Build Status](https://travis-ci.org/adolgert/CTDE.jl.svg?branch=master)](https://travis-ci.org/adolgert/CTDE.jl)
