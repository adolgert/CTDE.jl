
*********************
Examples
*********************

* `Well-mixed SIR <https://github.com/adolgert/CTDE.jl/tree/master/example/sir>`_ as a test of Ball and Nasell's predictions of SIR final size distributions.
* `Bouncing Rabbits <https://github.com/adolgert/CTDE.jl/blob/master/example/rabbits.jl>`_ which move on a two-dimensional board, without stepping on each other, but they infect each other, sadly.
* `Susceptible-Infectious-Susceptible <https://github.com/adolgert/CTDE.jl/blob/master/example/sis.jl>`_ This is a good test of the samplers because results are well-known. It plots the master probabilities for a regenerative process.
* `Object-based Simulation <https://github.com/adolgert/CTDE.jl/blob/master/example/sir.jl>`_ This doesn't use dictionary keys in order to access substates of the physical state. It uses the objects themselves, which lends it to more free-form simulation.
