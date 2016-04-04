*********************
CTDE API
*********************

Physical State
---------------------

The physical state is composed of a set of disjoint substates. Each substate can be read and written with a key. That means anything you want to be state should have a `getindex` and `setindex!`.
A simple example is `state=zeros(Int, 10)`.
A more complicated example, this state is a count of bees at different plants.

.. code-block:: julia

    type PhysicalState
        plants::Array{Int,2}
        bees::Array{Int,2}
    end
    
    function getindex(state::PhysicalState, x, y)
        state.bees[x, y]
    end
    
    function setindex!(state::PhysicalState, value, x, y)
    	state.bees[x, y]=value
    end

When specifying the process, intensities and firing functions can depend on the state only through a certain set of indexes, which you specify. In this example, the indices will be of the form of a tuple, `(x, y)` of type `Tuple{Int, Int}`.

Partial Process
----------------------
The library defines the partial process, `CTDE.PartialProcess`, which is used to specify the simulation and compute trajectories.

.. code-block:: julia

    process=PartialProcess(physical_state)

Then specify the process by adding transitions.


Transitions
-----------------

Each transition has two logical parts, an intensity function which describes when the transition is enabled and how long after enabling it will fire, and a firing function which describes what the transition does to the state when it fires.

.. function:: AddTransition!(process::PartialProcess, hazard::Intensity, hazard_dependencies, firing::Function, firing dependencies, name, keywords)

  * `process::PartialProcess` is the partial process just defined above.
  * `hazard::Intensity` is a hazard rate, derived from the abstract class Intensity, as will be described below.
  * `hazard_dependencies` is a list of indices into the physical state that will be passed to the intensity for this transition.
  * `firing::Function` is a function that changes the state, as described below.
  * `firing_dependencies` is a list of indices into the physical state that will be passed to the intensity for this transition.
  * `name` is a friendly name to use when the process tells you what transition fired.
  * `keywords` means that you can pass values such as `index=3` to the transition as messages to the samplers. Some samplers need hints about how to organize their work.
  * Returns nothing.


Intensity
------------------
An intensity says when the transition is enabled or disabled and, when it is enabled, the distribution of times at which it may fire.

.. class:: Intensity

    An *abstract* class for transition intensities. There are methods
    defined on this abstract class which help implement intensities.

.. function:: Enabled(intensity::Intensity)

    Returns a boolean to say whether the intensity is currently enabled.
    The helper method returns `intensity.enabled` for the object passed.

.. function:: Reset!(intensity::RecoverIntensity, time::Float64, state, keys..)

    When a transition fires, this is called to tell the intensity that it must forget all past observations of the state and determine, from the state at the values specified by the `keys`, what is the new distribution going forward.

.. function:: Update!(intensity::RecoverIntensity, time, state, keys...)

    This is the workhorse of the intensity distribution. Given the state at the given set of keys, the intensity chooses its current distribution for the hazard rate. It returns a symbol to report what happened. That symbol is either `:Unmodified`, `:Disabled`, `:Enabled`, or `:Modified`. The last choice, `:Modified`, means that the hazard was nonzero and is now nonzero but with a different distribution.

.. function:: Sample(intensity::Intensity, when::Float64, rng::MersenneTwister)

    samples the current distribution for the hazard, given that it has not yet fired by time `when`. This method calls `Sample(intensity.distribution)`, so a type which defines a `distribution` member doesn't need to reimpliment this method.

.. function:: Putative(intensity::Intensity, when::Float64, exponential_interval::Float64)

    integrates the current distribution for the hazard to determine at what time it will have used up an integrated hazard equal to `exponential_interval`. This is a way to sample distributions for Gibson and Bruck's Next Reaction Method or Anderson's method. This method calls `Putative(intensity.distribution)`, so a type which defines a `distribution` member doesn't need to reimpliment this method.

Implementing an intensity is simpler, thanks to the helper methods.

.. code-block:: julia

    type InfectIntensity <: Intensity
        distribution::TransitionDistribution
        enabled::Bool
        InfectIntensity(dist)=new(dist, 0.0, false)
    end
    
    function Reset!(intensity::InfectIntensity, time, state, who, whom)
        distribution.enabling_time=time
        Update!(intensity, time, state, who, whom)
    end
    
    function Update!(intensity::InfectIntensity, time, state, who, whom)
        modified=:Undefined
        enabled=(state[who]==1 && state[whom]==0)
        if enabled!=intensity.enabled
            if enabled
                intensity.distribution.enabling_time=time
                modified=:Enabled
            else
                modified=:Disabled
            end
            intensity.enabled=enabled
        else
            modified=:Unmodified
        end
        modified
    end

In general, the intensity can depend on any state since it last fired or the start of the simulation. In practice, an intensity will examine the state to create parameters for the distribution. The `WrappedDistribution` is a good example of the interface distributions support.

Transition Distributions
--------------------------
The distributions an intensity needs have different methods from
distributions in Julia's `Distributions` module.

.. class:: TransitionDistribution

    This **abstract** class is a base class for the continuous univariate distributions used by intensities.

WrappedDistribution
^^^^^^^^^^^^^^^^^^^^^^^^^
.. class:: WrappedDistribution

   This class uses the available distributions in the `Distributions` package to meet the API needed by the simulation. It's likely less efficient and possibly numerically inaccurate for some distributions, but here goes. Its members are `relative_distribution` and `enabling_time`.

.. function:: WrappedDistribution(dist::Distributions.ContinuousUnivariateDistribution, enabling_time::Float64)

    The constructor. Pass in a distribution from the Distributions package.

.. function:: Sample(distribution::WrappedDistribution, now::Float64,
        rng::MersenneTwister)

   This samples the distribution using the given random number generator. It calls `quantile(distribution, now, rand(rng))`.

.. function:: HazardIntegral(dist::WrappedDistribution, t1, t2)

    This integrates the hazard from time `t1` to time `t2` using
    `logccdf(dist, t1-te)-logccdf(dist, t2-te)` where `te` is the enabling time.

.. function:: ImplicitHazardIntegral(dist::WrappedDistribution, cumulative_hazard::Float64, when::Float64)

    This is the inverse of the hazard integral, so
    `ImplicitHazardIntegral(dist, HazardIntegral(dist, t1, t2), t1)=t2`.

.. function:: EnablingTime(dist::WrappedDistribution)

   Return the enabling time. It's a common parameter of all of these distributions.

.. function:: EnablingTime!(dist::WrappedDistribution, time::Float64)

   Set the enabling time.

.. function:: Parameters(dist::WrappedDistribution)

    This returns the parameters for the distribution. The exact set depends on the underlying distribution.

Exponential
^^^^^^^^^^^^^^^^^

.. function:: TransitionExponential(λ::Real, enabling_time::Real)

   The rate, λ, is the hazard rate, not a scale.

   .. math:: F(T)=1-exp\left[-\lambda(T-Te)\right]

Weibull
^^^^^^^^^^^^^^
.. function:: TransitionWeibull(θ::Float64, k::Float64)

    The scale is :math:`\theta`. :math:`k` is the exponent.

    .. math:: F(T)=1-exp\left[-\left(\frac{T-Te}{\theta}\right)^k\right]

Gamma
^^^^^^^^^^^^^^^^
.. function:: TransitionGamma(α, β, enabling_time)

    α is the shape parameter, β the inverse scale parameter, also called a rate parameter.

    .. math:: f(t)=\frac{\beta^\alpha}{\Gamma(\alpha)}x^{\alpha-1} e^{-\beta x}

LogLogistic
^^^^^^^^^^^^^^^^^^^^
.. function:: TransitionLogLogistic(α, β)

   .. math:: F(t)=\left[1 + \left(\frac{t-t_e}{\alpha}\right)^{-\beta}\right]^{-1}

NelsonAalen
^^^^^^^^^^^^^^^^
.. class:: NelsonAalenDistribution

    Given a list of times the distribution either fired or was right-censored, meaning it failed to fire, this constructs an estimator of the distribution that can be sampled.

    Because each transition in the process will either fire or be interrupted, this estimator can be used to ask whether the firing of each transition matches the expected distribution.

.. function:: cdf(dist::NelsonAalenDistribution, bypoint::Int)

Empirical
^^^^^^^^^^^^^
.. class:: EmpiricalDistribution

    This class estimates a distribution given samples of the times at which that distribution fired. First make the object and then use `push!` to add values. Finally, call `build!` before sampling from it.
    This is useful for testing distributions.

.. function:: EmpiricalDistribution()

    Constructor.

.. function:: push!(ed::EmpiricalDistribution, value::Float64)

   Add a sample to the list.

.. function:: build!(ed::EmpiricalDistribution)

   Internally, it needs to sort the list of samples.

.. function:: mean(ed::EmpiricalDistribution)

.. function:: variance(ed::EmpiricalDistribution)


.. function:: kolmogorov_smirnov_statistic(ed::EmpiricalDistribution, other)

    The `other` is a distrubition for which `cdf` is defined. This returns two values, the maximum difference between the two distributions and whether that maximum difference meets the 0.05 confidence interval for the hypothesis that they are the same distribution.



Firing Function
-----------------

The firing function is a function that modifies state. Its signature has to be

.. function:: FiringFunction(state, keys...)

the firing function returns a list of all substates which could be affected by having fired. While the list of hazard dependencies and firing dependencies is state, the list of what was affected by firing is not. As an example, recovery and infection for a disease model could look like
If the firing function reads or writes mutable state other than that specified by the indices in keys, then the simulation will be incorrect.

.. code-block:: julia

    function Recover!(state, who)
        state[who]=0
        [who]
    end
    
    function Infect!(state, who)
        state[who]=1
        [who]
    end

Simulation Observer
--------------------
Every time the simulation determines the next time and transition, it changes the state and then calls an observer with the signature

.. function:: StateObserver(state, affected_keys, clock_name, time::Float64)

    * `state` is the state of the system.
    * `affected_keys` are those substates which were affected when the last transition fired.
    * `clock_name` is the name given to that transition.
    * `time` is the time at which that transition happened.
    * Returns: The observer returns a boolean indicating whether the simulation may continue.

In practice, this function is a closure which adds data to a list of data, or writes that list to disk. For instance,

.. code-block:: julia

    function Observer(out::ScreenObserver)
        function sobserve(state::Array{Int,1}, affected_keys, clock_name, time::Float64))
            AddEntry!(out, state[affected_keys[1], time])
        end
    end

Sampler
--------

Define a sampler. There are a few to choose among.

.. class:: NextReactionHazards

    This sampler uses a variant of Gibson and Bruck's next reaction algorithm, described by Anderson and Kuntz.

.. function:: NextReactionHazards()

    Constructor.

.. class:: FixedDirect

    This is an optimized direct reaction sampler. It assumes there are a fixed number of transitions in the system and that every transition is given an ordinal with keywords of the form `index=<Int>`.

.. function:: FixedDirect(N::Int)

    Constructor. `N` is the number of transitions in the system.

.. class:: NaiveSampler

    This is appropriate only for simulations where no transition is ever reenabled after it fires or is disabled. It is equivalent to Next Reaction in this case.

Running a Simulation
----------------------

Once the process is created and sampler chosen, a single function runs the simulation. In this example, the `MakeProcess` function creates state and adds transitions to the process.

.. code-block:: julia

    rng=MersenneTwister(333333)
    N=3
    parameters=Dict(:Gamma =>1.0, :Beta => 1.0)
    process, state=MakeProcess(N, parameters, rng)
    observer=SamplingObserver(N, 1000)
    sampler=NextReactionHazards()

    RunSimulation(process, sampler, Observer(observer), rng)

