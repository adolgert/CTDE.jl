******************************************
Detailed Description of Distributions
******************************************
.. |nbsp| unicode:: 0xA0 
   :trim:

This describes implementation of distributions for sampling of
distributions for stochastic simulation in continuous time. This kind of
simulation makes specific demands on what calls a distribution must
support, and those calls are different from what libraries provide.
This is a guide for implementation of new distributions
and a way to ensure that those implemented look correct.
If something is wrong here, it matters, so file a bug report.

Notation
===========
First, let's affix notation. The cumulative distribution function
of every regular distribution can be written as an integral over
its hazard rate, :math:`\lambda`

.. math:: F(t)=1-e^{-\int_{0}^t \lambda(s)ds}.
   :label: basic-distribution

All algorithms for stochastic simulation treat distributions as
being defined in absolute time, specified as an enabling time,
:math:`t_e`,

.. math:: F(t, t_e)=1-e^{-\int_{0}^{t-t_e} \lambda(s)ds}.
   :label: absolute-distribution

Working with distributions in absolute time is a simple shift of the
time scale and will be ignored in further discussions, although
the enabling time, :math:`t_e`, will certainly appear in code.
Fig. |nbsp| :eq:`basic-distribution`.

The density function is the derivative of the cumulative distribution
function,

.. math:: f(t)=\frac{dF(t)}{dt}=\lambda(t)e^{-\int_{0}^t \lambda(s)ds}.
   :label: density-function

The survival is

.. math:: G(t)=1-F(t)=e^{-\int_{0}^t \lambda(s)ds}.
   :label: survival

Because survival is multiplicative, we further label the survival
from time :math:`t_0` to :math:`t_1` as

.. math:: G(t_0, t_1)=\frac{G(t_1)}{G(t_0)}=e^{-\int_{t_0}^{t_1} \lambda(s)ds}
   :label: survival-fromto

Requirements for a Continuous-Time Simulation
==============================================

Shifted Sample
-----------------
The First Reaction method requires that we sample a distribution
given that we known it has not yet fired by a time :math:`t_0`.
The statement that it hasn't fired by time :math:`t_0` creates a
new distribution from which to sample. If the old distribution
had the hazard :math:`G(t)=G(0, t)`, it could be written as

.. math:: G(0, t)=G(0, t_0)G(t_0, t).
   :label: survival-relationship

It is this partial survival, since :math:`t_0`, that we want to sample.
Solving for :math:`G(t_0, t)` and subtracting both sides from 1,

.. math:: 1-G(t_0, t)=\frac{G(0, t_0)-G(0, t)}{G(0, t_0)}.

Written in terms of the cumulative distribution functions, the
cdf of the new distribution, which we'll call :math:`F(t, t_0)`, is

.. math:: F(t, t_0)=\frac{F(t)-F(t_0)}{1-F(t_0)}
   :label: gibson-shifted

This kind of distribution could be sampled by a rejection method,
but the default way to sample it is by inversion, which
means generating a uniform random value between :math:`[0,1]`
and solving :math:`U=F(t)` for :math:`t`. For Eq. |nbsp| :eq:`gibson-shifted`,
this becomes

.. math::
   :label: sample-shifted

   \begin{eqnarray}
    U&=&F(t,t_0,t_e) \\
     &=&\frac{F(t,t_e)-F(t_0,t_e)}{1-F(t_0,t_e)} \\
   U(1-F(t_0,t_e))&=&F(t,t_e)-F(t_0,t_e) \\
   F(t,t_e)&=&U(1-F(t_0,t_e))+F(t_0,t_e) \\
   F(t-t_e)&=&U(1-F(t_0-t_e))+F(t_0-t_e) \\
   t-t_e &=& F^{-1}\left[U(1-F(t_0-t_e))+F(t_0-t_e)\right] \\
   t &=& t_e+F^{-1}\left[U(1-F(t_0-t_e))+F(t_0-t_e)\right]
   \end{eqnarray}

We will call this operation **SampleShifted.**

Hazard Rate for Next Reaction
--------------------------------
The Next Reaction method requires sampling a distribution such
that the quantile is saved, so that later adjustments to the
distribution can use the same quantile.

During a simulation, the hazard rate, :math:`\lambda`, is a function
of the state of the system, :math:`X(t)`. The state of the system
only changes in jumps, so the hazard rate is effectively
a series of distributions in time. For instance, a hazard rate,
from enabling time :math:`T_0` to firing time :math:`T_3`, might have three parts.

.. math::
   :label: hazard-steps

   \begin{eqnarray}
     \lambda(\{X_0, T_0\}, t)=h_0(t) & &\qquad T_0 \le t < T_1 \\
     \lambda(\{X_0, T_0, X_1, T_1\}, t)=h_1(t) & &\qquad T_1 \le t < T_2 \\
     \lambda(\{X_0, T_0, X_1, T_1, X_2, T_2\}, t)=h_2(t) & &\qquad T_2 \le t < T_3
   \end{eqnarray}

The algorithm therefore samples for a firing time from :math:`h_0(t)`
as soon as the transition is enabled, but that time will turn out to 
be wrong (we call it a putative time). Later, the algorithm will
resample using :math:`h_1(t)` using the original sample's quantile
and taking time off the clock.
If the first sample were by inversion, it would look like solving this
equation for :math:`t` (still ignoring enabling times),

.. math::

    U=1-\exp\left(\int_0^{t}h_0(s)ds\right).

Then a later sample would use the same :math:`U`, but with knowledge that
the distribution now contains a new part, :math:`h_1(t)`,

.. math::
   :label: inversion-bite

    U=1-\exp\left(-\int_0^{t_1}h_0(s)ds\right)\exp\left(-\int_{t_1}^{t}h_1(s)ds\right)

Anderson had the bright idea to write the quantile as an
exponential quantile,
:math:`1-U=e^{\ln (1-U)}`, so that the equation requires only
addition of integrated hazards,

.. math::
  :label: putative-implicit

  \begin{eqnarray}
    \ln(1-U)&=&-\int_0^{t_1}h_0(s)ds-\int_{t_1}^{t}h_1(s)ds \\
    \int_{t_1}^{t}h_1(s)ds&=&-\ln(1-U)-\int_0^{t_1}h_0(s)ds.
  \end{eqnarray}

As the underlying distribution, :math:`h_i(t)`, changes, the
right hand side gets smaller and smaller. Let's call the sum
of all consumed hazard :math:`\gamma`,

.. math::

   \gamma=\sum_i \int_{t_i}^{t_{i+1}}h_i(s)ds

The algorithm therefore needs three operations from the distribution.

1.  **MeasuredSample**---Sample the distribution, returning the exponential quantile.
    Calling the random number generator, "rng," and the putative
    time :math:`t_p`, it's

    .. math:: \left(\mbox{rng}, h_0(t)\right) \mapsto \left(t_p, -\ln(1-U)\right).

2.  **ConsumeSample**---Consume remaining quantile for the next sample. If the sum of
    hazard in the past is called :math:`\gamma`, then

    .. math:: \left(\gamma, h_i(t), t_i\right) \mapsto \left(\gamma', t_{i+1}\right)

3.  **Putative**---Generate a new putative time from the exponential quantile and
    the consumed hazard,

    .. math:: \left(-\ln(1-U), \gamma, t_i, h_i(t)\right) \mapsto p_t

The nice part about the first step is that there is no need
to sample by inversion. Any sampling method will do, as long
as the exponential quantile is calculated.



Cumulative Distributions for Next Reaction
-------------------------------------------
The original form of the Next Reaction, by Gibson and Bruck,
was written in terms, not of the hazards, but of the
cumulative distribution functions. This form remains useful because
some distributions are much simpler, or more accurate, to sample
as cdfs instead of sampling from their hazard rates.

Returning to Eq. |nbsp| :eq:`inversion-bite`, this can be rewritten
as 

.. math::

    1-U=\exp\left(-\int_0^{t_1}h_0(s)ds\right)\exp\left(-\int_{t_1}^{t}h_1(s)ds\right)

In terms of the survival functions, this becomes

.. math::

    1-U=G_0(0, t_1)G_1(t_1, t)

If we wish to solve this for :math:`t`, then, in terms of the survival,
it looks like

.. math::

   G_1(t_1, t)=\frac{1-U}{G_0(0, t_1)}

Writing the left-hand side as a cumulative distribution function
requires the transformation

.. math::

   G(a, b)=\frac{G(b)}{G(a)}=\frac{1-F(b)}{G(a)}

so we have

.. math::

   F_1(t)=1-G_1(t_1) \frac{1-U}{G_0(0, t_1)}

This generalizes to many steps as

.. math::

   F_j(t)=1-G_j(t_j) (1-U) \prod_i^j \frac{G_i(t_i)}{G_i(t_{i+1})}

Let's call the running product on the right :math:`\delta`,

.. math::

   \delta=\frac{G_i(t_i)}{G_i(t_{i+1})}

Then the algorithm requires three operations

1.  **MeasuredSample**---Sample the distribution, returning the quantile.
    Calling the random number generator, "rng," and the putative
    time :math:`t_p`, it's

    .. math:: \left(\mbox{rng}, h_0(t)\right) \mapsto \left(t_p, 1-U\right).

2.  **ConsumeSample**---Consume remaining quantile for the next sample. 

    .. math:: \left(\delta_i, h_i(t), t_i\right) \mapsto \left(\delta_{i+1}, t_{i+1}\right)

3.  **Putative**---Generate a new putative time from the exponential quantile and
    the consumed hazard,

    .. math:: \left(1-U, \delta_i, t_i, h_i(t)\right) \mapsto p_t

As you can see by comparison with the hazards version, it's simple
to write the algorithm to accommodate either method of sampling.
Therefore, each distribution can choose which interface to support.


Using Julia's Distributions
==============================

Julia's continuous univariate distributions support a fixed
interface. In this section, we look at how to translate any
distribution into the operations above.

In this table, ``d`` is the distribution.

+------------------------+-------------------------------------------------------------------------+
| ``cdf(d,t)``           | :math:`F(t)`                                                            |
+------------------------+-------------------------------------------------------------------------+
| ``quantile(d,q)``      | :math:`F^{-1}(q)`                                                       |
+------------------------+-------------------------------------------------------------------------+
| ``logcdf(d,t)``        | :math:`\ln(F(t))`                                                       |
+------------------------+-------------------------------------------------------------------------+
| ``ccdf(d,t)``          | :math:`G(t)`                                                            |
+------------------------+-------------------------------------------------------------------------+
| ``logccdf(d,t)``       | :math:`-\int_0^t \lambda(s)ds`                                          |
+------------------------+-------------------------------------------------------------------------+
| ``quantile(d,q)``      | :math:`F^{-1}(q)`                                                       |
+------------------------+-------------------------------------------------------------------------+
| ``cquantile(d,q)``     | :math:`F^{-1}(1-q)=G^{-1}(q)`                                           |
+------------------------+-------------------------------------------------------------------------+
| ``invlogcdf(d,lp)``    | :math:`F^{-1}(e^{l_p})`                                                 |
+------------------------+-------------------------------------------------------------------------+
| ``invlogccdf(d,lp)``   | :math:`G^{-1}(e^{l_p})` or :math:`-\int_0^{t(l_p)}\lambda(s)ds=l_p`     |
+------------------------+-------------------------------------------------------------------------+
| ``randexp(rng)``       | :math:`-\ln(1-U)`                                                       |
+------------------------+-------------------------------------------------------------------------+

A shifted sample, from Eq. |nbsp| :eq:`sample-shifted`, which
ends with

.. math::

   t = t_e+F^{-1}\left[U(1-F(t_0-t_e))+F(t_0-t_e)\right]

transliterates to

.. literalinclude:: ../../src/wrappeddistribution.jl
   :language: julia
   :lines: 58-64
   :name: sample-shifted-wrapped-code
   :caption: wrappeddistribution.jl


The next two pieces concern the hazard. The goal is to find the integral
of the hazard between two absolute times, :math:`t_1` and :math:`t_2`,
where both are :math:`t_{1,2}\ge t_0`. This is

.. math::

   \int_{t_1-t_e}^{t_2-t_e} \lambda(s)ds=\int_{0}^{t_2-t_e} \lambda(s)ds
       -\int_{0}^{t_1-t_e} \lambda(s)ds.

In terms of the given methods, this would be, noting the minus sign in
the table,

.. literalinclude:: ../../src/wrappeddistribution.jl
   :language: julia
   :lines: 79-85, 93-98
   :name: consume-wrapped-code
   :caption: wrappeddistribution.jl

Looking back to Eq. |nbsp| :eq:`putative-implicit`,

.. math::

  \begin{equation}
    \int_{t_1}^{t}h_1(s)ds=-\ln(1-U)-\int_0^{t_1}h_0(s)ds,
  \end{equation}

we can label ``xa`` the combination of the exponential quantile
and the sums of integrals on the right-hand side.

.. literalinclude:: ../../src/wrappeddistribution.jl
   :language: julia
   :lines: 101-105, 111-117
   :name: putative-wrapped-code
   :caption: wrappeddistribution.jl


Exponential
==============================================

The exponential distribution is constructed with a hazard
rate, even though the internal distributions object
uses a scale, which is :math:`\theta =1/\lambda`,

.. literalinclude:: ../../src/exponentialdistribution.jl
   :language: julia
   :lines: 13-16
   :name: exponential-constructor-code
   :caption: exponentialdistribution.jl

It doesn't matter how we sample the distribution, as long
as we return its quantile. This samples using
``Base.randexp``, which uses the ziggurat method for a sample
that's much faster than inversion. The value returned by
``randexp`` is equivalent to :math:`-\ln(1-U)`.

.. literalinclude:: ../../src/exponentialdistribution.jl
   :language: julia
   :lines: 47-51
   :name: exponential-measured-code
   :caption: exponentialdistribution.jl

The hazard integral for constant hazards is :math:`(t_2-t_1)\lambda`.

.. literalinclude:: ../../src/exponentialdistribution.jl
   :language: julia
   :lines: 53-57, 59-64
   :name: exponential-consume-code
   :caption: exponentialdistribution.jl

Even inverting the hazard integral is an increment with a multiplication.

.. literalinclude:: ../../src/exponentialdistribution.jl
   :language: julia
   :lines: 71-76, 78-81
   :name: exponential-putative-code
   :caption: exponentialdistribution.jl

Weibull
==========
Like the exponential distribution, the Weibull distribution
has an integrable hazard rate, which makes implementation
straightforward. Unfortunately, the use of the parameter
:math:`\lambda` in the definition of the Weibull is at odds
with our use of it as a hazard rate, but it's just a scale parameter
here.

.. math::
   :label: weibull-cdf

    F(t)=1-\exp\left[\left(\frac{t-t_e}{\lambda}\right)^k\right]

The constructor uses this cdf.

.. literalinclude:: ../../src/weibulldistribution.jl
   :language: julia
   :lines: 17-19
   :name: weibull-constructor-code
   :caption: weibulldistribution.jl

From the cdf, the hazard rate is

.. math::

   \Lambda(t)=\int_0^t\lambda(s)ds=\left(\frac{t-t_e}{\lambda}\right)^k

The inverse, where we ask when the integral equals :math:`l_u=-\ln(1-U)`,
is

.. math::

   t=t_e+ \lambda l_u^(1/k)

The version in the code is overachieving because it allows for shifting
the distribution.

.. literalinclude:: ../../src/weibulldistribution.jl
   :language: julia
   :lines: 49-60
   :name: weibull-sample-code
   :caption: weibulldistribution.jl

Given that the hazard is already integrated in Eq. |nbsp| :eq:`weibull-cdf`,
integrating the hazard is algebraic.

.. literalinclude:: ../../src/weibulldistribution.jl
   :language: julia
   :lines: 62-70, 72-77
   :name: weibull-consume-code
   :caption: weibulldistribution.jl

.. literalinclude:: ../../src/weibulldistribution.jl
   :language: julia
   :lines: 85-95, 97-100
   :name: weibull-putative-code
   :caption: weibulldistribution.jl


Log-Logistic
============

Working from wikipedia, because Gradstein and Ryzhik is too heavy to lift.

.. math:: F(x;\alpha, \beta)=\frac{1}{1+(x/\alpha)^{-\beta}}.

We shift this to

.. math:: F(t, t_e)=\frac{1}{1+((t-t_e)/\alpha)^{-\beta}}.

The pdf is

.. math::

   f(x;\alpha, \beta)=\frac{(\beta/\alpha)(x/\alpha)^{\beta-1}}
     {(1+(x/\alpha)^\beta)^2}.

The quantile is

.. math:: F^{-1}(p; \alpha, \beta)=\alpha \left(\frac{p}{1-p}\right)^{1/\beta}.

Survival

.. math:: G(t)=1-F(t)=\frac{1}{1+(t/\alpha)^\beta}.

Hazard

.. math::

   \lambda(t)=\frac{f(t)}{G(t)}=\frac{(\beta/\alpha)(t/\alpha)^{\beta-1}}
     {1+(t/\alpha)^\beta}

Lastly, we need ``invlogccdf(d,lp)``, which is
:math:`G_d^{-1}(e^{l_p})`, or :math:`-\int_0^t(l_p)\lambda(s)ds=l_p`.

.. math::

   \begin{aligned}
     l_p&=&\ln(G(t)) \\
     e^{l_p}&=&G(t) \\
     e^{l_p}&=&\frac{1}{1+(t/\alpha)^\beta} \\
     e^{-l_p}&=&1+(t/\alpha)^\beta \\
     (t/\alpha)^\beta&=&  1-e^{-l_p}\\
     t/\alpha&=& (1-e^{-l_p})^{1/\beta}\\
      t&=&\alpha(1-e^{-l_p})^{1/\beta}\\\end{aligned}

Gamma
=====

We will define paramaters from the shape :math:`\alpha` and rate
:math:`\beta`.

.. math:: f(x)=\frac{\beta^\alpha}{\Gamma(\alpha)}x^{\alpha-1}e^{-\beta x}

where

.. math:: \Gamma(t)=\int_0^\infty x^{t-1}e^{-x}dx.

The CDF is

.. math:: F(x;\alpha,\beta)=\frac{\gamma(\alpha,\beta x)}{\Gamma(\alpha)}

where :math:`\gamma` is the (lower) incomplete gamma function,

.. math:: \gamma(x;\alpha)=\int_0^x t^{\alpha-1}e^{-t}dt

In our back pocket, from ``Boost::Math``, are :math:`\Gamma(x)`,
:math:`\ln(|\Gamma(x)|)`, digamma, which is

.. math:: \psi(x)=\frac{d}{dx}\ln(\Gamma(x))=\frac{\Gamma'(x)}{\Gamma(x)},

gamma ratio, which is :math:`\Gamma(a)/\Gamma(b)`, gamma delta ratio,
which is :math:`\Gamma(a)/\Gamma(a+\Delta)`, and the set of incomplete
gamma functions. In order, they are normalized lower incomplete,
normalized upper, incomplete full (non-normalized) lower incomplete, and
full (non-normalized) upper incomplete gamma functions.

.. math::

   \begin{aligned}
     \mbox{gamma\_p}(a,z)&=&\frac{\gamma(a,z)}{\Gamma(a)}=\frac{1}{\Gamma(a)}
        \int_0^zt^{a-1}e^{-t}dt \\
     \mbox{gamma\_q}(a,z)&=&\frac{\Gamma(a,z)}{\Gamma(a)}=\frac{1}{\Gamma(a)}
        \int_z^0t^{a-1}e^{-t}dt \\
     \mbox{tgamma\_lower}(a,z)&=&\gamma(a,z)=
        \int_0^zt^{a-1}e^{-t}dt \\
     \mbox{tgamma}(a,z)&=&\Gamma(a,z)=\frac{1}{\Gamma(a)}
        \int_z^0t^{a-1}e^{-t}dt \\\end{aligned}

There are a set of inverses of incomplete gamma functions and
derivatives of incomplete gamma functions. OK, back to what we need.

.. math::

   \begin{aligned}
     F(x;\alpha,\beta)&=&\mbox{gamma\_p}(\alpha, \beta x) \\
     F^{-1}(y;\alpha,\beta)&=&\mbox{gamma\_p\_inv}(\alpha, y)/\beta\end{aligned}

The hazard integral, in terms of the cdf, is

.. math::

   \begin{aligned}
    \int_{t_1-t_e}^{t_2-t_e}\lambda(s)ds&=&-\ln(1-F(t_2-t_e))+\ln(1-F(t_1-t_e)) \\
    &=& \ln\left[\frac{1-F(t_1-t_e)}{1-F(t_2-t_e)}\right].\end{aligned}

Can we simplify this into something provided?

.. math::

   \begin{aligned}
   \int_{t_1-t_e}^{t_2-t_e}\lambda(s)ds & = & \ln\left[\frac{1-\frac{\gamma(\alpha,\beta (t_1-t_e))}{\Gamma(\alpha)}}{1-\frac{\gamma(\alpha,\beta (t_2-t_e))}{\Gamma(\alpha)}}\right] \\
    & = & \ln\left[\frac{\Gamma(\alpha)-\gamma(\alpha,\beta (t_1-t_e))}
    {\Gamma(\alpha)-\gamma(\alpha,\beta (t_2-t_e))} \right] \\
   \gamma(\alpha,\beta (t_1-t_e)) & = & \int_0^{\beta(t_1-t_e)} t^{\alpha-1}e^{-t}dt\end{aligned}

It looks like we might do best just with

::

    Ga=tgamma(a)
    hazint(te, t1, t2)=log((Ga-tgamma_lower(a,b*(t1-te)))/
        (Ga-tgamma_lower(a,b*(t2-te))))

Our other goal for Gamma distributions is to get the inverse hazard.
This can be seen as two steps. First find the integral

.. math:: l_p=-x+\left[\int_0^{t0-t_e}\lambda(s)ds\right].

Then solve for :math:`t'` in

.. math:: l_p=-\int_0^{t'-t_e}\lambda(s)ds.

Or, we could write this as

.. math:: l_e =e^{-x}e^{-\int_0^{t0-t_e}\lambda(s)ds}=e^{-x}(1-F(t_0-t_e))

and

.. math:: l_e=e^{-\int_0^{t'-t_e}\lambda(s)ds}=1-F(t'-t_e).

All at once,

.. math::

   \begin{aligned}
     F(t'-t_e)&=&1-e^{-x}(1-F(t_0-t_e)) \\
    t'&=&t_e+F^{-1}\left(1-e^{-x}(1-F(t_0-t_e))\right). \\
    F(t_0-t_e)&=&\mbox{gamma\_p}(\alpha,\beta(t_0-t_e)) \\
    F^{-1}(y)&=&\mbox{gamma\_p\_inv}(\alpha, y)/\beta\end{aligned}

So here is our inverse hazard integral.

::

      quad=1-exp(-x)*(1-gamma_p(a,b*(t0-te)))
      tp=te + gamma_p_inv(a, quad)/b

Uniform Distribution
====================

Maybe this one will be easier. This distribution has two parameters, a
start time and an end time, :math:`t_a` and :math:`t_b`. The pdf is
constant, :math:`f(t)=1/(t_b-t_a)` between :math:`t_a\le t<t_b`. The CDF
is just the integral of that, :math:`F(t)=(t-t_a)/(t_b-t_a)`. The
integrated hazard will have nonzero cases for for
:math:`t_1<t_a<t_2<t_b`, :math:`t_1<t_a<t_b<t_2`,
:math:`t_a<t_1<t_2<t_b`, :math:`t_a<t_1<t_b<t_2`. It is zero for
:math:`t_1<t_2<t_a` and :math:`t_a<t_b<t_1<t_2`

.. math::

   \int_{t_1-t_e}^{t_2-t_e}\lambda(s)ds=
         \ln\left[\frac{1-F(t_1-t_e)}{1-F(t_2-t_e)}\right]

If :math:`t_a\le t_n-t_e<t_b`, then
:math:`F(t_n-t_e)=(t_n-t_e-t_a)/(t_b-t_a)`. Otherwise it is :math:`0` or
:math:`1`. It should never be the case that a uniform distribution does
not fire before :math:`t_b`. The hazard integral always sums over time
already past in the simulation. Nevertheless, it will be necessary to
check for overflow near :math:`t_b`, and it would help to keep the two
logs separated, instead of in the fraction.

What about the inverse of the hazard integral?
:math:`F^{-1}(x)=t_a+(t_b-t_a)x` Therefore, for :math:`t_a\le t_0-t_e`,

.. math:: t'=t_e+t_a+(t_b-t_a)\left[1-e^{-x}\left(1-\frac{t_0-t_e-t_a}{t_b-t_a}\right)\right]

and for :math:`t_0-t_e< t_a`,

.. math:: t'=t_e+t_a+(t_b-t_a)\left[1-e^{-x}\right]

Triangular Distribution
=======================

The cumulative distribution function for the triangular distribution
with endpoints :math:`a` and :math:`b` and midpoint :math:`m` is

.. math::

   \begin{aligned}
     \frac{(t-a)^2}{(b-a)(m-a)} & & a\le t \le m \\
     1-\frac{(b-t)^2}{(b-a)(b-m)} & & m<t\le b.\end{aligned}

This makes the survival

.. math::

   \begin{aligned}
     1-\frac{(t-a)^2}{(b-a)(m-a)} & & a\le t \le m \\
     \frac{(b-t)^2}{(b-a)(b-m)} & & m<t\le b.\end{aligned}

Simple Sample
-------------

The cutoff is at :math:`t=m`, which is

.. math::

   \begin{aligned}
     U'&=&\frac{(m-a)^2}{(b-a)(m-a)} \\
     &=&\frac{m-a}{b-a}\end{aligned}

so first check whether :math:`U` is greater than that. Then, for
:math:`U` less than that,

.. math::

   \begin{aligned}
     t = a + \left[U(b-a)(m-a)\right]^{1/2} & & U\le U' \\
     t = b- \left[(1-U)(b-a)(b-m)\right]^{1/2} & & U'<U \\\end{aligned}

Shifted Sample
--------------

If this is sampled after some time, :math:`x`, then the thing we want to
invert is

.. math:: U=\frac{F(t)-F(x)}{G(x)}

so

.. math:: F(t)=UG(x)+F(x)

which, for :math:`a<t\le m` and :math:`a<x\le m`, is

.. math::

   \begin{aligned}
     \frac{(t-a)^2}{(b-a)(m-a)} & = & U\left[1-\frac{(x-a)^2}{(b-a)(m-a)}\right]
     +\frac{(x-a)^2}{(b-a)(m-a)} \\
     (t-a)^2 & = & U(b-a)(m-a)+ (1-U)(x-a)^2 \\
     t & = & a + \left[U(b-a)(m-a)+ (1-U)(x-a)^2\right]^{1/2}\end{aligned}

For :math:`m<t\le b` and :math:`m<x\le b`, this is

.. math::

   \begin{aligned}
     1-\frac{(b-t)^2}{(b-a)(b-m)} & = & U\frac{(b-x)^2}{(b-a)(b-m)}+
       1-\frac{(b-x)^2}{(b-a)(b-m)} \\
     -(b-t)^2 & = & U(b-x)^2-(b-x)^2 \\
     t & = & b-(b-x)\sqrt{1-U}\end{aligned}

In the case that :math:`m<t\le b` and :math:`a<x\le m`, the result is

.. math::

   \begin{aligned}
     1-\frac{(b-t)^2}{(b-a)(b-m)} & = & U\left[1-\frac{(x-a)^2}{(b-a)(m-a)}\right]
      + \frac{(x-a)^2}{(b-a)(m-a)}\\
     1-\frac{(b-t)^2}{(b-a)(b-m)} & = &U+\frac{(x-a)^2(1-U)}{(b-a)(m-a)} \\
     \frac{(b-t)^2}{(b-a)(b-m)} & = &(1-U)-\frac{(x-a)^2(1-U)}{(b-a)(m-a)} \\
     (b-t)^2& = &(1-U)(b-a)(b-m)-\frac{(x-a)^2(1-U)(b-m)}{(m-a)} \\
     t & = & b-\left[(1-U)(b-a)(b-m)-\frac{(x-a)^2(1-U)(b-m)}{(m-a)}\right]^{1/2}\end{aligned}

Sampling from Quantile
----------------------

The equation we have to solve is

.. math:: F(t)=1-G(t_j)(1-u)\prod_i\frac{G_i(t_i)}{G_i(t_{i+1})},

given to us as

.. math:: F(t)=1-G(t_j)(1-u)\gamma,

so, in terms of survivals, it’s

.. math:: G(t)=G(t_j)(1-u)\gamma

The value on the right is all known. If :math:`t_j<m`, the cutoff is

.. math:: G'=\frac{b-m}{(b-a)(b-m)}

Below that,
