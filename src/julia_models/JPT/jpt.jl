"""
```
JPT{T} <: AbstractRepModel{T}
```

The `JPT` type defines the structure of JPT, an implementation
of the DSGE model in Justiniano et al. (2010)
"Investment Shocks and Business Cycles".

### Fields

#### Parameters and Steady-States
* `parameters::Vector{AbstractParameter}`: Vector of all time-invariant model
  parameters.

* `steady_state::Vector{AbstractParameter}`: Model steady-state values, computed
  as a function of elements of `parameters`.

* `keys::OrderedDict{Symbol,Int}`: Maps human-readable names for all model
  parameters and steady-states to their indices in `parameters` and
  `steady_state`.

#### Inputs to Measurement and Equilibrium Condition Equations

The following fields are dictionaries that map human-readable names to row and
column indices in the matrix representations of of the measurement equation and
equilibrium conditions.

* `endogenous_states::OrderedDict{Symbol,Int}`: Maps each state to a column in
  the measurement and equilibrium condition matrices.

* `exogenous_shocks::OrderedDict{Symbol,Int}`: Maps each shock to a column in
  the measurement and equilibrium condition matrices.

* `expected_shocks::OrderedDict{Symbol,Int}`: Maps each expected shock to a
  column in the measurement and equilibrium condition matrices.

* `equilibrium_conditions::OrderedDict{Symbol,Int}`: Maps each equlibrium
  condition to a row in the model's equilibrium condition matrices.

* `endogenous_states_augmented::OrderedDict{Symbol,Int}`: Maps lagged states to
  their columns in the measurement and equilibrium condition equations. These
  are added after `gensys` solves the model.

* `observables::OrderedDict{Symbol,Int}`: Maps each observable to a row in the
  model's measurement equation matrices.

* `pseudo_observables::OrderedDict{Symbol,Int}`: Maps each pseudo-observable to
  a row in the model's pseudo-measurement equation matrices.

#### Model Specifications and Settings

* `spec::String`: The model specification identifier, \"m1002\", cached here for
  filepath computation.

* `subspec::String`: The model subspecification number, indicating that
  some parameters from the original model spec (\"ss10\") are initialized
  differently. Cached here for filepath computation.

* `settings::Dict{Symbol,Setting}`: Settings/flags that affect computation
  without changing the economic or mathematical setup of the model.

* `test_settings::Dict{Symbol,Setting}`: Settings/flags for testing mode

#### Other Fields

* `rng::MersenneTwister`: Random number generator. Can be is seeded to ensure
  reproducibility in algorithms that involve randomness (such as
  Metropolis-Hastings).

* `testing::Bool`: Indicates whether the model is in testing mode. If `true`,
  settings from `m.test_settings` are used in place of those in `m.settings`.

* `observable_mappings::OrderedDict{Symbol,Observable}`: A dictionary that
  stores data sources, series mnemonics, and transformations to/from model
  units. DSGE.jl will fetch data from the Federal Reserve Bank of St. Louis's
  FRED database; all other data must be downloaded by the user. See `load_data`
  and `Observable` for further details.

* `pseudo_observable_mappings::OrderedDict{Symbol,PseudoObservable}`: A
  dictionary that stores names and transformations to/from model units. See
  `PseudoObservable` for further details.
"""
mutable struct JPT{T} <: AbstractRepModel{T}
    parameters::ParameterVector{T}                         # vector of all time-invariant model parameters
    steady_state::ParameterVector{T}                       # model steady-state values
    keys::OrderedDict{Symbol,Int}                          # human-readable names for all the model
                                                           # parameters and steady-states

    endogenous_states::OrderedDict{Symbol,Int}             # these fields used to create matrices in the
    exogenous_shocks::OrderedDict{Symbol,Int}              # measurement and equilibrium condition equations.
    expected_shocks::OrderedDict{Symbol,Int}               #
    equilibrium_conditions::OrderedDict{Symbol,Int}        #
    endogenous_states_augmented::OrderedDict{Symbol,Int}   #
    observables::OrderedDict{Symbol,Int}                   #
    pseudo_observables::OrderedDict{Symbol,Int}            #

    spec::String                                           # Model specification number (eg "m990")
    subspec::String                                        # Model subspecification (eg "ss0")
    settings::Dict{Symbol,Setting}                         # Settings/flags for computation
    test_settings::Dict{Symbol,Setting}                    # Settings/flags for testing mode
    rng::MersenneTwister                                   # Random number generator
    testing::Bool                                          # Whether we are in testing mode or not

    observable_mappings::OrderedDict{Symbol, Observable}
    pseudo_observable_mappings::OrderedDict{Symbol, PseudoObservable}
end

description(m::JPT) = "Justiniano, Primiceri, and Tambalotti (2010), $(m.subspec)."

"""
`init_model_indices!(m::JPT)`

Arguments:
`m:: JPT`: a model object

Description:
Initializes indices for all of `m`'s states, shocks, and equilibrium conditions.
"""
function init_model_indices!(m::JPT)
    # Endogenous states
    endogenous_states = [:y_t, :k_t, :L_t, :Rk_t, :w_t, :π_t, :s_t, :λ_t,       # sticky prices
                         :c_t, :R_t, :u_t, :ϕ_t, :i_t, :kbar_t, :wgap_t,
                         :gdp_t, :z_t, :g_t, :μ_t, :λ_p_t, :λ_w_t,
                         :b_t, :mp_t, :Eπ_t, :Eλ_t, :Eϕ_t,
                         :Rk_t, :Ei_t, :Ew_t, :gdp_t1, :c_t1, :i_t1, :w_t1,     # CHECK BACK IF ANY CAN BE AUGMENTED VARIABLES
                         :y_f_t, :k_f_t, :L_f_t, :Rk_f_t, :w_f_t, :s_f_t,       # flexible prices
                         :λ_f_t, :c_f_t, :R_f_t, :u_f_t, :ϕ_f_t,
                         :i_f_t, :kbar_f_t, :wgap_f_t, :Ec_f_t, :Eλ_f_t,
                         :Eϕ_f_t, :ERk_f_t, :Ei_f_t]

    # Exogenous shocks
    exogenous_shocks = [:mp_sh, :z_sh, :g_sh, :μ_sh, :λ_p_sh, :λ_w_sh, :b_sh]

    # Expectations shocks
    expected_shocks = [:Eπ_sh, :Ec_sh, :Eλ_sh, :Eϕ_sh, :ERk_sh, :Ei_sh, :Ew_sh, # sticky prices
                       :Ec_f_sh, :Eλ_f_sh, :Eϕ_f_sh, :ERk_f_sh, :Ei_f_sh]       # flexible prices

    # Equilibrium conditions
    equilibrium_conditions = [:eq_y, :eq_k, :eq_L, :eq_Rk, :eq_w, :eq_π, :eq_s, :eq_λ,     # sticky prices
                              :eq_c, :eq_R, :eq_u, :eq_ϕ, :eq_i, :eq_kbar, :eq_wgap,
                              :eq_gdp, :eq_z, :eq_g, :eq_μ, :eq_λ_p, :eq_λ_w,
                              :eq_b, :eq_mp, :eq_Eπ, :eq_Eλ, :eq_Eϕ,
                              :eq_Rk, :eq_Ei, :eq_Ew, :eq_gdp1, :eq_c1, :eq_i1, :eq_w1,    # CHECK BACK IF ANY CAN BE AUGMENTED VARIABLES
                              :eq_y_f, :eq_k_f, :eq_L_f, :eq_Rk_f, :eq_w_f, :eq_s_f,       # flexible prices
                              :eq_λ_f, :eq_c_f, :eq_R_f, :eq_u_f, :eq_ϕ_f,
                              :eq_i_f, :eq_kbar_f, :eq_wgap_f, :eq_Ec_f, :eq_Eλ_f,
                              :eq_Eϕ_f, :eq_ERk_f, :eq_Ei_f]

    # Additional states added after solving model
    # Lagged states and observables measurement error
    endogenous_states_augmented = []

    # Observables
    observables = keys(m.observable_mappings)

    # Pseudo-observables
    pseudo_observables = keys(m.pseudo_observable_mappings)

    for (i,k) in enumerate(endogenous_states);           m.endogenous_states[k]           = i end
    for (i,k) in enumerate(exogenous_shocks);            m.exogenous_shocks[k]            = i end
    for (i,k) in enumerate(expected_shocks);             m.expected_shocks[k]             = i end
    for (i,k) in enumerate(equilibrium_conditions);      m.equilibrium_conditions[k]      = i end
    for (i,k) in enumerate(endogenous_states);           m.endogenous_states[k]           = i end
    for (i,k) in enumerate(endogenous_states_augmented); m.endogenous_states_augmented[k] = i+length(endogenous_states) end
    for (i,k) in enumerate(observables);                 m.observables[k]                 = i end
    for (i,k) in enumerate(pseudo_observables);          m.pseudo_observables[k]          = i end
end

function JPT(subspec::String="ss1";
                   custom_settings::Dict{Symbol, Setting} = Dict{Symbol, Setting}(),
                   testing = false)

    # Model-specific specifications
    spec               = split(basename(@__FILE__),'.')[1]
    subspec            = subspec
    settings           = Dict{Symbol,Setting}()
    test_settings      = Dict{Symbol,Setting}()
    rng                = MersenneTwister(0)

    # Initialize empty model
    m = JPT{Float64}(
            # model parameters and steady state values
            Vector{AbstractParameter{Float64}}(), Vector{Float64}(), OrderedDict{Symbol,Int}(),

            # model indices
            OrderedDict{Symbol,Int}(), OrderedDict{Symbol,Int}(), OrderedDict{Symbol,Int}(), OrderedDict{Symbol,Int}(), OrderedDict{Symbol,Int}(), OrderedDict{Symbol,Int}(), OrderedDict{Symbol,Int}(),

            spec,
            subspec,
            settings,
            test_settings,
            rng,
            testing,
            OrderedDict{Symbol,Observable}(),
            OrderedDict{Symbol,PseudoObservable}())

    # Set settings
    model_settings!(m)
    default_test_settings!(m)
    for custom_setting in values(custom_settings)
        m <= custom_setting
    end

    # Set observable and pseudo-observable transformations
    init_observable_mappings!(m)
    init_pseudo_observable_mappings!(m)

    # Initialize parameters
    init_parameters!(m)

    init_model_indices!(m)
    init_subspec!(m)
    steadystate!(m)

    return m
end

"""
```
init_parameters!(m::JPT)
```

Initializes the model's parameters, as well as empty values for the steady-state
parameters (in preparation for `steadystate!(m)` being called to initialize
those).
"""
function init_parameters!(m::JPT)
    # Calibrated parameters
    m <= parameter(:δ, 0.025, fixed=true,
                   description="δ: The capital depreciation rate.",
                   tex_label = "\\delta" )
    m <= parameter(:g_ss_ratio, 0.22, fixed=true,
                   description="g_ss_ratio: steady state government spending to GDP ratio.",
                   tex_label = "g^{ss}" )

    # Non standard devation parameters
    m <= parameter(:α, 0.1596, (0., 1.), (1e-5, 0.999), ModelConstructors.SquareRoot(), Normal(0.30, 0.05), fixed = false,
                   description="α: Capital elasticity in the intermediate goods sector's production function (also known as the capital share).",
                   tex_label = "\\alpha")

    m <= parameter(:ι_p, 0.1865, (0., 1.), (1e-5, 0.999), ModelConstructors.SquareRoot(), BetaAlt(0.5, 0.15), fixed = false,
                   description = "ι_p: The weight attributed to last period's inflation in price indexation. " *
                   "(1-ι_p) is the weight attributed to steady-state inflation.",
                   tex_label = "\\iota_p")

    m <= parameter(:ι_w, 0.2992, (0., 1.), (1e-5, 0.999), ModelConstructors.SquareRoot(), BetaAlt(0.5, 0.15), fixed = false,
                   description="ι_w: The weight attributed to last period's wage in wage indexation. "
                    * "(1-ι_w) is the weight attributed to steady-state wages.",
                   tex_label = "\\iota_w")

    m <= parameter(:γ100, 0.3673, (-5.0, 5.0), (-5., 5.), ModelConstructors.Untransformed(), Normal(0.4, 0.1), fixed = false,
                   scaling = x -> x/100,
                   description="γ: The log of the steady-state growth rate of technology.",
                   tex_label = "100\\gamma")

    m <= parameter(:h, 0.3673, (-5.0, 5.0), (-5., 5.), ModelConstructors.Untransformed(), Normal(0.4, 0.1), fixed = false,
                   description="h: habit formation parameter.",
                   tex_label = "h")

    m <= parameter(:λ_p_ss, 1.5000, fixed = false,
                   description="λ_p_ss: The steady state net price markup.",
                   tex_label = "\\lambda_p")

    m <= parameter(:λ_w_ss, 1.5000, fixed = false,
                   description = "λ_w_ss: The steady state net wage markup, which affects the elasticity of substitution between differentiated labor services.",
                   tex_label = "\\lambda_w")

    m <= parameter(:L_ss, -45.9364, (-1000., 1000.), (-1e3, 1e3), ModelConstructors.Untransformed(), Normal(-45., 5.), fixed=false,
                   description="L_ss: The steady state for log hours.",
                   tex_label = "\\log(L)^{ss}")

    m <= parameter(:π_ss100, 0.5000, (1e-5, 10.), (1e-5, 10.), ModelConstructors.Exponential(), GammaAlt(0.75, 0.4), fixed = false,
                   scaling = x -> 1 + x / 100,
                   description="π_ss100: The steady-state rate of net inflation multiplied by 100.",
                   tex_label = "100 \\pi^{ss}")

    m <= parameter(:Fβ, 0.1402, (1e-5, 10.), (1e-5, 10.), ModelConstructors.Exponential(), GammaAlt(0.25, 0.1), fixed=false,
                   scaling = x -> 1/(1 + x/100),
                   description = "Fβ: Discount rate transformed.",
                   tex_label = "100(\\beta^{-1} - 1)")

    m <= parameter(:ν, 2.5975, (1e-5, 10.), (1e-5, 10.), ModelConstructors.Exponential(), Normal(2, 0.75), fixed = false,
                   description="ν_l: The inverse Frisch elasticity.",
                   tex_label = "\\nu")

    m <= parameter(:ξ_p, 0.8940, (1e-5, 0.999), (1e-5, 0.999), ModelConstructors.SquareRoot(), BetaAlt(0.5, 0.1), fixed=false,
                   description="ξ_p: The Calvo parameter. In every period, intermediate goods producers optimize prices with probability (1-ξ_p). "
                   * "With probability ξ_p, prices are adjusted according to a weighted average of the previous period's inflation "
                   * "(π_t1) and steady-state inflation (π_ss).",
                   tex_label = "\\xi_p")

    m <= parameter(:ξ_w, 0.9291, (1e-5, 0.999), (1e-5, 0.999), ModelConstructors.SquareRoot(), BetaAlt(0.5, 0.1), fixed=false,
                   description="ξ_w: (1-ξ_w) is the probability with which households can freely choose wages in each period. "
                   * "With probability ξ_w, wages increase at a geometrically weighted average of the steady state rate of "
                   * "wage increases and last period's productivity times last period's inflation.",
                   tex_label = "\\xi_w")

    m <= parameter(:χ, 1.000, (0., 10.), (1e-5, 0.), ModelConstructors.Exponential(), GammaAlt(1., 0.5), fixed = false,
                   description="χ: The elasticity of the capital utilization cost function.",
                   tex_label = "\\chi")

    m <= parameter(:S′′, 2.7314, (-15., 15.), (-15., 15.), ModelConstructors.Untransformed(), Normal(4., 1.5), fixed = false,
                   description="S′′: The investment adjust cost.",
                   tex_label = "S^{\\prime\\prime}")

    m <= parameter(:ψ1, 1.3679, (1e-5, 10.), (1e-5, 10.00), ModelConstructors.Exponential(), Normal(1.5, 0.25), fixed = false,
                   description="ψ₁: Weight on inflation gap in monetary policy rule.",
                   tex_label = "\\psi_1")

    m <= parameter(:ψ2, 0.0388, (-0.5, 0.5), (-0.5, 0.5), ModelConstructors.Untransformed(), Normal(0.12, 0.05), fixed = false,
                   description="ψ₂: Weight on output gap in monetary policy rule.",
                   tex_label = "\\psi_2")

    m <= parameter(:ψ3, 0.2464, (-0.5, 0.5), (-0.5, 0.5), ModelConstructors.Untransformed(), Normal(0.12, 0.05), fixed = false,
                   description="ψ₃: Weight on rate of change of output gap in the monetary policy rule.",
                   tex_label = "\\psi_3")

    m <= parameter(:η_λ_p, 0.7892, (1e-5, 0.999), (1e-5, 0.999), ModelConstructors.SquareRoot(), BetaAlt(0.50, 0.20), fixed = false,
                   description="η_λ_p: Moving average component in the price markup shock.",
                   tex_label = "\\eta_{\\lambda_p}")

    m <= parameter(:η_λ_w, 0.4226, (1e-5, 0.999), (1e-5, 0.999), ModelConstructors.SquareRoot(), BetaAlt(0.50, 0.20), fixed = false,
                   description="η_λ_w: Moving average component in the wage markup shock.",
                   tex_label = "\\eta_{\\lambda_w}")

    # exogenous processes - autocorrelation
    m <= parameter(:ρ_R, 0.2135, (0., 1.), (1e-5, 0.999), ModelConstructors.SquareRoot(), BetaAlt(0.5, 0.2), fixed = false,
                   description="ρ_R: persistence in the monetary policy rule.",
                   tex_label = "\\rho_{R}")

    m <= parameter(:ρ_z, 0.9446, (1e-5, 0.999), (1e-5, 0.999), ModelConstructors.SquareRoot(), BetaAlt(0.5, 0.2), fixed = false,
                   description="ρ_z: AR(1) coefficient in the technology process.",
                   tex_label = "\\rho_{z}")

    m <= parameter(:ρ_g, 0.9863, (1e-5, 0.999), (1e-5, 0.999), ModelConstructors.SquareRoot(), BetaAlt(0.5, 0.2), fixed = false,
                   description="ρ_g: AR(1) coefficient in the government spending process.",
                   tex_label = "\\rho_g")

    m <= parameter(:ρ_μ, 0.8735, (1e-5, 0.999), (1e-5, 0.999), ModelConstructors.SquareRoot(), BetaAlt(0.5, 0.2), fixed = false,
                   description="ρ_μ: AR(1) coefficient in capital adjustment cost process.",
                   tex_label = "\\rho_{\\mu}")

    m <= parameter(:ρ_λ_p, 0.8827, (1e-5, 0.999), (1e-5, 0.999), ModelConstructors.SquareRoot(), BetaAlt(0.5, 0.2), fixed = false,
                   description="ρ_λ_p: AR(1) coefficient in the price mark-up shock process.",
                   tex_label = "\\rho_{\\lambda_p}")

    m <= parameter(:ρ_λ_w, 0.3884, (1e-5, 0.999), (1e-5, 0.999), ModelConstructors.SquareRoot(), BetaAlt(0.5, 0.2), fixed = false,
                   description="ρ_λ_w: AR(1) coefficient in the wage mark-up shock process.",
                   tex_label = "\\rho_{\\lambda_w}")

    m <= parameter(:ρ_b, 0.9410, (1e-5, 0.999), (1e-5, 0.999), ModelConstructors.SquareRoot(), BetaAlt(0.5, 0.2), fixed = false,
                   description="ρ_b: AR(1) coefficient in the intertemporal preference shifter process.",
                   tex_label = "\\rho_b")

    m <= parameter(:ρ_mp, 0.2135, (0., 1.), (1e-5, 0.999), ModelConstructors.SquareRoot(), BetaAlt(0.5, 0.2), fixed = false,
                   description="ρ_mp: AR(1) coefficient in the monetary policy shock process.",
                   tex_label = "\\rho_{mp}")

    # exogenous processes - standard deviation
    m <= parameter(:σ_mp, 0.2380, (1e-8, 5.), (1e-8, 5.), ModelConstructors.Exponential(), RootInverseGamma(2, 0.10), fixed = false,
                   description="σ_mp: The standard deviation of the monetary policy shock.",
                   tex_label = "\\sigma_{mp}")

    m <= parameter(:σ_z, 0.6742, (1e-8, 5.), (1e-8, 5.), ModelConstructors.Exponential(), RootInverseGamma(2, 0.10), fixed = false,
                   description="σ_z: The standard deviation of the technology process.",
                   tex_label = "\\sigma_{z}")

    m <= parameter(:σ_g, 2.5230, (1e-8, 5.), (1e-8, 5.), ModelConstructors.Exponential(), RootInverseGamma(2, 0.10), fixed = false,
                   description="σ_g: The standard deviation of the government spending process.",
                   tex_label = "\\sigma_{g}")

    m <= parameter(:σ_μ, 0.4559, (1e-8, 5.), (1e-8, 5.), ModelConstructors.Exponential(), RootInverseGamma(2, 0.10), fixed = false,
                   description="σ_μ: The standard deviation of the exogenous marginal efficiency of investment shock process.",
                   tex_label = "\\sigma_{\\mu}")

    m <= parameter(:σ_λ_p, 0.1314, (1e-8, 5.), (1e-8, 5.), ModelConstructors.Exponential(), RootInverseGamma(2, 0.10), fixed = false,
                   description="σ_λ_p: The mean of the process that generates the price elasticity of the composite good. " *
                   "Specifically, the elasticity is (1+λ_{f,t})/(λ_{p_t}).",
                   tex_label = "\\sigma_{\\lambda_p}")

    m <= parameter(:σ_λ_w, 0.3864, (1e-8, 5.), (1e-8, 5.), ModelConstructors.Exponential(), RootInverseGamma(2, 0.10), fixed = false,
                   tex_label = "\\sigma_{\\lambda_w}")

    m <= parameter(:σ_b, 0.0292, (1e-8, 5.), (1e-8, 5.), ModelConstructors.Exponential(), RootInverseGamma(2, 0.10), fixed = false,
                   description="σ_b: The standard deviation of the intertemporal preference shifter process.",
                   tex_label = "\\sigma_{b}")

    # steady states
    m <= SteadyStateParameter(:γ, NaN, tex_label = "\\gamma")
    m <= SteadyStateParameter(:β, NaN, tex_label = "\\beta")
    m <= SteadyStateParameter(:r_ss, NaN, tex_label = "r^{ss}")
    m <= SteadyStateParameter(:r_ss100, NaN, tex_label = "100 r^{ss}")
    m <= SteadyStateParameter(:π_ss, NaN, tex_label = "\\pi^{ss}")
    m <= SteadyStateParameter(:g_ss, NaN, tex_label = "G^{ss}")
    m <= SteadyStateParameter(:expL_ss, NaN, tex_label = "L^{ss}")
    m <= SteadyStateParameter(:Rk_ss, NaN, description = "Steady-state short-term rate of return on capital.", tex_label = "R^{k, ss}")
    m <= SteadyStateParameter(:s_ss, NaN, description = "Steady-state marginal cost", tex_label = "s^{ss}")
    m <= SteadyStateParameter(:kL_ss, NaN, tex_label = "k^{ss}/L^{ss}")
    m <= SteadyStateParameter(:FL_ss, NaN, description = "Steady-state ratio of fixed costs to labor", tex_label = "F^{ss}/L^{ss}")
    m <= SteadyStateParameter(:yL_ss, NaN, description = "Steady-state ratio of output to labor", tex_label = "y^{ss}/L^{ss}")
    m <= SteadyStateParameter(:k_ss, NaN, tex_label = "k^{ss}")
    m <= SteadyStateParameter(:i_ss, NaN, tex_label = "i^{ss}")
    m <= SteadyStateParameter(:F_ss, NaN, tex_label = "F^{ss}")
    m <= SteadyStateParameter(:y_ss, NaN, tex_label = "y^{ss}")
    m <= SteadyStateParameter(:c_ss, NaN, tex_label = "c^{ss}")
end

"""
```
steadystate!(m::JPT)
```

Calculates the model's steady-state values. `steadystate!(m)` must be called whenever
the parameters of `m` are updated.
"""
function steadystate!(m::JPT)

    m[:γ]        = m[:γ100] ./ 100.
    m[:β]        = 100. / (m[:Fβ] + 100.)
    m[:r_ss]     = exp(m[:γ]) / m[:β] - 1.
    m[:r_ss100]  = 100. * m[:r_ss]
    m[:π_ss]     = m[:π_ss100] ./ 100.
    m[:g_ss]     = 1. / (1. - m[:g_ss_ratio])
    m[:expL_ss]  = exp(m[:L_ss])
    m[:Rk_ss]    = exp(m[:γ]) / m[:β] - 1. + m[:δ]
    m[:s_ss]     = 1. / (1. + m[:λ_p_ss])
    m[:w_ss]     = (m[:s_ss] * ((1. - m[:α]) ^ (1. - m[:α])) / ((m[:α] ^ (-m[:α])) * m[:Rk_ss] ^ m[:α])) ^ (1. / (1. - m[:α]))
    m[:kL_ss]    = (m[:w_ss] / m[:Rk_ss]) ^ (m[:α] / (1. - m[:α]))
    m[:FL_ss]    = m[:kL_ss] ^ m[:α] - m[:Rk_ss] * m[:kL_ss] - m[:w_ss]
    m[:yL_ss]    = m[:kL_ss] ^ m[:α] - m[:FL_ss]
    m[:k_ss]     = m[:kL_ss] * m[:expL_ss]
    m[:i_ss]     = (1. - (1. - m[:δ]) * exp(-m[:γ])) * m[:k_ss] * exp(m[:γ])
    m[:F]        = m[:FL_ss] * m[:expL_ss]
    m[:y_ss]     = m[:yL_ss] * m[:expL_ss]
    m[:c_ss]     = m[:y_ss] / m[:g_ss] -  m[:i_ss]

    return m
end

function model_settings!(m::JPT)

    default_settings!(m)

    # Anticipated shocks
    m <= Setting(:n_anticipated_shocks, 0,
                 "Number of anticipated policy shocks")
    m <= Setting(:n_anticipated_shocks_padding, 0,
                 "Padding for anticipated policy shocks")

    # Data
    m <= Setting(:data_id, 3, "Dataset identifier")
    if get_setting(m, :cond_id) in collect(1:5)
        m <= Setting(:cond_full_names, [:obs_gdp, :obs_corepce, :obs_spread, :obs_nominalrate, :obs_longrate],
                     "Observables used in conditional forecasts")
    elseif get_setting(m, :cond_id) == 6
        m <= Setting(:cond_full_names, [:obs_gdp, :obs_corepce, :obs_spread, :obs_nominalrate, :obs_longrate,
                                        :obs_gdpdeflator], "Observables used in conditional forecasts")
    end
    m <= Setting(:cond_semi_names, [:obs_spread, :obs_nominalrate, :obs_longrate],
                 "Observables used in semiconditional forecasts")

    # Forecast
    m <= Setting(:use_population_forecast, true,
                 "Whether to use population forecasts as data")
    m <= Setting(:shockdec_startdate, Nullable(quartertodate("2007-Q1")),
                 "Date of start of shock decomposition output period. If null, then shockdec starts at date_mainsample_start")
end
