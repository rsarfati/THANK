"""
```
measurement(m::JPT{T}, TTT::Matrix{T}, RRR::Matrix{T},
            CCC::Vector{T}) where {T<:AbstractFloat}
```

Assign measurement equation

```
y_t = ZZ*s_t + DD + u_t
```

where

```
Var(ϵ_t) = QQ
Var(u_t) = EE
Cov(ϵ_t, u_t) = 0
```
"""
function measurement(m::JPT{T},
                     TTT::Matrix{T},
                     RRR::Matrix{T},
                     CCC::Vector{T}) where {T<:AbstractFloat}
    endo      = m.endogenous_states
    endo_addl = m.endogenous_states_augmented
    exo       = m.exogenous_shocks
    obs       = m.observables

    _n_observables = n_observables(m)
    _n_states = n_states_augmented(m)
    _n_shocks_exogenous = n_shocks_exogenous(m)

    ZZ = zeros(_n_observables, _n_states)
    DD = zeros(_n_observables)
    EE = zeros(_n_observables, _n_observables)
    QQ = zeros(_n_shocks_exogenous, _n_shocks_exogenous)

    ## Output growth - Quarterly!
    ZZ[obs[:obs_gdp], endo[:gdp_t]]       = 1.0
    ZZ[obs[:obs_gdp], endo_addl[:gdp_t1]] = -1.0
    ZZ[obs[:obs_gdp], endo[:z_t]]         = 1.0
    DD[obs[:obs_gdp]]                     = m[:γ100]

    ## Consumption Growth
    ZZ[obs[:obs_consumption], endo[:c_t]]       = 1.0
    ZZ[obs[:obs_consumption], endo_addl[:c_t1]] = -1.0
    ZZ[obs[:obs_consumption], endo[:z_t]]       = 1.0
    DD[obs[:obs_consumption]]                   = m[:γ100]

    ## Investment Growth
    ZZ[obs[:obs_investment], endo[:i_t]]       = 1.0
    ZZ[obs[:obs_investment], endo_addl[:i_t1]] = -1.0
    ZZ[obs[:obs_investment], endo[:z_t]]       = 1.0
    DD[obs[:obs_investment]]                   = m[:γ100]

    ## Hours growth
    ZZ[obs[:obs_hours], endo[:L_t]] = 1.0
    DD[obs[:obs_hours]]             = m[:L_ss]

    ## Real wage growth
    ZZ[obs[:obs_wages], endo[:w_t]]       = 1.0
    ZZ[obs[:obs_wages], endo_addl[:w_t1]] = -1.0
    ZZ[obs[:obs_wages], endo[:z_t]]       = 1.0
    DD[obs[:obs_wages]]                   = m[:γ100]

    ## Inflation (GDP Deflator)
    ZZ[obs[:obs_gdpdeflator], endo[:π_t]] = 1.0
    DD[obs[:obs_gdpdeflator]]             = m[:π_ss100]

    ## Nominal interest rate
    ZZ[obs[:obs_nominalrate], endo[:mp_t]] = 1.0
    DD[obs[:obs_nominalrate]]             = m[:π_ss100] + m[:r_ss100]

    # Variance of innovations
    QQ[exo[:mp_sh], exo[:mp_sh]]         = m[:σ_mp]^2
    QQ[exo[:z_sh], exo[:z_sh]]           = m[:σ_z]^2
    QQ[exo[:g_sh], exo[:g_sh]]           = m[:σ_g]^2
    QQ[exo[:μ_sh], exo[:μ_sh]]           = m[:σ_μ]^2
    QQ[exo[:λ_p_sh], exo[:λ_p_sh]]       = m[:σ_λ_p]^2
    QQ[exo[:λ_w_sh], exo[:λ_w_sh]]       = m[:σ_λ_w]^2
    QQ[exo[:b_sh], exo[:b_sh]]           = m[:σ_b]^2

    # Adjustment to DD because measurement equation assumes CCC is the zero vector
    if any(CCC .!= 0)
        DD += ZZ*((UniformScaling(1) - TTT)\CCC)
    end

    return Measurement(ZZ, DD, QQ, EE)
end
