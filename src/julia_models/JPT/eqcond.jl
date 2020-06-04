"""
```
eqcond(m::JPT)
```
 Expresses the equilibrium conditions in canonical form using Γ0, Γ1, C, Ψ, and Π matrices.
Using the mappings of states/equations to integers defined in jpt.jl, coefficients are
specified in their proper positions.

Γ0 (n_states x n_states) holds coefficients of current time states.
Γ1 (n_states x n_states) holds coefficients of lagged states.
C  (n_states x 1) is a vector of constants
Ψ  (n_states x n_shocks_exogenous) holds coefficients of iid shocks.
Π  (n_states x n_states_expectational) holds coefficients of expectational states.
"""
function eqcond(m::JPT)
    endo = m.endogenous_states
    exo  = m.exogenous_shocks
    ex   = m.expected_shocks
    eq   = m.equilibrium_conditions

    Γ0 = zeros(n_states(m), n_states(m))
    Γ1 = zeros(n_states(m), n_states(m))
    C  = zeros(n_states(m))
    Ψ  = zeros(n_states(m), n_shocks_exogenous(m))
    Π  = zeros(n_states(m), n_shocks_expectational(m))

    ### ENDOGENOUS STATES ###

    ### 1. Production function

    # Sticky prices
    Γ0[eq[:eq_y], endo[:y_t]] = 1.
    Γ0[eq[:eq_y], endo[:k_t]] = -((m[:y_ss] + m[:F]) / m[:y_s]) * m[:α]
    Γ0[eq[:eq_y], endo[:L_t]] = -((m[:y_ss] + m[:F]) / m[:y_s]) * (1. - m[:α])

    # Flexible prices
    Γ0[eq[:eq_y_f], endo[:y_f_t]] = 1.
    Γ0[eq[:eq_y_f], endo[:k_f_t]] = -((m[:y_ss] + m[:F]) / m[:y_s]) * m[:α]
    Γ0[eq[:eq_y_f], endo[:L_f_t]] = -((m[:y_ss] + m[:F]) / m[:y_s]) * (1. - m[:α])

    ### 2. Cost minimization

    # Sticky prices
    Γ0[eq[:eq_L], endo[:Rk_t]] = 1.
    Γ0[eq[:eq_L], endo[:k_t]]  = 1.
    Γ0[eq[:eq_L], endo[:w_t]]  = -1.
    Γ0[eq[:eq_L], endo[:L_t]]  = -1.

    # Flexible prices
    Γ0[eq[:eq_L_f], endo[:Rk_f_t]] = 1.
    Γ0[eq[:eq_L_f], endo[:k_f_t]]  = 1.
    Γ0[eq[:eq_L_f], endo[:w_f_t]]  = -1.
    Γ0[eq[:eq_L_f], endo[:L_f_t]]  = -1.

    ### 3. Marginal cost

    # Sticky prices
    Γ0[eq[:eq_s], endo[:s_t]]  = 1.
    Γ0[eq[:eq_s], endo[:Rk_t]] = -m[:α]
    Γ0[eq[:eq_s], endo[:w_t]]  = -(1. - m[:α])

    # Flexible prices
    Γ0[eq[:eq_s_f], endo[:s_f_t]]  = 1.
    Γ0[eq[:eq_s_f], endo[:Rk_f_t]] = -m[:α]
    Γ0[eq[:eq_s_f], endo[:w_f_t]]  = -(1. - m[:α])

    ### 4. Phillips Curve

    # Sticky prices
    Γ0[eq[:eq_π], endo[:π_t]]   = 1.
    Γ0[eq[:eq_π], endo[:Eπ_t]]  = -m[:β] / (1 + m[:ι_p] * m[:β])
    Γ1[eq[:eq_π], endo[:π_t]]   = m[:ι_p] / (1 + m[:ι_p] * m[:β])
    Γ0[eq[:eq_π], endo[:s_t]]   = -(1. - m[:β] * m[:ξ_p]) * (1. - m[:ξ_p]) / ((1 + m[:ι_p] * m[:β]) * m[:ξ_p])
    Γ0[eq[:eq_π], endo[:λ_p_t]] = -1.

    # Flexible prices
    Γ0[eq[:eq_R_f], endo[:s_f_t]] = 1.

    ### 5. Consumption FOC

    # Sticky prices
    expγ = exp(m[:γ])
    Γ0[eq[:eq_c], endo[:λ_t]]  = (expγ - m[:h] * m[:β]) * (expγ - m[:h])
    Γ0[eq[:eq_c], endo[:b_t]]  = -(expγ - m[:h] * m[:β] * m[:ρ_b]) * (expγ - m[:h]) /
        ((1. - m[:ρ_b]) * (expγ - m[:h] * m[:β] * m[:ρ_b]) * (expγ - m[:h]) /
         (expγ * m[:h] + expγ ^ 2 + m[:β] * m[:h] ^ 2))
    Γ0[eq[:eq_c], endo[:z_t]]  = -(m[:β] * m[:h] * expγ * m[:rho_z] - m[:h] * expγ)
    Γ0[eq[:eq_c], endo[:c_t]]  = expγ ^ 2 + m[:β] * m[:h] ^ 2
    Γ0[eq[:eq_c], endo[:Ec_t]] = -m[:β] * m[:h] * expγ
    Γ1[eq[:eq_c], endo[:c_t]]  = m[:h] * expγ

    # Flexible prices
    expγ = exp(m[:γ])
    Γ0[eq[:eq_c_f], endo[:λ_f_t]]  = (expγ - m[:h] * m[:β]) * (expγ - m[:h])
    Γ0[eq[:eq_c_f], endo[:b_t]]  = -(expγ - m[:h] * m[:β] * m[:ρ_b]) * (expγ - m[:h]) /
        ((1. - m[:ρ_b]) * (expγ - m[:h] * m[:β] * m[:ρ_b]) * (expγ - m[:h]) /
         (expγ * m[:h] + expγ ^ 2 + m[:β] * m[:h] ^ 2))
    Γ0[eq[:eq_c_f], endo[:z_t]]  = -(m[:β] * m[:h] * expγ * m[:rho_z] - m[:h] * expγ)
    Γ0[eq[:eq_c_f], endo[:c_f_t]]  = expγ ^ 2 + m[:β] * m[:h] ^ 2
    Γ0[eq[:eq_c_f], endo[:Ec_f_t]] = -m[:β] * m[:h] * expγ
    Γ1[eq[:eq_c_f], endo[:c_f_t]]  = m[:h] * expγ

    ### 6. Euler equation

    # Sticky prices
    Γ0[eq[:eq_λ], endo[:λ_t]]  = 1.
    Γ0[eq[:eq_λ], endo[:R_t]]  = -1.
    Γ0[eq[:eq_λ], endo[:Eλ_t]] = -1.
    Γ0[eq[:eq_λ], endo[:Eπ_t]] = 1.
    Γ0[eq[:eq_λ], endo[:z_t]]  = m[:ρ_z]

    # Flexible prices
    Γ0[eq[:eq_λ_f], endo[:λ_f_t]]  = 1.
    Γ0[eq[:eq_λ_f], endo[:R_f_t]]  = -1.
    Γ0[eq[:eq_λ_f], endo[:Eλ_f_t]] = -1.
    Γ0[eq[:eq_λ_f], endo[:z_f_t]]  = m[:ρ_z]

    ### 7. Capital utilization FOC

    # Sticky prices
    Γ0[eq[:eq_Rk], endo[:Rk_t]] = 1.
    Γ0[eq[:eq_Rk], endo[:u_t]]  = -m[:χ]


    # Flexible prices

    return Γ0, Γ1, C, Ψ, Π
end
