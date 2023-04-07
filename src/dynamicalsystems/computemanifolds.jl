"""
Numerically computes stable and unstable manifolds around each steadystate. It requires the user to provide the Jacobian, DF!, of the system. Follows the same syntax as `DynamicalSystems.jl`, such that 

```julia
F!(dx, x, p, t), DJ!(J, x, p, t)
```

# Keword Arguments

- T = 100: Number of time steps to compute the manifold
- limit = 1000: If the norm of the state vector |x| > limit, stop the computation and return NaN
- dt = 0.01: Time step size
- h = 1e-3: Step size in the direction of the eigenvectors around the steady states

The manifolds are stored in a 5-dimensional array of size (m, 2, 2, T, n), where m is the number of steady states, 2 is the number of eigenvectors, 2 is the number of directions (forward and backward), T is the number of time steps, and n is the number of state variables. The first index is the steady state, the second index is the eigenvector, the third index is the direction, the fourth index is the time step, and the fifth index is the state variable.
"""
function computemanifolds(
    F!::Function, DF!::Function,
    steadystates::Vector{Vector{Float64}},
    p::Vector{<:Real};
    T = 100, limit = 1000, dt = 0.01, h = 1e-3,
    kwargs...)

    m = length(steadystates)
    n = length(first(steadystates))

    function Finv!(dx, x, p, t)
        F!(dx, x, p, t)
        dx .= -dx
    end

    fwds = CoupledODEs(F!, zeros(n), p)
    bckds = CoupledODEs(Finv!, zeros(n), p)

    manifolds = NaN * ones(m, 2, 2, T, n)
        
    for (j, x̄) ∈ enumerate(steadystates)
        J = zeros(n, n); DF!(J, x̄, p, 0.0)
        λ, V = eigen(J)

        for (i, vᵢ) ∈ enumerate(eachcol(V))
            isstable = real(λ[i]) < 0 # Stable if real part of eigenvalue is negative

            ds = isstable ? bckds : fwds
            
            for (k, op) ∈ enumerate([-, +])
                reinit!(ds, op(x̄, vᵢ * h))

                for t ∈ 1:T
                    step!(ds, dt)
                    xₙ = ds.integ.u

                    manifolds[j, i, k, t, :] = xₙ

                    if norm(xₙ) > limit break end
                end         
            end
        end
    end

    return manifolds
end