using ComplexEconomicDynamics
using Test

@testset "ComplexEconomicDynamics.jl" begin

    function F!(dx, x, p, t)
        A = [ -2.0 p[1]; 1.0 2.0]
        dx .= A * x
    end

    function DF!(J, x, p, t)
        A = [ -2.0 p[1]; 1.0 2.0]
        J .= A
    end

    steadystates = [[0.0, 0.0]]
    manifolds = ComplexEconomicDynamics.computemanifolds(F!, DF!, steadystates, [1.0]; T = 10, limit = 10)

    @assert size(manifolds) == (1, 2, 2, 10, 2)
end
