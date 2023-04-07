module ComplexEconomicDynamics

import DynamicalSystemsBase: CoupledODEs, reinit!, step!
import LinearAlgebra: eigen, norm

export computemanifolds
include("dynamicalsystems/computemanifolds.jl")

import Plots: quiver!

export plotvectorfield, plotvectorfield!
include("plotting/plotvectorfield.jl")

end
