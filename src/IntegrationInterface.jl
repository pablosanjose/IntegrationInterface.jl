module IntegrationInterface

export integral, Backend, Domain

# function stubs for extensions to extend
function domainname end
function convert_domain end
function check_domain_solver end

include("serialization.jl")
include("domains.jl")
include("solvers.jl")
include("quadrature.jl")  # implementation of Backend.Quadrature, no extension required
include("integral.jl")

end
