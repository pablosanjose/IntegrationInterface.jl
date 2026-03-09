module IntegrationInterface

export integral, Backends, Domains

# function stubs for extensions to extend
function domainname end
function convert_domain end
function check_domain_solver end

include("serialization.jl")
include("domains.jl")
include("solvers.jl")
include("quadrature.jl")  # implementation of Backends.Quadrature, no extension required
include("multi.jl")       # multi-solver for nested integrals, no extension required
include("integral.jl")

end
