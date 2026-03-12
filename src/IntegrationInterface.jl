module IntegrationInterface

export integral, Backend, Domain, Infinity

# function stubs for extensions to extend
function domainname end
function convert_domain end
function convert_integrand end
function ungroup end
function sum_domains end

include("types.jl")
include("serialization.jl")
include("domains.jl")
include("backends.jl")
include("integral.jl")
include("quadrature.jl")  # implementation of Backend.Quadrature, no extension required
include("infinity.jl")

end
