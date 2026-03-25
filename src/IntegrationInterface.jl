module IntegrationInterface

export integral, Integral, Backend, Domain, Infinity, witherror

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
include("quadrature.jl")       # implementation of Backend.Quadrature, no extension required
include("integral.jl")
include("interface.jl")
include("changeofvariables.jl")
include("docstrings.jl")

end
