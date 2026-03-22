## API ##
Integral(f::F, domain::AbstractDomain; result = nothing, backend::AbstractBackend = default_backend(domain), zerofastpath = false) where {F} =
	Integral(f, result, domain, backend, zerofastpath)
Integral(domain::AbstractDomain; kw...) = f -> Integral(f, domain; kw...)  # currying version

# autoevaluated integral
integral(f, domain::AbstractDomain, args...; result = nothing, backend::AbstractBackend = default_backend(domain), zerofastpath = false, kw...) =
    Integral(f, domain; result, backend, zerofastpath)(args...; kw...)
integral(domain::AbstractDomain, args...; kw...) = f -> integral(f, domain, args...; kw...)  # currying version

# Can be overridden for user-defined f types and domains
default_backend(::Domain.Box{1}) = Backend.QuadGK()
default_backend(::Domain.Sum{<:NTuple{<:Any,Domain.Box{1}}}) = Backend.QuadGK()
default_backend(::Domain.Functional{<:Domain.Box{1}}) = Backend.QuadGK()
default_backend(::Domain.Box) = Backend.HCubature()
default_backend(::Domain.Simplex) = Backend.HAdaptiveIntegration()
default_backend(::Domain.Functional{<:Domain.Simplex}) = Backend.HAdaptiveIntegration()
default_backend(::Domain.Sum{<:NTuple{<:Any,Domain.Box}}) = Backend.HCubature()
default_backend(::Domain.Functional{<:Domain.Box}) = Backend.HCubature()
default_backend(d) = throw(ArgumentError("No default backend exists for domain $(domainname(d)), please specify one explicitly."))

ismutating(i::Integral) = !isnothing(i.result)

integrand(i::Integral) = i.integrand

backend(i::Integral) = i.backend

result(i::Integral) = i.result

zerofastpath(i::Integral) = i.zerofastpath

backendname(i::Integral) = backendname(backend(i))
backendname(s::AbstractBackend) = nameof(typeof(s))
backendname(s::Type{<:AbstractBackend}) = nameof(s)

domain(i::Integral) = i.domain

domainname(i::Integral) = domainname(domain(i))

## Show ##

Base.summary(i::Integral) = "Integral"

function Base.show(io::IO, J::Integral)
    i = get(io, :indent, "")
    ioindent = IOContext(io, :indent => i * "  ")
    print(io, summary(J), "\n",
"$i  Mutating   : $(ismutating(J))
$i  Domain     : $(domainname(J))
$i  Backend    : $(backendname(J))
$i  Integrand  : ")
  print(ioindent, integrand(J))
end
