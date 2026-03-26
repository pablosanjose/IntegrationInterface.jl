## API ##
Integral(f::F, domain; result = nothing, backend = Backend.default(domain)) where {F} =
	Integral(maybe_mutating(f, result), result, maybe_functional(domain), backend)
Integral(domain; kw...) = f -> Integral(f, domain; kw...)  # currying version

# autoevaluated integral (need two methods due to ambiguities)
integral(f, domain, args...; result = nothing, backend = Backend.default(domain), kw...) =
    Integral(f, domain; result, backend)(args...; kw...)

maybe_functional(domain::AbstractDomain) = domain
maybe_functional(domain) = Domain.Functional(domain)  # assumed callable

# if the integrand J is an Integral, mutation is handled externally to J
# (we should not pass `out` to J and its domains)
maybe_mutating(J::Integral, ::Nothing) = J
maybe_mutating(J::Integral, _) = Mutating(J)
maybe_mutating(f, _) = f

ismutating(i::Integral) = !isnothing(i.result)

integrand(i::Integral) = i.integrand

backend(i::Integral) = i.backend

result(i::Integral) = i.result

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

Base.show(io::IO, J::Mutating) = (print(io, "[Mutating] "); show(io, J.f))
