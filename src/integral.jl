struct Integral{R,S<:AbstractBackend,F,D}
    integrand::F	# usually a function, but can be other kind of object
    result::R		# will be Missing if not in-place
    domain::D		# Tuple of AbstractDomains or a single AbstractDomain
    solver::S		# object representing the integration solver backend(s)
end

## API ##
function integral(f::F, domain; result = missing, solver = default_solver(f, domain)) where {F}
    domain´ = sanitize_domain(domain)
    solver´ = sanitize_solver(solver, domain´)
    check_domain_solver(domain´, solver´)
	return Integral(f, result, domain´, solver´)
end

# currying version
integral(domain; kw...) = f -> integral(f, domain; kw...)

# fallback. Can be overridden for user-defined f types and domains
default_solver(_...) = Backend.QuadGK()

ismutating(i::Integral) = !ismissing(i.result)

integrand(i::Integral) = i.integrand

solver(i::Integral) = i.solver

solvername(i::Integral) = solvername(solver(i))
solvername(s::AbstractBackend) = nameof(typeof(s))
solvername(s::Type{<:AbstractBackend}) = nameof(s)

domain(i::Integral) = i.domain

domainname(i::Integral) = domainname(domain(i))

## sanitization ##

sanitize_domain(domain::AbstractDomain) = domain
sanitize_domain(domain) = throw(ArgumentError("Invalid domain specification $(domainname(domain))"))

sanitize_solver(solver::AbstractBackend, ::AbstractDomain) = solver
sanitize_solver(solver, domain) =
    throw(ArgumentError("Invalid solver specification $(solvername(solver)) for the given domain $(domainname(domain)), or solver backend not loaded."))

# error by default. Solvers must declare they understand the domain.

# we translate to typeof(domain) to allow Functional to be checked automatically
check_domain_solver(domain::AbstractDomain, solver::AbstractBackend) = check_domain_solver(typeof(domain), solver)
check_domain_solver(domain::Domain.Functional, solver::AbstractBackend) = check_domain_solver(domain.type, solver)
check_domain_solver(domain::Type{<:AbstractDomain}, solver::AbstractBackend) =
    error("The integral solver $(solvername(solver)) does not support $(nameof(domain)) domains, or solver backend not loaded.")

## call syntax (scalar and in-place) ##
(i::Integral{Missing})(args...; params...) = call!(i, args...; params...)
(i::Integral)(args...; params...) = copy(call!(i, args...; params...))

call!(i::Integral, args...; params...) = i.solver(i.integrand, i.domain, i.result, args; params...)

(s::AbstractBackend)(f, domain, result, args; params...) =
    error("The integral solver solver for $(nameof(typeof(s))) is not loaded.")

## Show ##

Base.summary(i::Integral) = "Integral"

function Base.show(io::IO, J::Integral)
    i = get(io, :indent, "")
    ioindent = IOContext(io, :indent => i * "  ")
    print(io, summary(J), "\n",
"$i  Mutating   : $(ismutating(J))
$i  Domain     : $(domainname(J))
$i  Solver     : $(solvername(J))
$i  Integrand  : ")
  print(ioindent, integrand(J))
end
