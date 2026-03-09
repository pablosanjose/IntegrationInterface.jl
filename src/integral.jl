struct Integral{R,S<:AbstractBackend,F,D}
    integrand::F	# usually a function, but can be other kind of object
    result::R		# will be Missing if not in-place
    domain::D		# Tuple of AbstractDomains or a single AbstractDomain
    solver::S		# object representing the integration solver backend(s)
end

## API ##
function integral(f::F, domains...; result = missing, solver = default_solver(f, domains...)) where {F}
    domain´ = sanitize_domain(domains)
    solver´ = sanitize_solver(solver, domain´)
    check_domain_solver(domain´, solver´)
	return Integral(f, result, domain´, solver´)
end
# fallback. Can be overridden for user-defined f types and domains
default_solver(_...) = Backend.QuadGK()

ismutating(i::Integral) = !ismissing(i.result)

solver(i::Integral) = i.solver

solvername(i::Integral) = solvername(solver(i))
solvername(s::AbstractBackend) = nameof(typeof(s))

domain(i::Integral) = i.domain

domainname(i::Integral) = domainname(domain(i))
domainname(d::Tuple) = string("(", join(domainname.(d), ", "), ")")

## sanitization ##

sanitize_domain(domainfunc::Function) = Domain.Functional(domainfunc)
sanitize_domain(domain::Tuple{AbstractDomain}) = only(domain)
sanitize_domain(domains::Tuple) = sanitize_domain.(domains)   # nested domains
sanitize_domain(domain::AbstractDomain) = domain
sanitize_domain(domain) = throw(ArgumentError("Invalid domain specification $(domainname(domain))"))

sanitize_solver(solver::AbstractBackend, ::AbstractDomain) = solver
sanitize_solver(solver::AbstractBackend, ::NTuple{N,AbstractDomain}) where {N} =
    Nested(solver, Val(N))
sanitize_solver(solvers::NTuple{N,AbstractBackend}, ::NTuple{N,AbstractDomain}) where {N} =
    Nested(solvers)
sanitize_solver(solver, _) =
    throw(ArgumentError("Invalid solver specification $(solvername(solver)) for the given domain $(domainname(domain))."))

# error by default. Solvers must declare they understand the domain.
check_domain_solver(domain::AbstractDomain, solver::AbstractBackend) =
    error("The integral solver $(solvername(solver)) does not support the domain $(domainname(domain)), or solver backend not loaded.")

# for Nested, check each domain with corresponding solver.
function check_domain_solver(domains::NTuple{N,AbstractDomain}, solver::Nested{N}) where {N}
    last(domains) isa Domain.Functional &&
        throw(ArgumentError("Outermost (last) domain in nested integral cannot be a Domain.Functional."))
    # Functional non-outer domains are allowed for any solver in Nested (we cannot know its type until runtime)
    foreach(zip(domains, solvers(solver))) do (d, s)
        d isa Domain.Functional || check_domain_solver(d, s)
    end
    return nothing
end

## call syntax (scalar and in-place) ##
(i::Integral{Missing})(args...; params...) = call!(i, args...; params...)
(i::Integral)(args...; params...) = copy(call!(i, args...; params...))

call!(i::Integral, args...; params...) = i.solver(i.integrand, i.domain, i.result, args; params...)

(s::AbstractBackend)(f, domain, result, args; params...) =
    error("The integral solver solver for $(nameof(typeof(s))) is not loaded.")

## Show ##

Base.summary(i::Integral) = "Integral: callable object representing a numerical integral over a domain"

Base.show(io::IO, i::Integral) = print(io, summary(i), "\n",
"  Mutating   : $(ismutating(i))
  Domain     : $(domainname(i))
  Solver     : $(solvername(i))")
