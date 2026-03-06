struct Integral{R,S<:AbstractIntegralSolver,F,D}
    integrand::F	# usually a function, but can be other kind of object
    domain::D		# the type of domain object should be tailored to the chosen solver
	result::R		# will be Missing if not in-place
    solver::S		# object representing theintegration algorithm
end

## API ##
integral(f, result = missing; domain = (0,1), solver = default_integrator(f)) =
	Integral(f, sanitize_domain(domain, solver), result, solver)

default_integrator(_) = IS.QuadGK()

# by default, do nothing. We can override for specific domain/solver types
sanitize_domain(domain, _) = domain

ismutating(i::Integral) = !ismissing(i.result)

solver(i::Integral) = i.solver

solvername(i::Integral) = solvername(solver(i))

solvername(s::AbstractIntegralSolver) = nameof(typeof(s))

domain(i::Integral) = i.domain

## call syntax (scalar and in-place) ##
(i::Integral{Missing})(args...; params...) = call!(i, args...; params...)
(i::Integral)(args...; params...) = copy(call!(i, args...; params...))

call!(i::Integral, args...; params...) = i.solver(i.integrand, i.domain, i.result, args; params...)

(s::AbstractIntegralSolver)(f, domain, result, args; params...) =
    error("The integral solver backend for $(nameof(typeof(s))) is not loaded.")

## Show ##

Base.summary(i::Integral) = "Integral: callable object representing a numerical integral over a domain"

Base.show(io::IO, i::Integral) = print(io, summary(i), "\n",
"  Mutating   : $(ismutating(i))
  Domain     : $(domain(i))
  Solver     : $(solvername(i))")
