## Multisolver for nested integrals ##

struct Nested{N,S<:NTuple{N,Union{Function,AbstractBackend}}} <: AbstractBackend
  solvers::S
end

Nested(solver::AbstractBackend, N) = Nested(ntuple(Returns(solver), N))

solvers(s::Nested) = s.solvers

solvername(s::Nested) = string("Nested(", join(solvername.(s.solvers), ", "), ")")

# required for recursion
Base.front(s::Nested) = Nested(Base.front(s.solvers))
Base.last(s::Nested) = last(s.solvers)
Base.only(s::Nested) = only(s.solvers)
Base.length(s::Nested) = length(s.solvers)

sanitize_domain(domain, s::Nested) = length(domain) == length(s) ? domain :
    throw(ArgumentError("Wrong number of domains for Nested solver, expected $(length(s)), got $(length(domain))"))

## Call ##

# add an extra argument with outer point variables
(s::Nested)(f, domains, result, args; params...) = s(f, domains, result, args, (); params...)
# recursive solver call, starting from last. Only innermost is potentially mutating
(s::Nested)(f, domains, result, args, point; params...) =
	last(s)(maybe_evaluate_domain(last(domains), point), missing, args; params...) do z
		Base.front(s)(f, Base.front(domains), result, args, (point..., z...); params...)
	end
# Innermost integral
(s::Nested{1})(f, domains, ::Missing, args, point; params...) =
	only(s)((x, as...; ps...) -> f((x..., point...), as...; ps...), maybe_evaluate_domain(only(domains), point), missing, args; params...)
(s::Nested{1})(f!, domains, result, args, point; params...) =
	only(s)((out, x, as...; ps...) -> f!(out, (x..., point...), as...; ps...), maybe_evaluate_domain(only(domains), point), result, args; params...)

maybe_evaluate_domain(domain, _) = domain
maybe_evaluate_domain(domain::Domain.Functional, point) = domain.f(point...)
