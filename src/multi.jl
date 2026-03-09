## Multisolver for nested integrals ##

struct Multi{N,S<:NTuple{N,Union{Function,AbstractBackend}}} <: AbstractBackend
  solvers::S
end

Multi(solver::AbstractBackend, N) = Multi(ntuple(Returns(solver), N))

solvers(s::Multi) = s.solvers

solvername(s::Multi) = string("Multi(", join(solvername.(s.solvers), ", "), ")")

# required for recursion
Base.front(s::Multi) = Multi(Base.front(s.solvers))
Base.last(s::Multi) = last(s.solvers)
Base.only(s::Multi) = only(s.solvers)
Base.length(s::Multi) = length(s.solvers)

sanitize_domain(domain, s::Multi) = length(domain) == length(s) ? domain :
    throw(ArgumentError("Wrong number of domains for Multi solver, expected $(length(s)), got $(length(domain))"))

## Call ##

# add an extra argument with outer point variables
(s::Multi)(f, domains, result, args; params...) = s(f, domains, result, args, (); params...)
# recursive solver call, starting from last. Only innermost is potentially mutating
(s::Multi)(f, domains, result, args, point; params...) =
	last(s)(maybe_evaluate_domain(last(domains), point), missing, args; params...) do z
		Base.front(s)(f, Base.front(domains), result, args, (point..., z...); params...)
	end
# Innermost integral
(s::Multi{1})(f, domains, ::Missing, args, point; params...) =
	only(s)((x, as...; ps...) -> f((x..., point...), as...; ps...), maybe_evaluate_domain(only(domains), point), missing, args; params...)
(s::Multi{1})(f!, domains, result, args, point; params...) =
	only(s)((out, x, as...; ps...) -> f!(out, (x..., point...), as...; ps...), maybe_evaluate_domain(only(domains), point), result, args; params...)

maybe_evaluate_domain(domain, _) = domain
maybe_evaluate_domain(domain::Domains.Functional, point) = domain.f(point...)
