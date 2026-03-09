## Quadrature ##
# this, unlike other solvers, does not require an extension

check_domain_solver(::Union{Domain.Segment, Domain.SegmentGroup}, ::Backend.Quadrature) = nothing

convert_domain(d::Domain.Segment, ::Backend.Quadrature) = (d.x1, d.x2)

(s::Backend.Quadrature)(f, domain::Domain.SegmentGroup, result, args; params...) =
    sum(ungroup(domain)) do segment
        s(f, segment, missing, args; params...)
    end

function (s::Backend.Quadrature)(f, domain, ::Missing, args; params...)
    xmin, xmax = convert_domain(domain, s)
    Δx = 0.5 * (xmax - xmin)
    result = sum(zip(s.nodes, s.weights); init = zero(float(Δx))) do (x, w)
        f(xmin + Δx * (x + 1), args...; params...) * w * Δx
    end
    return result
end

function (s::Backend.Quadrature)(f!, domain, result, args; params...)
    xmin, xmax = convert_domain(domain, s)
    Δx = 0.5 * (xmax - xmin)
    fill!(result, 0)
    out = similar(result)
    result = foreach(zip(s.nodes, s.weights)) do (x, w)
        result .+= f!(out, xmin + Δx * (x + 1), args...; params...) .* w .* Δx
    end
    return result
end
