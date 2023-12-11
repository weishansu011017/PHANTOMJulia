using LinearAlgebra
using Interpolations
using QuadGK

function cart2cylin(point::Vector)
    s = sqrt(point[1]^2 + point[2]^2)
    phi = atan(point[2],point[1])
    if size(point) == 2
        return [s,phi]
    else
        return [s,phi,point[3]]
    end
end

function cylin2cart(point::Vector)
    x = point[1]*cos(point[2])
    y = point[1]*sin(point[2])
    if size(point) == 2
        return [x,y]
    else
        return [x,y,point[3]]
    end
end

function Integral_1d(x::AbstractVector ,y::AbstractVector, inteval::Vector)
    spline = CubicSplineInterpolation(x, y, extrapolation_bc=Line())
    f_interp(x) = spline(x)
    integral, error = quadgk(f_interp, inteval[1], inteval[2])
    return integral, error
end

