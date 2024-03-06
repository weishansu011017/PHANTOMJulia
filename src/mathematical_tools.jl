function _cart2cylin(point::Vector)
    s = sqrt(point[1]^2 + point[2]^2)
    theta = atan(point[2],point[1])
    if length(point) == 2
        return [s,theta]
    else
        return [s,theta,point[3]]
    end
end

function _cylin2cart(point::Vector)
    x = point[1]*cos(point[2])
    y = point[1]*sin(point[2])
    if length(point) == 2
        return [x,y]
    else
        return [x,y,point[3]]
    end
end

function _Integral_1d(x::AbstractVector ,y::AbstractVector, inteval::Vector)
    spline = CubicSplineInterpolation(x, y, extrapolation_bc=Line())
    f_interp(x) = spline(x)
    integral, error = quadgk(f_interp, inteval[1], inteval[2])
    return integral
end

function value2closestvalueindex(array::AbstractVector, target::Float64)
    target_index = argmin(abs.(target .- array))
    return target_index
end

function find_array_max_index(y::AbstractVector)
    arange = 1:length(y)
    spline(x) = CubicSplineInterpolation(arange, y, extrapolation_bc=Line())(x)
    dspline(x) = ForwardDiff.derivative(spline, x)
    ddspline(x) = ForwardDiff.derivative(dspline, x)

    closest_to_zero = Inf
    target_index = 0
    for i in arange
        dspline_value = dspline(i)
        if abs(dspline_value) < closest_to_zero
            ddspline_value = ddspline(i)
            signature = sign(ddspline_value)
            if signature <= 0  # 包括极大值点和拐点
                closest_to_zero = abs(dspline_value)
                target_index = i
            end
        end
    end

    if target_index == 0
        error("SearchMaxError: The maximum value is not found.")
    end

    return target_index
end

function generate_weights(n)
    center = ceil(n / 2)
    weights = [max(abs(i - center), abs(j - center)) + 1 for i in 1:n, j in 1:n]
    max_weight = maximum(weights)
    normalized_weights = 1 .- (weights .- 1) / (max_weight - 1)
    return normalized_weights
end
function weighted_mean(matrix::Array)
    if length(matrix) == 1
        return matrix[1]
    end
    n = size(matrix, 1)
    weights = generate_weights(n)
    weighted_sum = sum(matrix .* weights)
    total_weight = sum(weights)
    return weighted_sum / total_weight
end
