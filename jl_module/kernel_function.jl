# Some useful package


# Functional struct
struct KernelFunctionValid
    truncated_radius::Dict{Function, Float64}
end

function KernelFunctionValid()
    truncated_radius = Dict(
        M4_spline => 2.0,
        M5_spline => 2.5,
        M6_spline => 3.0,
        C2_Wendland => 2.0,
        C4_Wendland => 2.0,
        C6_Wendland => 2.0
    )
    return KernelFunctionValid(truncated_radius)
end

# Kernel Functions
function M4_spline(q::Float64)
    pir::Float64 = 1/pi
    Cnorm::Float64 = pir
    if 0 <= q < 1
        return Cnorm*(1 - (1.5 * q) + 0.75*(q^3))
    elseif 1 <= q < 2
        return Cnorm*(0.25 * (2 - q)^3)
    elseif q >= 2
        return 0.0
    else
        error("Kernal Error: Invaild numbers!")
    end
end

function M5_spline(q::Float64)
    pir::Float64 = 1/pi
    Cnorm::Float64 = 0.05*pir
    if 0 <= q < 0.5
        return Cnorm*((2.5 - q)^4 - 5*(1.5 - q)^4 + 10*(0.5 - q)^4)
    elseif 0.5 <= q < 1.5
        return Cnorm*((2.5 - q)^4 - 5*(1.5 - q)^4)
    elseif 1.5 <= q < 2.5
        return Cnorm*((2.5 - q)^4)
    elseif q >= 2.5
        return 0.0
    else
        error("Kernal Error: Invaild numbers!")
    end
end

function M6_spline(q::Float64)
    pir::Float64 = 1/pi
    Cnorm::Float64 = (120^-1)*pir
    if 0 <= q < 1
        return Cnorm*((3 - q)^5 - 6*(2 - q)^5 + 15*(1 - q)^5)
    elseif 1 <= q < 2
        return Cnorm*((3 - q)^5 - 6*(2 - q)^5)
    elseif 2 <= q < 3
        return Cnorm*((3 - q)^5)
    elseif q >= 3
        return 0.0
    else
        error("Kernal Error: Invaild numbers!")
    end
end

function C2_Wendland(q::Float64)
    pir::Float64 = 1/pi
    Cnorm::Float64 = 1.3125*pir #(21/16pi)
    if q < 2
        return Cnorm*(((1-0.5*q)^4)*(2*q + 1))
    elseif q >= 2
        return 0.0
    else
        error("Kernal Error: Invaild numbers!")
    end
end

function C4_Wendland(q::Float64)
    pir::Float64 = 1/pi
    Cnorm::Float64 = 1.9375*pir #(21/16pi)
    if q < 2
        return Cnorm*(((1-0.5*q)^6)*((35/12)*(q^2) + 3*q + 1))
    elseif q >= 2
        return 0.0
    else
        error("Kernal Error: Invaild numbers!")
    end
end

function C6_Wendland(q::Float64)
    pir::Float64 = 1/pi
    Cnorm::Float64 = 2.666015625*pir #(1365/512pi)
    if q < 2
        return Cnorm*(((1-0.5*q)^8)*(4*(q^3) + 6.25*(q^2) + 4*q + 1))
    elseif q >= 2
        return 0.0
    else
        error("Kernal Error: Invaild numbers!")
    end
end


function Smoothed_kernel_function(f::Function, h::Float32, ra::Vector, rb::Vector)
    """Note that the input function should be normalized properly."""
    r::Float64 = norm(ra - rb)
    q::Float64 = r/h
    h3r::Float64 = h^(-3)
    influence::Float64 = h3r * f(q)
    return influence
end

function Smoothed_kernel_function(f::Function, h::Float32, r::Float64)
    """Note that the input function should be normalized properly."""
    q::Float64 = r/h
    h3r::Float64 = h^(-3)
    influence::Float64 = h3r * f(q)
    return influence
end

function Smoothed_greident_kernel_function(f::Function, h::Float32, ra::Vector, rb::Vector)
    """Note that the input function should be normalized properly."""
    rab::Vector = ra - rb
    hatrab::Vector = rab./norm(rab)
    q::Float64 = norm(rab)/h
    h4r::Float64 = h^(-4)
    F_ab::Float64 = h4r * (ForwardDiff(f,q))
    influence::Vector = hatrab.*F_ab
    return influence
end

function Smoothed_dh_kernel_function(f::Function, h::Float32, ra::Vector, rb::Vector)
    """Note that the input function should be normalized properly."""
    r::Float64 = norm(ra - rb)
    q::Float64 = r/h
    h4r::Float64 = h^(-4)
    influence::Float64 = -h4r * (3*f(q) + q*ForwardDiff(f,q))
    return influence
end
