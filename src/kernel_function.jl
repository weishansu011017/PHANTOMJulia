# Functional constant
const TRUNCATED_RADIUS = Dict(
    :M4_spline => 2.0,
    :M5_spline => 2.5,
    :M6_spline => 3.0,
    :C2_Wendland => 2.0,
    :C4_Wendland => 2.0,
    :C6_Wendland => 2.0
)

const CNORM = Dict(
    :M4_spline => [2/3, 10/(7*pi), pi],
    :M5_spline => [1/24, 96/(1199*pi), 0.05*pi],
    :M6_spline => [120^(-1), 7/(478*pi), 120^(-1)*pi],
    :C2_Wendland => [5/8, 7/(4*pi), 21/(16*pi)],
    :C4_Wendland => [3/4, 9/(4*pi), 495/(256*pi)],
    :C6_Wendland => [64/55, 78/(28*pi), 1365/(512*pi)]
)
function KernelFunctionValid()
    return TRUNCATED_RADIUS
end

function KernelFunctionnorm()
    return CNORM
end

# Kernel Functions
function M4_spline(q::T) where T
    if 0 <= q < 1
        return (1 - (1.5 * q) + 0.75*(q^3))
    elseif 1 <= q < 2
        return (0.25 * (2 - q)^3)
    elseif q >= 2
        return 0.0
    else
        error("Kernal Error: Invaild numbers!")
    end
end

function M5_spline(q::T) where T
    if 0 <= q < 0.5
        return ((2.5 - q)^4 - 5*(1.5 - q)^4 + 10*(0.5 - q)^4)
    elseif 0.5 <= q < 1.5
        return ((2.5 - q)^4 - 5*(1.5 - q)^4)
    elseif 1.5 <= q < 2.5
        return ((2.5 - q)^4)
    elseif q >= 2.5
        return 0.0
    else
        error("Kernal Error: Invaild numbers!")
    end
end

function M6_spline(q::T) where T
    if 0 <= q < 1
        return ((3 - q)^5 - 6*(2 - q)^5 + 15*(1 - q)^5)
    elseif 1 <= q < 2
        return ((3 - q)^5 - 6*(2 - q)^5)
    elseif 2 <= q < 3
        return ((3 - q)^5)
    elseif q >= 3
        return 0.0
    else
        error("Kernal Error: Invaild numbers!")
    end
end


function C2_Wendland(q::T) where T
    if q < 2
        return (((1-0.5*q)^4)*(2*q + 1))
    elseif q >= 2
        return 0.0
    else
        error("Kernal Error: Invaild numbers!")
    end
end

function C4_Wendland(q::T) where T
    if q < 2
        return (((1-0.5*q)^6)*((35/12)*(q^2) + 3*q + 1))
    elseif q >= 2
        return 0.0
    else
        error("Kernal Error: Invaild numbers!")
    end
end

function C6_Wendland(q::T) where T
    if q < 2
        return (((1-0.5*q)^8)*(4*(q^3) + 6.25*(q^2) + 4*q + 1))
    elseif q >= 2
        return 0.0
    else
        error("Kernal Error: Invaild numbers!")
    end
end


function Smoothed_kernel_function(f::Function, h::Float32, ra::Vector, rb::Vector)
    r::Float64 = norm(ra - rb)
    q::Float64 = r/h
    dim::Int32 = length(ra)
    hr::Float64 = h^(-dim)
    influence::Float64 = hr * f(q) * KernelFunctionnorm()[nameof(f)][dim]
    return influence
end

function Smoothed_kernel_function(f::Function, h::Float32, r::Float64,dim::Int)
    q::Float64 = r/h
    hr::Float64 = h^(-dim)
    influence::Float64 = hr * f(q) * KernelFunctionnorm()[nameof(f)][dim]
    return influence
end

function Smoothed_greident_kernel_function(f::Function, h::Float32, ra::Vector, rb::Vector)
    rab::Vector = ra - rb
    hatrab::Vector = rab./norm(rab)
    q::Float64 = norm(rab)/h
    dim::Int32 = length(rab)
    hr::Float64 = h^(-(dim+1))
    F_ab::Float64 = hr * (ForwardDiff(f,q)) * KernelFunctionnorm()[nameof(f)][dim]
    influence::Vector = hatrab.*F_ab
    return influence
end

function Smoothed_greident_kernel_function(f::Function, h::Float32, rab::Vector)
    hatrab::Vector = rab./norm(rab)
    q::Float64 = norm(rab)/h
    dim::Int32 = length(rab)
    hr::Float64 = h^(-(dim+1))
    F_ab::Float64 = hr * (ForwardDiff.derivative(f,q)) * KernelFunctionnorm()[nameof(f)][dim]
    influence::Vector = hatrab.*F_ab
    return influence
end

function Smoothed_dh_kernel_function(f::Function, h::Float32, ra::Vector, rb::Vector)
    """Note that the input function should be normalized properly."""
    r::Float64 = norm(ra - rb)
    dim::Int32 = length(ra)
    q::Float64 = r/h
    hr::Float64 = h^(-(dim+1))
    influence::Float64 = -hr * (3*f(q) + q*ForwardDiff(f,q)) * KernelFunctionnorm()[nameof(f)][dim]
    return influence
end

