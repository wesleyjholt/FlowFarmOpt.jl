abstract type OptAlgorithm end


"""
    SnoptWECAlgorithm(wec, tol, withTI)

Container for parameters for the Snopt with WEC optimization algorithm

# Arguments
- `wec::Array{Float64,1}`: WEC values (one for each optimization run)
- `tol::Array{Float64,1}`: convergence tolerance for optimizer
- `withTI::Array{Bool,1}`: whether or not to include TI in the wind farm model
- `maxiter::Float64`: maximum number of iterations
- `checkgradients::Bool`: whether or not to check gradients
"""
struct SnoptWECAlgorithm{AF, AB, F, I, B} <: OptAlgorithm
    objscale::F
    wec::AF
    tol::AF
    withTI::AB
    maxiter::I
    checkgradients::B
    parallel_processing::B
end
SnoptWECAlgorithm() = SnoptWECAlgorithm(1e-7, [3.0, 2.6, 2.2, 1.8, 1.4, 1.0, 1.0], [1e-3, 1e-3, 1e-4, 1e-4, 1e-5, 1e-6, 1e-6], [false, false, false, false, false, false, true], Int(1e6), true, false)
SnoptWECAlgorithm(a::Float64) = SnoptWECAlgorithm(a, [3.0, 2.6, 2.2, 1.8, 1.4, 1.0, 1.0], [1e-3, 1e-3, 1e-4, 1e-4, 1e-5, 1e-6, 1e-6], [false, false, false, false, false, false, true], Int(1e6), true, false)


"""
    SnoptAlgorithm(withTI)

Container for parameters for the Snopt with WEC optimization algorithm

# Arguments
- `withTI::Bool`: whether or not to include TI in the wind farm model
- `maxiter::Float64`: maximum number of iterations
- `checkgradients::Bool`: whether or not to check gradients
"""
struct SnoptAlgorithm{F, B} <: OptAlgorithm
    objscale::F
    withTI::B
    maxiter::F
    checkgradients::B
end
SnoptAlgorithm() = SnoptAlgorithm(1.0, true, 1e6, true)
