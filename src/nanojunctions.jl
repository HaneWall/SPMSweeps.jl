abstract type Nanojunction end

mutable struct vdW_oscillator <: Nanojunction
  const k::Float64 # spring constant [N/m]
  const a_0::Float64 # intermolecular distance [m]
  const f::Float64 # excitation force [N]
  const H::Float64 # Hamaker constant [J]
  const R::Float64 # tip radius [m]
  const d::Float64 # distance point of mass and surface (base of cantilever) [m]
  const Q::Float64 # quality factor [1]
  Ω_s::Float64

  omegas::Array{Float64}

  targets::Array{Float64}
  target_idx::Integer 
  target_err::Float64 
  
  # rows will be target points, columns: [Ω_s + A.ctrl.output, X_0, X_1_s, X_1_c..., X_k_c]
  matrix_harmonic_state_converged_control::Matrix{Float64}
  matrix_harmonic_state_converged_sweep::Matrix{Float64}
  error::Array{Float64}
  current_ctrl::Array{Float64}
  
  ctrl::Controller
  filter::AdaptiveFilter
  checker::SteadyStateChecker

  # inner constructor

  function vdW_oscillator(k, f, Ω_s, a_0, d, H, R, Q, omegas, tar, err, S::Controller, T::AdaptiveFilter, U::SteadyStateChecker)
    return new(
      k, a_0, f, H, R, d, Q, Ω_s, omegas, tar, 1, err, 
      zeros(Float64, 2*length(T.harmonics) + 2,  length(tar)), 
      zeros(Float64, 2*length(T.harmonics) + 2,  length(omegas)), 
      Float64[], Float64[],
      S, T, U 
    )
  end

end

function f_RHS(x, A::vdW_oscillator, t)
  dx = x[2]
  dy = - 1/A.Q * x[2] - x[1] + A.f/A.k * sin(x[3])  
  if x[1] < A.d-A.a_0
    dy += A.H * A.R / (6 * A.k * (A.d - x[1])^2)
  else
    dy += A.H * A.R / (6 * A.a_0^2 * A.k) 
  end
  dz = A.Ω_s
  return SVector(dx, dy, dz)
end

function f_RHS_Ctrl(x, A::vdW_oscillator, t)
  dx = x[2]
  dy = - 1/A.Q * x[2] - x[1] + A.f/A.k * sin(x[3])  
  if x[1] < A.d-A.a_0
    dy += A.H * A.R / (6 * A.k * (A.d - x[1])^2)
  else
    dy += A.H * A.R / (6 * A.a_0^2 * A.k) 
  end
  dz = A.Ω_s + A.ctrl.output
  return SVector(dx, dy, dz)
end

mutable struct DMT_oscillator <: Nanojunction
  const k::Float64 # spring constant [N/m]
  const a_0::Float64 # intermolecular distance [m]
  const f::Float64 # excitation force [N]
  const H::Float64 # Hamaker constant [J]
  const R::Float64 # tip radius [m]
  const d::Float64 # distance point of mass and surface (base of cantilever) [m]
  const Q::Float64 # quality factor [1]
  Ω_s::Float64
  const E::Float64 # effective elasticity module between tip and sample 



  omegas::Array{Float64}

  targets::Array{Float64}
  target_idx::Integer 
  target_err::Float64 
  
  # rows will be target points, columns: [Ω_s + A.ctrl.output, X_0, X_1_s, X_1_c..., X_k_c]
  matrix_harmonic_state_converged_control::Matrix{Float64}
  matrix_harmonic_state_converged_sweep::Matrix{Float64}
  error::Array{Float64}
  current_ctrl::Array{Float64} 

  ctrl::Controller
  filter::AdaptiveFilter
  checker::SteadyStateChecker

  # inner constructor

  function DMT_oscillator(k, f, Ω_s, a_0, d, H, R, Q, nu_t, nu_s, E_t, E_s, omegas, tar, err, S::Controller, T::AdaptiveFilter, U::SteadyStateChecker)
    return new(
      k, a_0, f, H, R, d, Q, Ω_s, 
      inv((1-nu_s^2)/E_s + (1-nu_t^2)/E_t),
      omegas, tar, 1, err, 
      zeros(Float64, 2*length(T.harmonics) + 2,  length(tar)), 
      zeros(Float64, 2*length(T.harmonics) + 2,  length(omegas)), 
      Float64[], Float64[],
      S, T, U 
    )
  end
end


function f_RHS(x, A::DMT_oscillator, t)
  dx = x[2]
  dy = - 1/A.Q * x[2] - x[1] + A.f/A.k * sin(x[3])  
  if x[1] < A.d-A.a_0
    dy += A.H * A.R / (6 * A.k * (A.d - x[1])^2)
  else
    dy += A.H * A.R / (6* A.a_0^2 * A.k) - 4/3 * A.E * sqrt(A.R) * 1/A.k * (x[1] - (A.d - A.a_0))^(3/2)
  end
  dz = A.Ω_s
  return SVector(dx, dy, dz)
end

function f_RHS_Ctrl(x, A::DMT_oscillator, t)
  dx = x[2]
  dy = - 1/A.Q * x[2] - x[1] + A.f/A.k * sin(x[3])  
  if x[1] < A.d-A.a_0
    dy += A.H * A.R / (6 * A.k * (A.d - x[1])^2)
  else
    dy += A.H * A.R / (6 * A.a_0^2 * A.k) - 4/3 * A.E * sqrt(A.R) * 1/A.k * (x[1] - (A.d - A.a_0))^(3/2)
  end
  dz = A.Ω_s + A.ctrl.output
  return SVector(dx, dy, dz)
end


mutable struct Lennard_Jones_oscillator <: Nanojunction
  const k::Float64 # spring constant [N/m]
  const γ::Float64 # friction parameter
  const f::Float64 # excitation force [N]
  const V_0::Float64 # potential constant [J]
  const d::Float64 # distance point of mass and surface (base of cantilever) [m]
  const Q::Float64 # quality factor [1]
  const δx::Float64 # softening parameter [m]
  const σ::Float64 # determines minimum [m]
  const ω_0::Float64 # eigenfrequncy [rad/s]
  Ω_s::Float64


  omegas::Array{Float64}

  targets::Array{Float64}
  target_idx::Integer 
  target_err::Float64 
  
  # rows will be target points, columns: [Ω_s + A.ctrl.output, X_0, X_1_s, X_1_c..., X_k_c]
  matrix_harmonic_state_converged_control::Matrix{Float64}
  matrix_harmonic_state_converged_sweep::Matrix{Float64}
  error::Array{Float64}
  current_ctrl::Array{Float64} 

  ctrl::Controller
  filter::AdaptiveFilter
  checker::SteadyStateChecker

  # inner constructor

  function Lennard_Jones_oscillator(k, f, V_0, γ, Ω_s, d, σ, Q, δx, ω_0, omegas, tar, err, S::Controller, T::AdaptiveFilter, U::SteadyStateChecker)
    return new(
      k, γ, f, V_0, d, Q, δx, σ, ω_0, Ω_s,
      omegas, tar, 1, err, 
      zeros(Float64, 2*length(T.harmonics) + 2,  length(tar)), 
      zeros(Float64, 2*length(T.harmonics) + 2,  length(omegas)), 
      Float64[], Float64[],
      S, T, U 
    )
  end
end

function f_RHS(x, A::Lennard_Jones_oscillator, t)
  h_x = (x[1] - A.d)^2 + (A.δx)^2
  dx = x[2]
  dy = -1/A.Q * x[2] - x[1] - 12*A.V_0/(A.k * sqrt(h_x)) * ((A.σ^2 / h_x)^6 - (A.σ^2 / h_x)^3) - x[2]*A.ω_0*A.γ/(A.d- x[1])^3 + A.f*sin(x[3])
  dz = A.Ω_s
  return SVector(dx, dy, dz)
end

function f_RHS_Ctrl(x, A::Lennard_Jones_oscillator, t)
  h_x = (x[1] - A.d)^2 + (A.δx)^2
  dx = x[2]
  dy = -1/A.Q * x[2] - x[1] - 12*A.V_0/(A.k * sqrt(h_x)) * ((A.σ^2 / h_x)^6 - (A.σ^2 / h_x)^3) - x[2]*A.ω_0*A.γ/(A.d- x[1])^3 + A.f*sin(x[3])
  dz = A.Ω_s + A.ctrl.output
  return SVector(dx, dy, dz)
end

mutable struct Duffing_oscillator <: Nanojunction
  const k_1::Float64
  const k_2::Float64
  const k_3::Float64
  const f::Float64
  Ω_s::Float64

  omegas::Array{Float64}

  targets::Array{Float64}
  target_idx::Integer
  target_err::Float64

  # rows will be target points, columns: [Ω_s + A.ctrl.output, X_0, X_1_s, X_1_c..., X_k_c]
  matrix_harmonic_state_converged_control::Matrix{Float64}
  matrix_harmonic_state_converged_sweep::Matrix{Float64}
  error::Array{Float64}
  current_ctrl::Array{Float64}
  
  ctrl::Controller
  filter::AdaptiveFilter
  checker::SteadyStateChecker

  function Duffing_oscillator(k_1, k_2, k_3, f, Ω_s, omegas, tar, err, S::Controller, T::AdaptiveFilter, U::SteadyStateChecker)
      return new(
        k_1, k_2, k_3, f, Ω_s, omegas, tar, 1, err, 
        zeros(Float64, 2*length(T.harmonics) + 2,  length(tar)), 
        zeros(Float64, 2*length(T.harmonics) + 2,  length(omegas)),
        Float64[], Float64[],
        S, T, U
      )
  end
end


function f_RHS(x::SVector, A::Duffing_oscillator, t)
  dx = x[2] # displacement
  dy = A.f * sin(x[3]) - A.k_1 * x[1] - A.k_2 * x[2] - A.k_3 * x[1]^3 # velocity
  dz = A.Ω_s # instantanoes phase
  return @SVector [dx, dy, dz]
end

function f_RHS_Ctrl(x::SVector, A::Duffing_oscillator, t)
  dx = x[2] # displacement
  dy = A.f * sin(x[3]) - A.k_1 * x[1] - A.k_2 * x[2] - A.k_3 * x[1]^3 # velocity
  dz = A.Ω_s + A.ctrl.output # instantanoes phase
  return @SVector [dx, dy, dz]
end

