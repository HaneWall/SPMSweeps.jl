abstract type SteadyStateChecker end

mutable struct Constant_Time_Check <: SteadyStateChecker
  converged::Bool
  function Constant_Time_Check()
    new(
      true
    )
  end
end

function check_converged!(S::Constant_Time_Check)
  S.converged = true
end

function check_converged!(S::Constant_Time_Check, measurement::Float64)
  S.converged = true
end


mutable struct Welford_Buffer <: SteadyStateChecker
  const N::Integer # window size of CircularBuffer
  mean::Float64
  variance_sum::Float64
  converged::Bool
  std_tolerance::Float64
  C::CircularBuffer{Float64}
  counter_call::Integer
  function Welford_Buffer(N, std_tolerance)
    c = CircularBuffer{Float64}(N)
    fill!(c, 0.0)
    new(
      N, 0.0, 0.0, false, std_tolerance, c, 1
    )
  end
end

function check_converged!(S::Welford_Buffer, measurement::Float64)
  new_mean = S.mean + (measurement - S.C[1]) / S.N
  S.variance_sum += (measurement + S.C[1] - S.mean - new_mean) * (measurement - S.C[1])
  S.mean = new_mean
  push!(S.C, measurement)
  S.counter_call += 1
  if sqrt(S.variance_sum / S.N) <= S.std_tolerance && S.counter_call > S.N
    S.converged = true
  else
    S.converged = false
  end
end
