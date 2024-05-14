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
  var_tolerance::Float64
  C::CircularBuffer{Float64}
  function Welford_Buffer(N, var_tolerance)
   new(
    N, 0., 0., 0., false, var_tolerance, CircularBuffer{Float64}(N)
   ) 
  end
end

function check_converged!(S::Welford_Buffer, measurement::Float64)
  new_mean = wb.mean + (measurement - S.C[1]) / S.window_size
  S.variance_sum += (measurement + S.C[1] - S.mean - new_mean) * (measurement - S.C[1])
  S.mean = new_mean
  push!(S.C, measurement)
  if sqrt(S.variance_sum/S.window_size) <= S.var_tolerance && S.mean 
    S.converged = true
  else
    S.converged = false
  end
end

