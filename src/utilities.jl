function create_bases_psi(harmonic_order::Array{Float64})
  bases = Function[]
  push!(bases, psi -> 1/2)
  @inbounds for harmonic in harmonic_order
    push!(bases, psi -> sin(harmonic * psi))
    push!(bases, psi -> cos(harmonic * psi))
  end
  return bases
end

@inline function update_basis!(b::Array{Float64}, bases::Array{Function}, psi::Float64)
  @inbounds for idx in eachindex(bases)
    b[idx] = bases[idx](psi)
  end
end

@inline function error_lms(data::Float64, weights::Array{Float64}, b::Array{Float64})
  return data - dot(weights, b)
end

@inline function update_weights!(w::Array{Float64}, b::Array{Float64}, data::Float64, mu::Float64)
  fact = 1. / 5. #inv(dot(b, b))
  error = error_lms(data, w, b)
  @inbounds for idx in eachindex(w)
    w[idx] += mu * fact * error * b[idx]
  end
end
