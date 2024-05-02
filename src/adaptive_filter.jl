abstract type AdaptiveFilter end

mutable struct LMS_Algorithm <: AdaptiveFilter
  const μ::Float64
  const harmonics::Array{Float64}
  const basis_functions::Array{Function}
  basis_state::Array{Float64}
  harmonic_state::Array{Float64}

  function LMS_Algorithm(μ::Float64, k::Array{Float64})
    b_f = create_bases_psi(k)
    b = zeros(Float64, 2 * length(k) + 1)
    w = zeros(Float64, 2 * length(k) + 1)
    new(
      μ, k, b_f, b, w
    )
  end
end

function step_demodulation!(p::LMS_Algorithm, signal::Float64, reference_phase::Float64)
  update_basis!(p.basis_state, p.basis_functions, reference_phase)
  update_weights!(p.harmonic_state, p.basis_state, signal, p.μ)
end
