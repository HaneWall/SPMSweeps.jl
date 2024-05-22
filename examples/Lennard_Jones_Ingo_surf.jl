using SPMsweeps
using StaticArrays, DifferentialEquations, CairoMakie, LinearAlgebra

initial_condition = SVector(0.0, 0.0, 0.0)

# params from Ingo sims / Annas dissertation / my master thesis
σ = 2.8e-9
V_0 = 4.2e-18
# softening parameter
δx = 0.0e-9
# distance betweeen cantilever and surface/sample
ds = collect(range(60e-9, 120e-9, length=20))
# Quality-factor of cantilever
Q = 400.0
# spring constant of cantilever
k = 0.7
# resonance frequency
ω_0 = 55.e3 * 2π
# friction parameter
γ = 0.6e-31
#forcing amplitude
f = 120.e-9 / Q




# time params sweep
Δ_t_filt = 0.05
Δ_t_sweep = 10_000.0
μ_sweep = 0.05

# time params control
Δ_t_ctrl = 0.08 # timestep for each control step
μ_ctrl = 0.05
Δ_t_checker = 15.3 # timestep in which we check convergence 


μ = Δ_t_ctrl # stepsize LMS (tends to be equal to sample time of control)
harms = [1.0, 2.0, 3.0, 4.0] # respected higher harmonics (DC always automatically included)
K_P = 0.1
K_I = 0.001
K_D = 1.5
τ = 5.0
int_min = -0.05
int_max = 0.05
ctrl_min = -0.051
ctrl_max = 0.051

# for control 
# targets for control (not relevant for sweeps)
TARGETS = collect(range(-0.2, -2.9, length=60))
err = 0.02

# frequncy array that we would like to sweep through (not relevant for control)
OMEGAS = collect(range(0.985, 1.015, length=60))

pll_problems = []
for d in ds
  CTRL = PID_Controller_Tustin(Δ_t_ctrl, K_P, K_I, K_D, τ, int_min, int_max, ctrl_min, ctrl_max)
  FILT = LMS_Algorithm(μ_ctrl, harms)
  CHECK = Welford_Buffer(40, 1e-3)
  pll_problem = Lennard_Jones_oscillator(k, f, V_0, γ, OMEGAS[1], d, σ, Q, δx, ω_0, OMEGAS, TARGETS, err, CTRL, FILT, CHECK)
  push!(pll_problems, pll_problem)
end

# for control 
control_cb = PeriodicCallback(ctrl_cb!, Δ_t_ctrl)
convergence_cb = PeriodicCallback(conv_wf_cb!, Δ_t_checker)
all_cb_control = CallbackSet(control_cb, convergence_cb)



t_span = (0, 200_000)

control_sols = []

for pll_problem in pll_problems
  control_prob = ODEProblem(f_RHS_Ctrl, initial_condition, t_span, pll_problem)
  control_sol = solve(control_prob, Tsit5(), callback=all_cb_control, save_everystep=false, maxiters=8_500_000)
  push!(control_sols, control_sol)
end

function plot_amplitude_surface(ps::Array{S}, alphas::Array{Float64}; k=1) where {S<:Any}
  CairoMakie.activate!(type="svg")
  fig = Figure(
    size=(254, 300),
    backgroundcolor=:white,
    framevisible=false
  )
  ax_a = Axis3(fig[1, 1],
    xlabelvisible=false,
    ylabelvisible=false,
    zlabelvisible=false,
    xgridvisible=false,
    ygridvisible=false,
    zgridvisible=false,
    yreversed=false,
    xypanelvisible=false, yzpanelvisible=false, xzpanelvisible=false,
    azimuth=-1.0
  )
  for (idx, p) in enumerate(ps)
    omega_drivings = zeros(Float64, length(p.targets))
    amplitude_k = zeros(Float64, length(p.targets))
    omega_drivings .= p.matrix_harmonic_state_converged_control[1, :]
    amplitude_k .= norm.(p.matrix_harmonic_state_converged_control[2*k+2, :], p.matrix_harmonic_state_converged_control[2*k+1, :]) * 1.e9
    alpha = ones(Float64, length(p.targets)) .* alphas[idx] .* 1.e9
    scatterlines!(ax_a, omega_drivings, alpha, amplitude_k, color=(:black, 0.4), markersize=2.0)
  end
  return fig
end


function plot_phaselag_surface(ps::Array{S}, alphas::Array{Float64}; k=1) where {S<:Any}
  CairoMakie.activate!(type="svg")
  fig = Figure(
    size=(254, 300),
    backgroundcolor=:white,
    framevisible=false
  )
  ax_phi = Axis3(fig[1, 1],
    xlabelvisible=false,
    ylabelvisible=false,
    zlabelvisible=false,
    xgridvisible=false,
    ygridvisible=false,
    zgridvisible=false,
    yreversed=false,
    xypanelvisible=false, yzpanelvisible=false, xzpanelvisible=false,
    azimuth=-1.0
  )
  for (idx, p) in enumerate(ps)
    omega_drivings = zeros(Float64, length(p.targets))
    phaselag_k = zeros(Float64, length(p.targets))
    omega_drivings .= p.matrix_harmonic_state_converged_control[1, :]
    phaselag_k .= atan.(p.matrix_harmonic_state_converged_control[2*k+2, :], p.matrix_harmonic_state_converged_control[2*k+1, :])
    alpha = ones(Float64, length(p.targets)) .* alphas[idx] .* 1.e9
    scatterlines!(ax_phi, omega_drivings, alpha, phaselag_k, color=(:black, 0.4), markersize=2.0)
  end
  zlims!(ax_phi, (-π, 0.0))
  return fig
end

