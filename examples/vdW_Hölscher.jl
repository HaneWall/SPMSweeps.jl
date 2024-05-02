using SPMsweeps, StaticArrays, DifferentialEquations, CairoMakie, LinearAlgebra, ProgressLogging

initial_condition = SVector(0., 0., 0.)

# params from Hölscher and Schwartz
H = 2.e-19
a_0 = 0.3e-9
d = 10.e-9
R = 10.18e-9
Q = 300
k = 40
f = 1.33e-9




# time params sweep
Δ_t_filt = 0.05
Δ_t_sweep = 10_000.
μ_sweep = 0.05

# time params control
Δ_t_ctrl = 0.05 # timestep for each control step
μ_ctrl = 0.05 
Δ_t_checker = 4000. # timestep in which we check convergence 


μ = Δ_t_ctrl # stepsize LMS (tends to be equal to sample time of control)
harms = [1., 2., 3., 4.] # respected higher harmonics (DC always automatically included)
K_P = 1.
K_I = 0.001
K_D = 0.
τ = 2.
int_min = -0.01
int_max =  0.01
ctrl_min = -0.011
ctrl_max = 0.011

# for control 
CTRL = PID_Controller_Tustin(Δ_t_ctrl, K_P, K_I, K_D, τ, int_min, int_max, ctrl_min, ctrl_max)
#CTRL = PID_Controller_Euler_FWD(Δ_t_ctrl, K_P, K_I, K_D)
FILT = LMS_Algorithm(μ_ctrl, harms)
CHECK = Constant_Time_Check()

# for sweeps 
CTRL_fwd = PID_Controller_Euler_FWD()
FILT_fwd = LMS_Algorithm(μ_sweep, harms) 
CHECK_fwd = Constant_Time_Check()

CTRL_bwd = PID_Controller_Euler_FWD()
FILT_bwd = LMS_Algorithm(μ_sweep, harms) 
CHECK_bwd = Constant_Time_Check()

# targets for control (not relevant for sweeps)
TARGETS = collect(range(-0.5, -2.6, length=120))
err = 0.02

# frequncy array that we would like to sweep through (not relevant for control)
OMEGAS = collect(range(0.997, 1.003, length=60))

pll_problem = vdW_oscillator(k, f, OMEGAS[1], a_0, d, H, R, Q, OMEGAS, TARGETS, err, CTRL, FILT, CHECK)
fwd_problem = vdW_oscillator(k, f, OMEGAS[1], a_0, d, H, R, Q,  OMEGAS, TARGETS, err, CTRL_fwd, FILT_fwd, CHECK_fwd)
bwd_problem = vdW_oscillator(k, f, OMEGAS[end], a_0, d, H, R, Q,  reverse(OMEGAS), TARGETS, err, CTRL_fwd, FILT_fwd, CHECK_fwd)

# for control 
control_cb = PeriodicCallback(ctrl_cb!, Δ_t_ctrl)
convergence_cb = PeriodicCallback(conv_cb!, Δ_t_checker)
all_cb_control = CallbackSet(control_cb, convergence_cb)

# for sweep
sweep_cb = PeriodicCallback(freq_sweep_cb!, Δ_t_sweep)
filter_cb = PeriodicCallback(lms_cb!, Δ_t_filt)
all_cb_sweep = CallbackSet(sweep_cb, filter_cb)


t_span = (0, 800_000)

sweep_fwd_prob = ODEProblem(f_RHS, initial_condition, t_span, fwd_problem)
sweep_fwd_sol = solve(sweep_fwd_prob, callback=all_cb_sweep, save_everystep=false, maxiters=100_000_000)


## sweep backward
sweep_bwd_prob = ODEProblem(f_RHS, initial_condition, t_span, bwd_problem)
sweep_bwd_sol = solve(sweep_bwd_prob, callback=all_cb_sweep, save_everystep=false, maxiters=100_000_000)


## control sweep 
control_prob = ODEProblem(f_RHS_Ctrl, initial_condition, (0, 600_000), pll_problem)
control_sol = solve(control_prob, callback=all_cb_control, save_everystep=false, maxiters=20_000_000)



function plot_sweeps(ps::Array{S}; k=1) where S <: Nanojunction
    CairoMakie.activate!(type="svg")
    fig = Figure(
        size = (508, 258), 
        fonts = (; regular = "Computer Modern"), 
        fontsize=12, 
        backgroundcolor=:white
    )
    ax_a = Axis(fig[1, 1],
        xlabel=L"\omega_d/\omega_0", 
        ylabel=L"amplitude $k$th harmonic $A_k$", 
        xgridvisible=false,
        ygridvisible=false)
    ax_phi = Axis(fig[1, 2],
        xlabel=L"\omega_d/\omega_0", 
        ylabel=L"phaselag $k$th harmonic $\phi_k$", 
        xgridvisible=false,
        ygridvisible=false)
    for p in ps
        omega_drivings_fwd = zeros(Float64, length(p.omegas))
        phaselag_1 = zeros(Float64, length(p.omegas))
        amplitude_1 = zeros(Float64, length(p.omegas))
        omega_drivings_fwd .= p.matrix_harmonic_state_converged_sweep[1, :]
        phaselag_1 .= atan.(p.matrix_harmonic_state_converged_sweep[2*k+2, :], p.matrix_harmonic_state_converged_sweep[2*k+1, :])
        amplitude_1 .= norm.(p.matrix_harmonic_state_converged_sweep[2*k+2, :], p.matrix_harmonic_state_converged_sweep[2*k+1, :])
        scatterlines!(ax_a, omega_drivings_fwd[1:end], amplitude_1[1:end])
        scatterlines!(ax_phi, omega_drivings_fwd[1:end], phaselag_1[1:end]) 
    end
    fig
  end

  function plot_control(ps::Array{S}; k=1) where S <: Nanojunction
    CairoMakie.activate!(type="svg")
    fig = Figure(
        size = (508, 258), 
        fonts = (; regular = "Computer Modern"), 
        fontsize=12, 
        backgroundcolor=:white
    )
    ax_a = Axis(fig[1, 1],
        xlabel=L"\omega_d/\omega_0", 
        ylabel=L"amplitude $k$th harmonic $A_k$", 
        xgridvisible=false,
        ygridvisible=false)
    ax_phi = Axis(fig[1, 2],
        xlabel=L"\omega_d/\omega_0", 
        ylabel=L"phaselag $k$th harmonic $\phi_k$", 
        xgridvisible=false,
        ygridvisible=false)
  
    for p in ps
        omega_drivings = zeros(Float64, length(p.targets))
        phaselag_k = zeros(Float64, length(p.targets))
        amplitude_k = zeros(Float64, length(p.targets))
        omega_drivings .= p.matrix_harmonic_state_converged_control[1, :]
        phaselag_k .= atan.(p.matrix_harmonic_state_converged_control[2*k+2, :], p.matrix_harmonic_state_converged_control[2*k+1, :])
        amplitude_k .= norm.(p.matrix_harmonic_state_converged_control[2*k+2, :], p.matrix_harmonic_state_converged_control[2*k+1, :])
        scatterlines!(ax_a, omega_drivings[1:end], amplitude_k[1:end])
        scatterlines!(ax_phi, omega_drivings[1:end], phaselag_k[1:end]) 
    end
    fig
  end
  
