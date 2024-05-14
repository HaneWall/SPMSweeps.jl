using SPMsweeps, StaticArrays, DifferentialEquations, CairoMakie, LinearAlgebra

initial_condition = SVector(0., 0., 0.)

# params duffing 
k_1 = 1. 
k_2 = 0.01
k_3s = collect(range(-1., 0., length=20))
Ω_s = 0.97
f = 0.0025


# time params sweep
Δ_t_filt = 0.05
Δ_t_sweep = 3_000.
μ_sweep = 0.05

# time params control
Δ_t_ctrl = 0.05 # timestep for each control step
μ_ctrl = 0.05 
Δ_t_checker = 800. # timestep in which we check convergence 


μ = Δ_t_ctrl # stepsize LMS (tends to be equal to sample time of control)
harms = [1., 2., 3., 4.] # respected higher harmonics (DC always automatically included)
K_P = 1.
K_I = 0.04
K_D = 2.
τ = 2.
int_min = -0.5
int_max =  0.5
ctrl_min = -0.6
ctrl_max = 0.6

# for control 
#CTRL = PID_Controller_Tustin(Δ_t_ctrl, K_P, K_I, K_D, τ, int_min, int_max, ctrl_min, ctrl_max)
#CTRL = PID_Controller_Euler_FWD(Δ_t_ctrl, K_P, K_I, K_D)
#FILT = LMS_Algorithm(μ_ctrl, harms)
#CHECK = Constant_Time_Check()

# for sweeps 
CTRL_fwd = PID_Controller_Euler_FWD()
FILT_fwd = LMS_Algorithm(μ_sweep, harms) 
CHECK_fwd = Constant_Time_Check()

CTRL_bwd = PID_Controller_Euler_FWD()
FILT_bwd = LMS_Algorithm(μ_sweep, harms) 
CHECK_bwd = Constant_Time_Check()

# targets for control (not relevant for sweeps)
TARGETS = collect(range(-0.1, -3., length=50))
err = 0.02

# frequncy array that we would like to sweep through (not relevant for control)
OMEGAS = collect(range(0.97, 1.03, length=40))

#fwd_problem = Duffing_oscillator(k_1, k_2, k_3, f, OMEGAS[1], OMEGAS, TARGETS, err, CTRL_fwd, FILT_fwd, CHECK_fwd)
#bwd_problem = Duffing_oscillator(k_1, k_2, k_3, f, OMEGAS[end], reverse(OMEGAS), TARGETS, err, CTRL_bwd, FILT_bwd, CHECK_bwd)

pll_problems = []

for k_3 in k_3s
    CTRL = PID_Controller_Tustin(Δ_t_ctrl, K_P, K_I, K_D, τ, int_min, int_max, ctrl_min, ctrl_max)
    #CTRL = PID_Controller_Euler_FWD(Δ_t_ctrl, K_P, K_I, K_D)
    FILT = LMS_Algorithm(μ_ctrl, harms)
    CHECK = Constant_Time_Check()
    push!(pll_problems, Duffing_oscillator(k_1, k_2, k_3, f, OMEGAS[1], OMEGAS, TARGETS, err, CTRL, FILT, CHECK))
end

# for control 
control_cb = PeriodicCallback(ctrl_cb!, Δ_t_ctrl)
convergence_cb = PeriodicCallback(conv_cb!, Δ_t_checker)
all_cb_control = CallbackSet(control_cb, convergence_cb)

# for sweep
sweep_cb = PeriodicCallback(freq_sweep_cb!, Δ_t_sweep)
filter_cb = PeriodicCallback(lms_cb!, Δ_t_filt)
all_cb_sweep = CallbackSet(sweep_cb, filter_cb)

t_span = (0, 300_000)

## sweep forward 
#sweep_fwd_prob = ODEProblem(f_RHS, initial_condition, t_span, fwd_problem)
#sweep_fwd_sol = solve(sweep_fwd_prob, callback=all_cb_sweep, save_everystep=false, maxiters=3_000_000)


## sweep backward
#sweep_bwd_prob = ODEProblem(f_RHS, initial_condition, t_span, bwd_problem)
#sweep_bwd_sol = solve(sweep_bwd_prob, callback=all_cb_sweep, save_everystep=false, maxiters=3_000_000)


## control sweep 
control_sols = []

for pll_problem in pll_problems
    control_prob = ODEProblem(f_RHS_Ctrl, initial_condition, (0, 200_000), pll_problem)
    push!(control_sols, solve(control_prob, callback=all_cb_control, save_everystep=false, maxiters=20_000_000))
end


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
  
function plot_control(ps::Array{S}; k=1) where S<:Any
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
        scatterlines!(ax_a, omega_drivings[1:end], amplitude_k[1:end], markersize=4.)
        scatterlines!(ax_phi, omega_drivings[1:end], phaselag_k[1:end], markersize=4.) 
    end
    fig
end
  

function plot_response_surface(ps::Array{S}, alphas::Array{Float64}; k=1) where S <: Any
    CairoMakie.activate!(type="svg")
    fig = Figure(
        size = (508, 300), 
        fonts = (; regular = "Computer Modern"), 
        fontsize=12, 
        backgroundcolor=:white,
        framevisible=false
    ) 
    ax_a = Axis3(fig[1, 1],
        xlabel=L"\omega_d/\omega_0", 
        ylabel=L"nonlinear coeffcient $\beta$",
        zlabel=L"amplitude first harmonic $A_1$",   
        xgridvisible = false, 
        ygridvisible = false, 
        zgridvisible = false,
        yreversed=false,
        xypanelvisible=false, yzpanelvisible=false, xzpanelvisible=false,
        azimuth=-1.
    )

    ax_phi = Axis3(fig[1, 2],
        xlabel=L"\omega_d/\omega_0", 
        ylabel=L"nonlinear coeffcient $\beta$",
        zlabel=L"phaselag first harmonic $\phi_1$",   
        xgridvisible = false, 
        ygridvisible = false, 
        zgridvisible = false,
        yreversed=false,
        xypanelvisible=false, yzpanelvisible=false, xzpanelvisible=false,
        azimuth=-1.
    ) 
    for (idx, p) in enumerate(ps)
        omega_drivings = zeros(Float64, length(p.targets))
        phaselag_k = zeros(Float64, length(p.targets))
        amplitude_k = zeros(Float64, length(p.targets))
        omega_drivings .= p.matrix_harmonic_state_converged_control[1, :]
        phaselag_k .= atan.(p.matrix_harmonic_state_converged_control[2*k+2, :], p.matrix_harmonic_state_converged_control[2*k+1, :])
        amplitude_k .= norm.(p.matrix_harmonic_state_converged_control[2*k+2, :], p.matrix_harmonic_state_converged_control[2*k+1, :])
        
        alpha = ones(Float64, length(p.targets)) .* alphas[idx]
        scatterlines!(ax_a, omega_drivings, alpha, amplitude_k, color=(:navy,0.4), markersize=3, linewidth=1.5)
        scatterlines!(ax_phi, omega_drivings, alpha, phaselag_k, color=(:navy, 0.4), markersize=3, linewidth=1.5)
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
        ylabel=L"amplitude first harmonic $A_1$[nm]", 
        xgridvisible=false,
        ygridvisible=false)
    ax_phi = Axis(fig[1, 2],
        xlabel=L"\omega_d/\omega_0", 
        ylabel=L"phaselag first harmonic $\phi_1$[rad]", 
        xgridvisible=false,
        ygridvisible=false)
  
    for p in ps
        omega_drivings = zeros(Float64, length(p.targets))
        phaselag_k = zeros(Float64, length(p.targets))
        amplitude_k = zeros(Float64, length(p.targets))
        omega_drivings .= p.matrix_harmonic_state_converged_control[1, :]
        phaselag_k .= atan.(p.matrix_harmonic_state_converged_control[2*k+2, :], p.matrix_harmonic_state_converged_control[2*k+1, :])
        amplitude_k .= norm.(p.matrix_harmonic_state_converged_control[2*k+2, :], p.matrix_harmonic_state_converged_control[2*k+1, :]) .* 1e9
        scatterlines!(ax_a, omega_drivings[1:end], amplitude_k[1:end], markersize=3, color=(:navy, 0.4))
        scatterlines!(ax_phi, omega_drivings[1:end], phaselag_k[1:end], markersize=3, color=(:navy, 0.4)) 
    end
    fig
  end
  