using SPMsweeps, StaticArrays, DifferentialEquations, CairoMakie, LinearAlgebra

initial_condition = SVector(0., 0., 0.)

# params from Hölscher and Schwartz
H = 2.e-19
a_0 = 0.3e-9
ds = collect(range(10.e-9, 10.5e-9, length=10))
R= 10.e-9
Q = 300
k = 40
f = 1.33e-9

E_s = 1.e9
nu_s = 0.3

E_t = 130.e9
nu_t = 0.3


# time params control
Δ_t_ctrl = 0.05 # timestep for each control step
μ_ctrl = 0.05 
Δ_t_checker = 700. # timestep in which we check convergence 


# targets for control (not relevant for sweeps)
TARGETS = collect(range(-0.5, -2.6, length=120))
err = 0.02

# frequncy array that we would like to sweep through (not relevant for control)
OMEGAS = collect(range(0.997, 1.003, length=60))


μ = Δ_t_ctrl # stepsize LMS (tends to be equal to sample time of control)
harms = [1., 2., 3., 4.] # respected higher harmonics (DC always automatically included)
K_P = 1.6
K_I = 0.001
K_D = 0.5
τ = 2.
int_min = -0.01
int_max =  0.01
ctrl_min = -0.011
ctrl_max = 0.011



pll_problems = []

for d in ds
    CTRL = PID_Controller_Tustin(Δ_t_ctrl, K_P, K_I, K_D, τ, int_min, int_max, ctrl_min, ctrl_max)
    FILT = LMS_Algorithm(μ_ctrl, harms)
    CHECK = Constant_Time_Check()
    push!(pll_problems, DMT_oscillator(k, f, OMEGAS[1], a_0, d, H, R, Q, nu_t, nu_s, E_t, E_s, OMEGAS, TARGETS, err, CTRL, FILT, CHECK))
end


# for control 
control_cb = PeriodicCallback(ctrl_cb!, Δ_t_ctrl)
convergence_cb = PeriodicCallback(conv_cb!, Δ_t_checker)
all_cb_control = CallbackSet(control_cb, convergence_cb)

control_sols = []

for pll_problem in pll_problems
    control_prob = ODEProblem(f_RHS_Ctrl, initial_condition, (0, 600_000), pll_problem)
    push!(control_sols, solve(control_prob, callback=all_cb_control, save_everystep=false, maxiters=20_000_000, progress=true))
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




function plot_response_surface(ps::Array{S}, alphas::Array{Float64}; k=1) where S <: Any
    CairoMakie.activate!(type="svg")
    fig = Figure(
        size = (508, 700), 
        fonts = (; regular = "Computer Modern"), 
        fontsize=12, 
        backgroundcolor=:white,
        framevisible=false
    ) 
    ax_a = Axis3(fig[1, 1],
        xlabel=L"\omega_d/\omega_0", 
        ylabel=L"seperation distance $d$[nm]",
        zlabel=L"amplitude first harmonic $A_1$[nm]",   
        xgridvisible = false, 
        ygridvisible = false, 
        zgridvisible = false,
        yreversed=false,
        xypanelvisible=false, yzpanelvisible=false, xzpanelvisible=false,
        azimuth=-1.
    )

    ax_phi = Axis3(fig[2, 1],
        xlabel=L"\omega_d/\omega_0", 
        ylabel=L"seperation distance $d$[nm]",
        zlabel=L"phaselag first harmonic $\phi_1$[rad]",   
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
        amplitude_k .= norm.(p.matrix_harmonic_state_converged_control[2*k+2, :], p.matrix_harmonic_state_converged_control[2*k+1, :]) * 1.e9
        
        alpha = ones(Float64, length(p.targets)) .* alphas[idx] .* 1.e9
        scatterlines!(ax_a, omega_drivings, alpha, amplitude_k, color=(:navy,0.4), markersize=3, linewidth=1.5)
        scatterlines!(ax_phi, omega_drivings, alpha, phaselag_k, color=(:navy, 0.4), markersize=3, linewidth=1.5)
    end
    fig
end