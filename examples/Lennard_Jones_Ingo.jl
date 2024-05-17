using SPMsweeps
using StaticArrays, DifferentialEquations, CairoMakie, LinearAlgebra

initial_condition = SVector(0., 0., 0.)

# params from Ingo sims / Annas dissertation / my master thesis
σ = 2.8e-9             
V_0 = 4.2e-18
# softening parameter
δx = 0.5e-9
# distance betweeen cantilever and surface/sample
d = 60.e-9
# Quality-factor of cantilever
Q = 400.
# spring constant of cantilever
k = 0.7
# resonance frequency
ω_0 =  50.e3 * 2π
# friction parameter
γ = 1.e-31
#forcing amplitude
f= 170.e-9 / Q 




# time params sweep
Δ_t_filt = 0.05
Δ_t_sweep = 10_000.
μ_sweep = 0.05

# time params control
Δ_t_ctrl = 0.08 # timestep for each control step
μ_ctrl = 0.05 
Δ_t_checker = 15.3 # timestep in which we check convergence 


μ = Δ_t_ctrl # stepsize LMS (tends to be equal to sample time of control)
harms = [1., 2., 3., 4.] # respected higher harmonics (DC always automatically included)
K_P = 0.1
K_I = 0.001
K_D = 1.5
τ = 5.
int_min = -0.05
int_max =  0.05
ctrl_min = -0.051
ctrl_max = 0.051

# for control 
CTRL = PID_Controller_Tustin(Δ_t_ctrl, K_P, K_I, K_D, τ, int_min, int_max, ctrl_min, ctrl_max)
#CTRL = PID_Controller_Euler_FWD(Δ_t_ctrl, K_P, K_I, K_D)
FILT = LMS_Algorithm(μ_ctrl, harms)
CHECK = Welford_Buffer(40, 1e-3)

# for sweeps 
CTRL_fwd = PID_Controller_Euler_FWD()
FILT_fwd = LMS_Algorithm(μ_sweep, harms) 
CHECK_fwd = Constant_Time_Check()

CTRL_bwd = PID_Controller_Euler_FWD()
FILT_bwd = LMS_Algorithm(μ_sweep, harms) 
CHECK_bwd = Constant_Time_Check()

# targets for control (not relevant for sweeps)
TARGETS = collect(range(-0.2, -2.8, length=50))
err = 0.02

# frequncy array that we would like to sweep through (not relevant for control)
OMEGAS = collect(range(0.985, 1.015, length=60))

pll_problem = Lennard_Jones_oscillator(k, f, V_0, γ, OMEGAS[1], d, σ, Q, δx, ω_0, OMEGAS, TARGETS, err, CTRL, FILT, CHECK)
fwd_problem = Lennard_Jones_oscillator(k, f, V_0, γ, OMEGAS[1], d, σ, Q, δx, ω_0, OMEGAS, TARGETS, err, CTRL_fwd, FILT_fwd, CHECK_fwd)
bwd_problem = Lennard_Jones_oscillator(k, f, V_0, γ, OMEGAS[end], d, σ, Q, δx, ω_0, reverse(OMEGAS), TARGETS, err, CTRL_fwd, FILT_fwd, CHECK_fwd)

# for control 
control_cb = PeriodicCallback(ctrl_cb!, Δ_t_ctrl)
convergence_cb = PeriodicCallback(conv_wf_cb!, Δ_t_checker)
all_cb_control = CallbackSet(control_cb, convergence_cb)

# for sweep
sweep_cb = PeriodicCallback(freq_sweep_cb!, Δ_t_sweep)
filter_cb = PeriodicCallback(lms_cb!, Δ_t_filt)
all_cb_sweep = CallbackSet(sweep_cb, filter_cb)

t_span = (0, 800_000)

## sweep forward 
sweep_fwd_prob = ODEProblem(f_RHS, initial_condition, t_span, fwd_problem)
@time sweep_fwd_sol = solve(sweep_fwd_prob, Tsit5(), callback=all_cb_sweep, save_everystep=false, maxiters=100_000_000)


## sweep backward
sweep_bwd_prob = ODEProblem(f_RHS, initial_condition, t_span, bwd_problem)
sweep_bwd_sol = solve(sweep_bwd_prob, Tsit5(), callback=all_cb_sweep, save_everystep=false, maxiters=100_000_000)


## control sweep 
control_prob = ODEProblem(f_RHS_Ctrl, initial_condition, (0, 200_000), pll_problem)
@time control_sol = solve(control_prob, Tsit5(), callback=all_cb_control, save_everystep=false, maxiters=8_500_000)
