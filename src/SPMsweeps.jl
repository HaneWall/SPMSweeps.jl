module SPMsweeps
# Write your package code here.
using DifferentialEquations, StaticArrays, DataStructures, LinearAlgebra

include("utilities.jl")
export create_bases_psi, error_lms, update_basis!, update_weights!, welford_cb!

include("controller.jl")
export Controller, PID_Controller_Euler_FWD, PID_Controller_Tustin, step_controller!

include("ODE_callbacks.jl")
export ctrl_cb!, freq_sweep_cb!, conv_cb!, lms_cb!, saving_cb!

include("adaptive_filter.jl")
export AdaptiveFilter, LMS_Algorithm, step_demodulation!

include("steadystate.jl")
export SteadyStateChecker, Constant_Time_Check, Welford_Buffer, check_converged!

include("nanojunctions.jl")
export f_RHS, f_RHS_Ctrl, Nanojunction, Duffing_oscillator, vdW_oscillator, DMT_oscillator

end
