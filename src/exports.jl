export
  #utils
  create_bases_psi, error_lms, update_basis!, update_weights!,
  #controller 
  Controller, PID_Controller_Euler_FWD, PID_Controller_Tustin, step_controller!,
  #callbacks
  ctrl_cb!, freq_sweep_cb!, conv_cb!, lms_cb!, saving_cb_control!, conv_wf_cb!,
  #AdaptiveFilter
  AdaptiveFilter, LMS_Algorithm, step_demodulation!,
  #SteadyStateChecker
  SteadyStateChecker, Constant_Time_Check, Welford_Buffer, check_converged!,
  #Right hand side
  f_RHS, f_RHS_Ctrl, Nanojunction, Duffing_oscillator, vdW_oscillator, DMT_oscillator, Lennard_Jones_oscillator,
  #colors
  generate_cmap, theme!, plot_control, plot_sweeps, plot_sweep, plot_sweeps_control, save_fig
