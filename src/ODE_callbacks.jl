
"""
This will be a periodic callback, that uses the LMS-algorithm to estimate the current harmonic state every sample timestep
"""
function lms_cb!(integrator)
  step_demodulation!(integrator.p.filter, integrator.u[1], integrator.u[3])
end

"""
This will be a periodic callback, that is applied every sample timestep, LMS is alrdy included and should not be applied again
"""
function ctrl_cb!(integrator)
  # here apply LMS to calcualte ϕ from measurement x[1], that is integrator.u[1]
  step_demodulation!(integrator.p.filter, integrator.u[1], integrator.u[3])
  ϕ = atan(integrator.p.filter.harmonic_state[3], integrator.p.filter.harmonic_state[2])

  step_controller!(integrator.p.ctrl, integrator.p.targets[integrator.p.target_idx], ϕ)
end


"""
This will be a periodic callback, that is applied to
 save some outputs for plots
"""
function saving_cb!(integrator)
  nothing
end

"""
This will be a periodic callback, that is applied to check convergence towards steady state and furthermore stores converged results
"""
function conv_cb!(integrator)
  # for now we use phaselag of first harmonic to ensure steady stateness
  ϕ = atan(integrator.p.filter.harmonic_state[3], integrator.p.filter.harmonic_state[2])
  check_converged!(integrator.p.checker, ϕ)
  if integrator.p.checker.converged && abs(integrator.p.targets[integrator.p.target_idx] - ϕ) < integrator.p.target_err
    #save to matrix with converged results
    integrator.p.matrix_harmonic_state_converged_control[1, integrator.p.target_idx] = integrator.p.Ω_s + integrator.p.ctrl.output
    integrator.p.matrix_harmonic_state_converged_control[2:end, integrator.p.target_idx] .= integrator.p.filter.harmonic_state
    integrator.p.target_idx += 1
  end
  
  if integrator.p.target_idx > length(integrator.p.targets)
    terminate!(integrator)
  end  

end

"""
This will be periodic callback, that changes the current angular frequncy for sweeps 
"""
function freq_sweep_cb!(integrator)
  ϕ = atan(integrator.p.filter.harmonic_state[3], integrator.p.filter.harmonic_state[2])
  check_converged!(integrator.p.checker, ϕ)

  if integrator.p.checker.converged 
    integrator.p.matrix_harmonic_state_converged_sweep[1, integrator.p.target_idx] = integrator.p.Ω_s
    integrator.p.matrix_harmonic_state_converged_sweep[2:end, integrator.p.target_idx] .= integrator.p.filter.harmonic_state
    integrator.p.target_idx += 1
    if integrator.p.target_idx > length(integrator.p.omegas)
      terminate!(integrator)
    else
      integrator.p.Ω_s = integrator.p.omegas[integrator.p.target_idx]
    end
  end
end
