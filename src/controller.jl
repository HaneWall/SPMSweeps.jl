abstract type Controller end

"""
Discrete PID controller, that can be obtained by using Tustin-Transformation. In contrast 
to Euler-forward or Euler-backward implementations, this discete PID essentially arrises 
by trapezoidal approximations.

PID controller, that contains Anti-Windup for the Integralterm and furthermore maximum and minimum output terms. 


"""
mutable struct PID_Controller_Tustin <: Controller
    const Δt::Float64
    const c_1::Float64
    const c_2::Float64
    const c_3::Float64
    const c_4::Float64

    # derivative low-pass filer time constant 
    const τ::Float64

    # Integrator and Differentiator values 
    proportional::Float64
    integral::Float64
    differential::Float64

    prev_error::Float64
    prev_measurement::Float64

    # anti-windup for integrator, maximal and minimal values that the integrator is allowed to accomplish
    const int_max::Float64
    const int_min::Float64

    # maximal and minimal outputs, that the controller is allowed to accomplish
    const output_max::Float64
    const output_min::Float64

    # output of the controller
    output::Float64

    function PID_Controller_Tustin(sample_time, K_P, K_I, K_D, time_const_low_pass, int_min, int_max, out_min, out_max)
        new(
            sample_time,
            K_P,
            K_I * sample_time / 2.0,
            2 * K_D / (2 * time_const_low_pass + sample_time),
            (2 * time_const_low_pass - sample_time) / (2 * time_const_low_pass + sample_time),
            time_const_low_pass,
            0.0,
            0.0,
            0.0,
            0.0,
            0.0,
            int_max,
            int_min,
            out_max,
            out_min,
            0.0
        )

    end
end

function step_controller!(p::PID_Controller_Tustin, target::Float64, measurement::Float64)
    error = measurement - target
    p.proportional = p.c_1 * error
    p.integral = p.integral + p.c_2 * (error + p.prev_error)


    p.integral = clamp(p.integral, p.int_min, p.int_max)


    p.differential = p.c_3 * (error - p.prev_error) + p.c_4 * p.differential
    p.output = p.proportional + p.integral + p.differential

    p.output = clamp(p.output, p.output_min, p.output_max)

    p.prev_error = error
end


mutable struct PID_Controller_Euler_FWD <: Controller
    const Δt::Float64
    const c_1::Float64
    const c_2::Float64
    const c_3::Float64
    storage::Array{Float64}
    output::Float64

    function PID_Controller_Euler_FWD(time_step, K_P, K_I, K_D)
        new(time_step, # Δt
            K_P + K_I * time_step + K_D / time_step,  # c_1
            -K_P - 2 * K_D / time_step, # c_2
            K_D / time_step, # c_3
            zeros(Float64, 3),
            0.0
        )
    end

    function PID_Controller_Euler_FWD()
        PID_Controller_Euler_FWD(1.0, 0.0, 0.0, 0.0)
    end
end

function step_controller!(p::PID_Controller_Euler_FWD, target::Float64, measurement::Float64)
    p.storage[3] = p.storage[2]
    p.storage[2] = p.storage[1]
    p.storage[1] = measurement - target
    p.output += p.c_1 * p.storage[1] + p.c_2 * p.storage[2] + p.c_3 * p.storage[3]
end


