include("colorscheme.jl")
using CairoMakie

function generate_cmap(n)
    if n > length(COLORS)
        return :viridis
    else
        return cgrad(COLORS[1:n], n; categorical = true)
    end
end

function new_theme!()
    set_theme!(;
        palette = (color = COLORS,), 
        fonts = (; regular = "Helvetica"),
        pt_per_unit = 1,
        fontsize = 8,
    )
end


function plot_control(ps::Array{S}; k=1) where S <: Any
    n_counter = 1
    CairoMakie.activate!(type="svg")
    fig = Figure(
        size = (508, 258), 
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
        scatterlines!(ax_a, omega_drivings[1:end], amplitude_k[1:end], marker = CMARKERS[n_counter], markersize=CMARKERSIZES[n_counter])
        scatterlines!(ax_phi, omega_drivings[1:end], phaselag_k[1:end], marker = CMARKERS[n_counter], markersize=CMARKERSIZES[n_counter])
        n_counter =+ 1 
    end
    fig
end

function plot_control(p::Nanojunction; k=1)
    plot_control([p]; k)
end



function plot_sweeps(ps::Array{S}; k=1) where S <: Nanojunction
    n_counter=1
    CairoMakie.activate!(type="svg")
    fig = Figure(
        size = (508, 258), 
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
        scatterlines!(ax_a, omega_drivings_fwd[1:end], amplitude_1[1:end],marker = CMARKERS[n_counter], markersize=CMARKERSIZES[n_counter])
        scatterlines!(ax_phi, omega_drivings_fwd[1:end], phaselag_1[1:end],marker = CMARKERS[n_counter], markersize=CMARKERSIZES[n_counter]) 
        n_counter += 1
    end
    fig
  end


function plot_sweep(ps::S; k=1) where S <: Nanojunction
   plot_sweeps([ps]; k)
end

function plot_sweeps_control(sweeps::Array{S}, controls::Array{S}; k=1) where S <: Nanojunction
    n_counter = 1
    CairoMakie.activate!(type="svg")
    fig = Figure(
        fontsize=7,
        size = (125, 220), 
        ppt=300,
        pt_per_unit = 1,
    )
    ax_a = Axis(fig[1, 1],
        #xlabel=L"\omega_d/\omega_0", 
        #ylabel=L"amplitude $A_k$", 
        xgridvisible=false,
        ygridvisible=false)
    ax_phi = Axis(fig[2, 1],
        xlabel=L"\omega_d/\omega_0", 
        #ylabel=L"phaselag $\phi_k$", 
        xgridvisible=false,
        ygridvisible=false)
    for p in sweeps
        omega_drivings_fwd = zeros(Float64, length(p.omegas))
        phaselag_1 = zeros(Float64, length(p.omegas))
        amplitude_1 = zeros(Float64, length(p.omegas))
        omega_drivings_fwd .= p.matrix_harmonic_state_converged_sweep[1, :]
        phaselag_1 .= atan.(p.matrix_harmonic_state_converged_sweep[2*k+2, :], p.matrix_harmonic_state_converged_sweep[2*k+1, :])
        amplitude_1 .= norm.(p.matrix_harmonic_state_converged_sweep[2*k+2, :], p.matrix_harmonic_state_converged_sweep[2*k+1, :])
        scatter!(ax_a, omega_drivings_fwd[1:end], amplitude_1[1:end],marker = CMARKERS[n_counter], markersize=CMARKERSIZES[n_counter])
        scatter!(ax_phi, omega_drivings_fwd[1:end], phaselag_1[1:end],marker = CMARKERS[n_counter], markersize=CMARKERSIZES[n_counter]) 
        n_counter += 1
    end

    for c in controls
        omega_drivings = zeros(Float64, length(c.targets))
        phaselag_k = zeros(Float64, length(c.targets))
        amplitude_k = zeros(Float64, length(c.targets))
        omega_drivings .= c.matrix_harmonic_state_converged_control[1, :]
        phaselag_k .= atan.(c.matrix_harmonic_state_converged_control[2*k+2, :], c.matrix_harmonic_state_converged_control[2*k+1, :])
        amplitude_k .= norm.(c.matrix_harmonic_state_converged_control[2*k+2, :], c.matrix_harmonic_state_converged_control[2*k+1, :])
        scatter!(ax_a, omega_drivings[1:end], amplitude_k[1:end],marker =CMARKERS[n_counter], markersize=CMARKERSIZES[n_counter])
        scatter!(ax_phi, omega_drivings[1:end], phaselag_k[1:end],marker = CMARKERS[n_counter], markersize=CMARKERSIZES[n_counter]) 
        n_counter+=1
    end
    return fig
end
