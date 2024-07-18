using CSV
using DelimitedFiles
using DataFrames
using CairoMakie

# read data from experimental data folder

# controlled sweeps

freq_control = CSV.read("./experimental_data/2DKontrolle/Control_f_Omega_kHz.csv", DataFrame
  ; header=false)

error_freq_control = CSV.read("./experimental_data/2DKontrolle/error_f_kHz.csv", DataFrame
  ; header=false)

phi_control = CSV.read("./experimental_data/2DKontrolle/Control_phi_Omega-1_125.csv", DataFrame
  ; header=false)

error_phi_control = CSV.read("./experimental_data/2DKontrolle/error_phi_kHz.csv", DataFrame
  ; header=false)
# open-loop sweeps

freq_forward = CSV.read("./experimental_data/2DKontrolle/fpre_f_OMEGA_kHz.csv", DataFrame
  ; header=false)

freq_backward = CSV.read("./experimental_data/2DKontrolle/bpre_f_OMEGA_kHz.csv", DataFrame
  ; header=false)

phi_forward = CSV.read("./experimental_data/2DKontrolle/ϕ_fpre-1_125.csv", DataFrame
  ; header=false)

phi_backward = CSV.read("./experimental_data/2DKontrolle/ϕ_bpre-1_125.csv", DataFrame
  ; header=false)

# transient frequency plot
data_range = 1:70000
data_range_long = 1:300000

phi_control_transient = CSV.read("./experimental_data/2DKontrolle/ϕ_Control_long-1_125.csv", DataFrame
  ; header=false)

phi_control_target = CSV.read("./experimental_data/2DKontrolle/ϕ_target_long_-1_125.csv", DataFrame
  ; header=false)


# three dimensional plots (phase and amplitude plots)
setpoints = CSV.read("./experimental_data/Serviette/SP.csv", DataFrame; header=false)
frequencies = CSV.read("./experimental_data/Serviette/f.csv", DataFrame; header=false)
frequencies_sg = CSV.read("./experimental_data/Serviette/SG_f.csv", DataFrame; header=false)
phases = CSV.read("./experimental_data/Serviette/ϕ.csv", DataFrame; header=false)
amplitudes = CSV.read("./experimental_data/Serviette/R.csv", DataFrame; header=false)


function plot_exp_sweep_control_amplitude()
  CairoMakie.activate!(type="svg", pt_per_unit=1)
  fig = Figure(
    size=(150, 200),
  )
  MARKERS = [:rect, :circle, :circle]
  MARKERSIZES = [2, 1, 3]

  ax_a = Axis(fig[1, 1],
    xlabelvisible=false,
    ylabelvisible=false,
    xgridvisible=false,
    ygridvisible=false,
    xticks=collect(55.89:0.04:55.98))

  scatterlines!(ax_a, freq_forward[:, 1], phi_forward[:, 1], linewidth=MARKERSIZES[1], markersize=MARKERSIZES[1])
  scatterlines!(ax_a, freq_backward[:, 1], phi_backward[:, 1], linewidth=MARKERSIZES[2], markersize=MARKERSIZES[2])
  errorbars!(ax_a, freq_control[:, 1], phi_control[:, 1], error_freq_control[:, 1], whiskerwidth=1, direction=:x, color=:black)
  scatter!(ax_a, freq_control[:, 1], phi_control[:, 1], marker=MARKERS[3], markersize=MARKERSIZES[3], color=:black)
  ylims!(ax_a, (0.5 - 1.125, 2.2 - 1.125))
  xlims!(ax_a, (55.879, 55.981))
  return fig
end

function global_frequency_plot()
  CairoMakie.activate!(type="svg", pt_per_unit=1)
  fig = Figure(
    size=(135, 100),
  )
  MARKERS = [:rect, :circle, :circle]
  MARKERSIZES = [2, 1, 3]
  ax_a = Axis(fig[1, 1],
    xlabelvisible=false,
    ylabelvisible=false,
    xgridvisible=false,
    ygridvisible=false)

end

function phase_three_dimensional()
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
  zlims!(ax_phi, (-1.9, -0.3))
  ylims!(ax_phi, (80, 110))
  for i in 1:19
    lines!(ax_phi, frequencies_sg[(i*27-27+1):(i*27), 1], setpoints[(i*27-27+1):(i*27), 1], phases[(i*27-27+1):(i*27), 1], color=(:black, 0.4))
    scatter!(ax_phi, frequencies[(i*27-27+1):(i*27), 1], setpoints[(i*27-27+1):(i*27), 1], phases[(i*27-27+1):(i*27), 1], markersize=2.0, color=(:black, 0.4))
  end
  return fig
end


function amplitude_three_dimensional()
  CairoMakie.activate!(type="svg")
  fig = Figure(
    size=(254, 300),
    backgroundcolor=:white,
    framevisible=false
  )
  ax_r = Axis3(fig[1, 1],
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
  zlims!(ax_r, (65, 115))
  ylims!(ax_r, (80, 110))
  for i in 1:19
    lines!(ax_r, frequencies_sg[(i*27-27+1):(i*27), 1], setpoints[(i*27-27+1):(i*27), 1], amplitudes[(i*27-27+1):(i*27), 1], color=(:black, 0.4))
    scatter!(ax_r, frequencies[(i*27-27+1):(i*27), 1], setpoints[(i*27-27+1):(i*27), 1], amplitudes[(i*27-27+1):(i*27), 1], markersize=2.0, color=(:black, 0.4))
  end
  return fig
end

function save_fig(f::Figure; name="tmp", DPI=300)
  path = "./figs/" * name * ".pdf"
  save(path, f, pt_per_unit=1, dpi=DPI)
end


