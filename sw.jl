using Oceananigans
using Oceananigans.Models: ShallowWaterModel

Nx, Ny = 4 * 16, 4 * 16
Lx, Ly = 4 * 16, 4 * 16

grid = RectilinearGrid(size = (Nx, Ny),
                            x=(0, Lx), y=(0, Ly),
                            topology = (Periodic, Periodic, Flat))

model = ShallowWaterModel(;
    grid=grid,
    gravitational_acceleration=10,
    coriolis=nothing,
    )

x₀, y₀ = Lx/2, Ly/2
H = 10 # unperturbed layer depth
h₀(x, y, z) = H + 0.2 * exp( - ((x - x₀)^2 + (y - y₀)^2)/2 )

set!(model, h=h₀)


uh, vh, h = model.solution
simulation = Simulation(model, Δt=0.01, stop_iteration=200)

using Oceananigans.OutputWriters: NetCDFOutputWriter, IterationInterval

simulation.output_writers[:fields] =
    NetCDFOutputWriter(model, (; uh, vh, h), filename="test.nc", schedule=IterationInterval(1))
        

run!(simulation)

#using JLD2

#file = jldopen(simulation.output_writers[:fields].filepath)
file = NCDataset(simulation.output_writers[:fields].filepath, "r")

iterations = parse.(Int, keys(file["timeseries/t"]))

xh, yh, zh = nodes(h)

using Plots

@info "Making a neat movie of height..."
anim = @animate for (i, iteration) in enumerate(iterations)

    @info "Plotting frame $i from iteration $iteration..."

    t = file["timeseries/t/$iteration"]
    h_snapshot = file["timeseries/h/$iteration"][:, :, 1]

    h_lim = 2.0
    h_levels = range(-h_lim, stop=h_lim, length=20)

    kwargs = (xlabel="x", ylabel="y", aspectratio=1, linewidth=0, colorbar=true,
        xlims=(0, model.grid.Lx), ylims=(0, model.grid.Ly))

    h_plot = contourf(xh, yh, clamp.(h_snapshot', -h_lim, h_lim);
                             color = :balance,
                            levels = h_levels,
                             clims = (-h_lim, h_lim),
                            kwargs...)

    plot(h_plot, title="Height", layout=(1), size=(600, 500))
end
mp4(anim, "Simple Gaussian Height3.mp4", fps = 8) # hide
