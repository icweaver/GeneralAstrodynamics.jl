### A Pluto.jl notebook ###
# v0.20.4

using Markdown
using InteractiveUtils

# ╔═╡ 16301f5f-4ea6-4cbb-8470-710b216e7202
begin
	import Pkg
	Pkg.add(["OrdinaryDiffEq", "LinearAlgebra", "CairoMakie", "Plots", "PlutoLinks"])
	Pkg.activate(Base.current_project())
	Pkg.instantiate()
end

# ╔═╡ 4d2b67f6-61ae-4952-8bc5-03a0042098f8
using PlutoLinks: @revise

# ╔═╡ c224b949-9484-4c19-8bf8-27ecea000aad
@revise using GeneralAstrodynamics

# ╔═╡ 4528d54f-c652-406c-aecc-41866e4eb882
using OrdinaryDiffEq, LinearAlgebra, CairoMakie

# ╔═╡ 43358f87-2d39-4fa9-8ddb-abd121336170
import Plots as Pl

# ╔═╡ 12521b0e-2878-4f45-aabf-3ed8bcfe1872
md"""
## Periodic orbits
"""

# ╔═╡ c491687a-50e3-4db1-80fb-bf4c571e2842
μ = 0.012150584395829193

# ╔═╡ 6c1e0f32-e572-4b4a-9394-53c2e27525c3
sol_planar = let
    ic = halo(μ, 1) # lyapunov (planar) orbit
    u = Vector(CartesianState(ic))
    problem = ODEProblem(CR3BFunction(), u, (0, ic.Δt), (μ,))
    solution = solve(problem, Vern9(), reltol=1e-14, abstol=1e-14)
end

# ╔═╡ 5c755747-b499-4c58-8161-2aeb3e01bd5a
sol_extraplanar = let
    ic = halo(μ, 2; amplitude=0.01) # halo (non-planar) orbit
    u = Vector(CartesianState(ic))
    problem = ODEProblem(CR3BFunction(), u, (0, ic.Δt), (μ,))
    solution = solve(problem, Vern9(), reltol=1e-14, abstol=1e-14)
end

# ╔═╡ 1a35590d-6a3e-4b38-ae02-7844af08b1e8
Pl.plot(
	Pl.plot(sol_planar, idxs=(:x,:y,:z), title = "Lyapunov Orbit", label=:none, size=(1600,900), dpi=400, aspect_ratio=1),
	Pl.plot(sol_extraplanar, idxs=(:x,:y,:z), title = "Halo Orbit", label=:none, size=(1600,900), dpi=400, aspect_ratio=1);
	layout = (1,2),
)

# ╔═╡ 976adbfa-53e6-4063-b0dd-a0efec832cda
let
	fig = Figure(size=(800, 400); fontsize=11)

	common_kwargs = (; aspect=:equal, azimuth=-π/3)

	ax_left = Axis3(fig[1, 1];
		title = "Lyapunov Orbit",
		limits = (0.78, 0.90, -0.09, 0.09, -0.02, 1.04),
		common_kwargs...,
	)
	ax_right = Axis3(fig[1, 2];
		title = "Halo Orbit",
		limits = (1.05, 1.26, -0.1, 0.1, -0.02, 0.01),
		protrusions = (30, 100, 0, 0),
		common_kwargs...,
	)
	
	plot!(ax_left, sol_planar; idxs=(:x, :y, :z))
	plot!(ax_right, sol_extraplanar; idxs=(:x, :y, :z))
	
	fig
end

# ╔═╡ 46abf696-1467-4d16-87d2-f0a76a79ea35
md"""
## Manifold Computations
"""

# ╔═╡ e995801d-f487-44c9-8ff1-cb4351582a35
unstable = let
    ic = halo(μ, 1; amplitude=0.005)

    u = CartesianState(ic)
    Φ = monodromy(u, μ, ic.Δt, CR3BFunction(stm=true))

    ics = let
        problem = ODEProblem(CR3BFunction(stm=true), vcat(u, vec(I(6))), (0, ic.Δt), (μ,))
        solution = solve(problem, Vern9(), reltol=1e-12, abstol=1e-12, saveat=(ic.Δt / 10))

        solution.u
    end

    perturbations = [
        diverge(ic[1:6], reshape(ic[7:end], 6, 6), Φ; eps=-1e-7)
        for ic in ics
    ]

    problem = EnsembleProblem(
        ODEProblem(CR3BFunction(), u, (0.0, 2 * ic.Δt), (μ,)),
        prob_func=(prob, i, repeat) -> remake(prob; u0=perturbations[i]),
    )

    solution = solve(problem, Vern9(), trajectories=length(perturbations), reltol=1e-14, abstol=1e-14)
end

# ╔═╡ 09ae1f89-a573-45ba-9a34-b4542bf607f8
stable = let
    ic = halo(μ, 2; amplitude=0.005)

    u = CartesianState(ic)
    Φ = monodromy(u, μ, ic.Δt, CR3BFunction(stm=true))

    ics = let
        problem = ODEProblem(CR3BFunction(stm=true), vcat(u, vec(I(6))), (0, ic.Δt), (μ,))
        solution = solve(problem, Vern9(), reltol=1e-12, abstol=1e-12, saveat=(ic.Δt / 10))

        solution.u
    end
    
    perturbations = [
        converge(ic[1:6], reshape(ic[7:end], 6, 6), Φ; eps=1e-7)
        for ic in ics
    ]

    problem = EnsembleProblem(
        ODEProblem(CR3BFunction(), u, (0.0, -2.1 * ic.Δt), (μ,)),
        prob_func=(prob, i, repeat) -> remake(prob; u0=perturbations[i]),
    )

    solution = solve(problem, Vern9(), trajectories=length(perturbations), reltol=1e-14, abstol=1e-14)
end

# ╔═╡ 8bee3ddc-4a10-4a1c-971b-1f2364a99891
let
	figure = Pl.plot(; 
	    aspect_ratio = 1.0,
	    background = :transparent,
	    grid = true,
	    title = "Unstable and Stable Invariant Manifolds",
	)
	Pl.plot!(figure, unstable, idxs=(:x, :y), aspect_ratio=1, label=:none, palette=:blues)
	Pl.plot!(figure, stable, idxs=(:x, :y), aspect_ratio=1, label=:none, palette=:blues)
	Pl.scatter!(figure, [1-μ], [0], label="Moon", xlabel="X (Earth-Moon Distance)", ylabel="Y (Earth-Moon Distance)", marker=:x, color=:black, markersize=10,)
end

# ╔═╡ cda6f545-a8cc-4d1c-a50c-75ac23d1a4f4
md"""
!!! warning "Workaround"

	Manually iterating over the colormap for now until point #5 in [SciML/SciMLBase.jl#697 (comment)](https://github.com/SciML/SciMLBase.jl/issues/697#issuecomment-2135801331) is addressed.
"""

# ╔═╡ f9db7e6d-5ba0-4c18-94ba-dc94afd641a0
let
	fig = Figure(size=(800, 400), fontsize=20)
	ax = Axis(fig[1, 1];
		xreversed = true,
		aspect = DataAspect(),
		xlabel = "X (Earth-Moon Distance)",
		ylabel = "Y (Earth-Moon Distance)",
		# labelsize = 12,
		title = "Unstable and Stable Invariant Manifolds",
		titlesize = 24,
	)
	
	cmap = resample_cmap(:blues, length(unstable))
	for (i, m) in enumerate(unstable)
		plot!(ax, m; idxs=(:x, :y), color=cmap[i])
	end

	cmap = resample_cmap(:blues, length(stable))
	for (i, m) in enumerate(stable)
		plot!(ax, m; idxs=(:x, :y), color=cmap[i])
	end

	scatter!(ax, [1-μ], [0]; marker=:x, color=:black, markersize=50, label="Moon")
	
	fig
end

# ╔═╡ Cell order:
# ╠═16301f5f-4ea6-4cbb-8470-710b216e7202
# ╠═4d2b67f6-61ae-4952-8bc5-03a0042098f8
# ╠═c224b949-9484-4c19-8bf8-27ecea000aad
# ╠═4528d54f-c652-406c-aecc-41866e4eb882
# ╠═43358f87-2d39-4fa9-8ddb-abd121336170
# ╟─12521b0e-2878-4f45-aabf-3ed8bcfe1872
# ╠═c491687a-50e3-4db1-80fb-bf4c571e2842
# ╠═6c1e0f32-e572-4b4a-9394-53c2e27525c3
# ╠═5c755747-b499-4c58-8161-2aeb3e01bd5a
# ╠═1a35590d-6a3e-4b38-ae02-7844af08b1e8
# ╠═976adbfa-53e6-4063-b0dd-a0efec832cda
# ╟─46abf696-1467-4d16-87d2-f0a76a79ea35
# ╠═e995801d-f487-44c9-8ff1-cb4351582a35
# ╠═09ae1f89-a573-45ba-9a34-b4542bf607f8
# ╠═8bee3ddc-4a10-4a1c-971b-1f2364a99891
# ╟─cda6f545-a8cc-4d1c-a50c-75ac23d1a4f4
# ╠═f9db7e6d-5ba0-4c18-94ba-dc94afd641a0
