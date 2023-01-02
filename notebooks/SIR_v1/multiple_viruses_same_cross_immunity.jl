### A Pluto.jl notebook ###
# v0.19.14

using Markdown
using InteractiveUtils

# ╔═╡ 1451a192-60f5-11ed-1626-95e06fb6e8ae
begin
	using Revise
	using Pkg; Pkg.activate("../../../")
	using SIR_partial_sweep
	using Plots
	using StatsBase
end

# ╔═╡ 9e42c547-d21d-4ff0-84c3-0e6bd9d8bad9
begin
	M = 1 # One region
	N = 3 # Many viruses
	md"Number of regions: $M - Number of viruses $N"
end

# ╔═╡ a9621ca5-b982-4ab9-9cae-88c45baf800b
md"Making the cross-immunity matrix"

# ╔═╡ 97f1fb5a-6509-4a51-a87b-0d09732e7025
begin
	b0 = 0.9
	f0 = 0.8
	K = let
		K = rand(N,N)
		for i in 1:N, j in 1:N 
			if i == j
				K[i,i] = 1 
			elseif i > j
				K[i,j] = f0
			elseif j > i
				K[i,j] = b0
			end
		end
		K
	end
	K
end

# ╔═╡ 39d7d19f-09a2-46ed-ac04-662d13b5c007
begin
	X_init = let
		params = SIRParameters(; M, N)
		region = SIRRegion(N; I0 = 1e-6, K)
		SIRState(; regions = [region], parameters = params)
	end
	md"Initial state of the SIR model"
end

# ╔═╡ 86d89520-723f-4f3e-a1f7-931a3e6bd7d2
begin
	sol = let
		tspan = (0, 15/X_init.parameters.γ)
		simulate(X_init, tspan)
	end
	md"Simulating the SIR for time 5/γ"
end

# ╔═╡ 0ae6d869-6823-4b7b-8ef0-86ac77b1d88b
X_eq = equilibrium(X_init); # Expected equilibrium

# ╔═╡ a7bb8647-4696-45c3-8355-2988adc8bed1
let
	g = :I
	
	tvals = range(sol.tspan..., length=100)
	I = [sol[tvals, 1, g, a] for a in 1:N]

	p = plot(ylabel = g)
	for a in 1:N
		plot!(tvals, I[a], label="$g - $a")
	end
	p
end

# ╔═╡ bd2b0fd2-4db8-4189-9b64-93392b6b953e
let 
	tvals = range(sol.tspan..., length=100)
	f = [frequency(sol, tvals, a) for a in 1:N]

	p = plot(ylabel = "frequency")
	for a in 1:N
		plot!(tvals, f[a], label="Freq. $a", color = a)
		hline!([frequency(X_eq, a)], color = a, line = :dash, label="")
	end
	p
end

# ╔═╡ 26d6ec55-87c3-4622-905c-8365a8540931
let 
	tvals = range(sol.tspan..., length=100)
	f = [frequency(sol, tvals, a) for a in 1:N]
	
	p = plot(ylabel = "frequency")

	# Relative freq of 2/1
	i,j = (1,2)
	f_rel = f[j] ./ (f[i] + f[j])
	plot!(tvals, f_rel, label="Freq. $j vs $i", color = 1)
	hline!([(1-f0)/(2-b0-f0)], color = 1, line = :dash, label="")

	# Relative freq of 3/2
	i,j = (2,3)
	f_rel = f[j] ./ (f[i] + f[j])
	plot!(tvals, f_rel, label="Freq. $j vs $i", color = 2)
	hline!([(1-f0)/(2-b0-f0)], color = 2, line = :dash, label="")

	# Relative freq of 3/1
	i,j = (1,3)
	f_rel = f[j] ./ (f[i] + f[j])
	plot!(tvals, f_rel, label="Freq. $j vs $i", color = 3)
	# hline!([(1-f0)/(2-b0-f0)], color = 3, line = :dash, label="")
	p
end

# ╔═╡ f2e08531-a9ca-4864-8fc5-2d00322826bf
let
	i,j = (1,3)
	f = [frequency(sol, [maximum(sol.tspan)], a) for a in 1:N]
	f_rel = f[j] ./ (f[i] + f[j])
end

# ╔═╡ 13da288e-0a44-4cd9-b096-74d1ae5b0d71
((1-f0)^2) / ((1-f0)^2 + (1-b0)^2)

# ╔═╡ c08cb3a3-8c3b-4516-974c-47bb5162d5d4
((1-f0)) / ((1-f0) + (1-b0))

# ╔═╡ Cell order:
# ╠═1451a192-60f5-11ed-1626-95e06fb6e8ae
# ╠═9e42c547-d21d-4ff0-84c3-0e6bd9d8bad9
# ╟─a9621ca5-b982-4ab9-9cae-88c45baf800b
# ╠═97f1fb5a-6509-4a51-a87b-0d09732e7025
# ╠═39d7d19f-09a2-46ed-ac04-662d13b5c007
# ╠═86d89520-723f-4f3e-a1f7-931a3e6bd7d2
# ╠═0ae6d869-6823-4b7b-8ef0-86ac77b1d88b
# ╟─a7bb8647-4696-45c3-8355-2988adc8bed1
# ╟─bd2b0fd2-4db8-4189-9b64-93392b6b953e
# ╠═26d6ec55-87c3-4622-905c-8365a8540931
# ╠═f2e08531-a9ca-4864-8fc5-2d00322826bf
# ╠═13da288e-0a44-4cd9-b096-74d1ae5b0d71
# ╠═c08cb3a3-8c3b-4516-974c-47bb5162d5d4
