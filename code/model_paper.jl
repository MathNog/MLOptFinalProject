using JuMP, Gurobi, CSV, DataFrames, Random, Distributions

path_data = pwd() * "/data/"
experiment = "Dataset - 3bus system/"
# experiment = "Dataset - 300 bus system"

branch     = CSV.read(path_data * experiment * "branch_data.csv", DataFrame)
bus        = CSV.read(path_data * experiment * "bus_data.csv", DataFrame)
demand     = CSV.read(path_data * experiment * "demand_data.csv", DataFrame)
gen        = CSV.read(path_data * experiment * "gen_data.csv", DataFrame)
simulation = CSV.read(path_data * experiment * "simulation_data.csv", DataFrame)

# Sets and indexes
Lines   = collect(axes(branch, 1))
B       = bus.Bus_id
L_plus  = [findall(x -> x == b, branch.From_Bus) for b in B]
L_minus = [findall(x -> x == b, branch.To_Bus) for b in B]
b_plus  = [[j] for j in branch.From_Bus]
b_minus = [[j] for j in branch.To_Bus]
G       = bus.Gen_id
Ub      = [findall(x -> x == b, gen.Gen_Bus) for b in B]
nL      = simulation.L[1]
L       = collect(1:nL)
T       = collect(1:simulation.T[1])
NΩ      = simulation.S[1]
Ω       = collect(1:NΩ)
Load    = bus.Load_id

# Parameters
d0    = simulation.d0[1]
c_g   = gen.Energy_Cost
c_up  = gen.Up_reserve_cost
c_dn  = gen.Dn_reserve_cost
G_max = gen.Max_gen
G_min = gen.Min_gen
c_δ   = ones(length(B))*(1.1*maximum(c_g))
c_γ   = ones(length(B))*(0.1*maximum(c_g))
F     = branch.Transmission_Capacity 
x     = branch.Reactance
p     = ones(NΩ)/NΩ
RD_dn = gen.Ramp_Limit
RD_up = gen.Ramp_Limit
R_dn  = gen.Max_gen
R_up  = gen.Max_gen
Gb    = [sum(G_max[g] for g in Ub[b]) for b in B]

d = zeros(T[end]+1, NΩ, B[end])
d̂ = zeros(T[end]+1, NΩ, B[end])
d[1, :, :] .= d0
d̂[1, :, :] .= d0

Random.seed!(1234)
ϵ = rand(Normal(0, 3), (T[end]+1, NΩ, B[end]))
ϕ = simulation.phi_1
for t in T
    d[t+1, :, :] .= (ϕ.*d[t,:,:] .+ ϵ[t+1, :, :]) .* Load'
    d̂[t+1,:,:]   .= (ϕ.*d[t,:,:]) .* Load'
end

d = d[2:end,:,:]
d̂ = d̂[2:end,:,:]

plot(d̂[:, 1, 3])
plot(d[:, 1, 3])

function DispachLinear(λ)
    # Definir o modelo de otimização
    model = Model(Gurobi.Optimizer)

    # Índices e conjuntos
    @variable(model, g[t in T, ω in Ω, i in G] >= 0)  # Variável g para geração
    @variable(model, r_up[t in T, ω in Ω, i in G] >= 0)  # Variável de reserva para r^up
    @variable(model, r_dn[t in T, ω in Ω, i in G] >= 0)  # Variável de reserva para r^dn

    @variable(model, f_RT[t in T, ω in Ω, j in Lines])  # Variável f^{RT} para fluxo
    @variable(model, f[t in T, ω in Ω, j in Lines])  # Variável f para fluxo

    @variable(model, γ[t in T, ω in Ω, b in B] >= 0)  # Variável γ para cargas
    @variable(model, γ_RT[t in T, ω in Ω, b in B] >= 0)  # Variável γ^{RT}

    @variable(model, δ[t in T, ω in Ω, b in B] >= 0)  # Variável δ para derivação
    @variable(model, δ_RT[t in T, ω in Ω, b in B] >= 0)  # Variável δ^{RT}

    @variable(model, θ[t in T, ω in Ω, b in B] >= 0)  # Variável θ para derivação
    @variable(model, θ_RT[t in T, ω in Ω, b in B] >= 0)  # Variável θ^{RT}

    @variable(model, Δ[t in T, ω in Ω, i in G])

    @variable(model, β0_g[t in T, i in G])
    @variable(model, β0_up[t in T, i in G])
    @variable(model, β0_dn[t in T, i in G])
    @variable(model, β_g[t in T, i in G, b in B, l in L])
    @variable(model, β_up[t in T, i in G, b in B, l in L])
    @variable(model, β_dn[t in T, i in G, b in B, l in L])

    @variable(model, Φ_g[t in T, i in G, b in B, l in L] >= 0)  # Variável Φ^{(g)}
    @variable(model, Φ_up[t in T, i in G, b in B, l in L] >= 0)  # Variável Φ^{(up)}
    @variable(model, Φ_dn[t in T, i in G, b in B, l in L] >= 0)  # Variável Φ^{(dn)}

    # Função objetivo
    @expression(model, expected_cost, sum(p[ω] *
                                        sum(sum(c_g[i] * g[t, ω, i] + c_up[i] * r_up[t, ω, i] + c_dn[i] * r_dn[t, ω, i] for i in G) +
                                            sum(c_δ[b] * (δ[t, ω, b] + δ_RT[t, ω, b]) + c_γ[b] * (γ[t, ω, b] + γ_RT[t, ω, b]) for b in B)
                                            for t in T) 
                                        for ω in Ω))
    @objective(model, Min, expected_cost + 
                            λ * sum(Φ_g[t, i, b, l] + Φ_up[t, i, b, l] + Φ_dn[t, i, b, l] for l in L, b in B, i in G, t in T))

    # Restrição (2)
    @constraint(model, Rest2[ω in Ω, t in T, b in B], sum(g[t, ω, i] for i in Ub[b]) + sum(f[t, ω, j] for j in L_plus[b]) - sum(f[t, ω, j] for j in L_minus[b]) +
                                                    δ[t, ω, b] - γ[t, ω, b] == d̂[t, ω, b])

    # Restrição (3)
    @constraint(model, Rest3[ω in Ω, t in T, b in B], sum(g[t, ω, i] +  Δ[t, ω, i] for i in Ub[b]) + 
                                                    sum(f_RT[t, ω, j] for j in L_plus[b]) - sum(f_RT[t, ω, j] for j in L_minus[b]) + 
                                                    δ_RT[t, ω, b] - γ_RT[t, ω, b] == d[t, ω, b])

    # Restrição (4)
    @constraint(model, [ω in Ω, t in T, j in Lines], f[t, ω, j] == (θ[t, ω, b_plus[j]...] - θ[t, ω, b_minus[j]...]) / x[j])

    # Restrição (5)
    @constraint(model, [ω in Ω, t in T, j in Lines], f_RT[t, ω, j] == (θ_RT[t, ω, b_plus[j]...] - θ_RT[t, ω, b_minus[j]...]) / x[j])

    # Restrição (6)
    @constraint(model, [t in T, ω in Ω, j in Lines], -F[j] <= f[t, ω, j] <= F[j])

    # Restrição (7)
    @constraint(model, [t in T, ω in Ω, j in Lines], -F[j] <= f_RT[t, ω, j] <= F[j])

    @constraint(model, [t in T, ω in Ω, i in G], G_min[i] <= g[t, ω, i] <= G_max[i])  # Restrição (8)
    @constraint(model, [t in T, ω in Ω, i in G], g[t, ω, i] + r_up[t, ω, i] <= G_max[i])  # Restrição (9)
    @constraint(model, [t in T, ω in Ω, i in G], g[t, ω, i] - r_dn[t, ω, i] >= G_min[i])  # Restrição (10)

    @constraint(model, [t in T[2:end], ω in Ω, i in G], -RD_dn[i] <= g[t, ω, i] + Δ[t, ω, i] - g[t-1, ω, i] - Δ[t-1, ω, i] <= RD_up[i])  # Restrição (11)
    @constraint(model, [t in T, ω in Ω, i in G], -r_dn[t, ω, i] <= Δ[t, ω, i])  # Restrição (12a)
    @constraint(model, [t in T, ω in Ω, i in G],  r_up[t, ω, i] >= Δ[t, ω, i])  # Restrição (12b)
    @constraint(model, [t in T, ω in Ω, i in G], r_up[t, ω, i] <= R_up[i])  # Restrição (13)
    @constraint(model, [t in T, ω in Ω, i in G], r_dn[t, ω, i] <= R_dn[i])  # Restrição (14)

    @constraint(model, [t in T, ω in Ω, b in B], δ[t, ω, b] <= d̂[t, ω, b])  # Restrição (15)
    @constraint(model, [t in T, ω in Ω, b in B], δ_RT[t, ω, b] <= d[t, ω, b])  # Restrição (16)
    @constraint(model, [t in T, ω in Ω, b in B], γ[t, ω, b] <= Gb[b])  # Restrição (17)
    @constraint(model, [t in T, ω in Ω, b in B], γ_RT[t, ω, b] <= Gb[b])  # Restrição (18)

    @constraint(model, [t in T[nL+1:end], ω in Ω, i in G], g[t, ω, i] == β0_g[t, i] + sum(sum(β_g[t, i, b, l] * d[t-l, ω, b] for b in B) for l in L))  # Restrição (19)
    @constraint(model, [t in T[nL+1:end], ω in Ω, i in G], r_up[t, ω, i] == β0_up[t, i] + sum(sum(β_up[t, i, b, l] * d[t-l, ω, b] for b in B) for l in L))  # Restrição (20)
    @constraint(model, [t in T[nL+1:end], ω in Ω, i in G], r_dn[t, ω, i] == β0_dn[t, i] + sum(sum(β_dn[t, i, b, l] * d[t-l, ω, b] for b in B) for l in L))  # Restrição (21)

    @constraint(model, [t in T, i in G, b in B, l in L], Φ_g[t, i, b, l] - β_g[t, i, b, l] >= 0)  # Restrição (24)
    @constraint(model, [t in T, i in G, b in B, l in L], Φ_g[t, i, b, l] + β_g[t, i, b, l] >= 0)  # Restrição (25)
    @constraint(model, [t in T, i in G, b in B, l in L], Φ_up[t, i, b, l] - β_up[t, i, b, l] >= 0)  # Restrição (26)
    @constraint(model, [t in T, i in G, b in B, l in L], Φ_up[t, i, b, l] + β_up[t, i, b, l] >= 0)  # Restrição (27)
    @constraint(model, [t in T, i in G, b in B, l in L], Φ_dn[t, i, b, l] - β_dn[t, i, b, l] >= 0)  # Restrição (28)
    @constraint(model, [t in T, i in G, b in B, l in L], Φ_dn[t, i, b, l] + β_dn[t, i, b, l] >= 0)  # Restrição (29)

    # Resolver o modelo
    optimize!(model)

    return value.(g), value(expected_cost)
end

function create_lifted_variable(x::Vector{Float64}, Z::Vector{Float64})

    N = length(x)
    R = length(Z)

    Z = sort(Z)
    ξr = zeros(N, R)

    for r in 1:R
        if r == 1
            ξr[:, r] = min.(x, Z[r])
        elseif r == R
            ξr[:, r] = max.(x .- Z[r-1], 0)
        else
            ξr[:, r] = max.(min.(x, Z[r]) .- Z[r-1], 0 )
        end
    end
    return ξr
end

function uniform_intervals(R, A, B)
    # Verifica se A < B, caso contrário, inverte os valores
    A, B = min(A, B), max(A, B)
    
    # Calcula o tamanho de cada intervalo
    step = (B - A) / R
    
    # Gera os limites superiores dos intervalos
    upper_limits = [A + (i+1)*step for i in 0:R-1]
    
    return upper_limits
end

function DispachPiecewiseLinear(λ)
    # Definir o modelo de otimização
    model = Model(Gurobi.Optimizer)

    # Índices e conjuntos
    @variable(model, g[t in T, ω in Ω, i in G] >= 0)  # Variável g para geração
    @variable(model, r_up[t in T, ω in Ω, i in G] >= 0)  # Variável de reserva para r^up
    @variable(model, r_dn[t in T, ω in Ω, i in G] >= 0)  # Variável de reserva para r^dn

    @variable(model, f_RT[t in T, ω in Ω, j in Lines])  # Variável f^{RT} para fluxo
    @variable(model, f[t in T, ω in Ω, j in Lines])  # Variável f para fluxo

    @variable(model, γ[t in T, ω in Ω, b in B] >= 0)  # Variável γ para cargas
    @variable(model, γ_RT[t in T, ω in Ω, b in B] >= 0)  # Variável γ^{RT}

    @variable(model, δ[t in T, ω in Ω, b in B] >= 0)  # Variável δ para derivação
    @variable(model, δ_RT[t in T, ω in Ω, b in B] >= 0)  # Variável δ^{RT}

    @variable(model, θ[t in T, ω in Ω, b in B] >= 0)  # Variável θ para derivação
    @variable(model, θ_RT[t in T, ω in Ω, b in B] >= 0)  # Variável θ^{RT}

    @variable(model, Δ[t in T, ω in Ω, i in G])

    @variable(model, β0_g[t in T, i in G])
    @variable(model, β0_up[t in T, i in G])
    @variable(model, β0_dn[t in T, i in G])
    @variable(model, β_g[t in T, i in G, b in B, l in L, r in 1:R])
    @variable(model, β_up[t in T, i in G, b in B, l in L, r in 1:R])
    @variable(model, β_dn[t in T, i in G, b in B, l in L, r in 1:R])

    @variable(model, Φ_g[t in T, i in G, b in B, l in L, r in 1:R] >= 0)  # Variável Φ^{(g)}
    @variable(model, Φ_up[t in T, i in G, b in B, l in L, r in 1:R] >= 0)  # Variável Φ^{(up)}
    @variable(model, Φ_dn[t in T, i in G, b in B, l in L, r in 1:R] >= 0)  # Variável Φ^{(dn)}

    # Função objetivo
    @expression(model, expected_cost, sum(p[ω] *
                                        sum(sum(c_g[i] * g[t, ω, i] + c_up[i] * r_up[t, ω, i] + c_dn[i] * r_dn[t, ω, i] for i in G) +
                                            sum(c_δ[b] * (δ[t, ω, b] + δ_RT[t, ω, b]) + c_γ[b] * (γ[t, ω, b] + γ_RT[t, ω, b]) for b in B)
                                            for t in T) 
                                        for ω in Ω))
    @objective(model, Min, expected_cost + 
                            λ * sum(Φ_g[t, i, b, l, r] + Φ_up[t, i, b, l, r] + Φ_dn[t, i, b, l, r] for l in L, b in B, i in G, t in T, r in 1:R))

    # Restrição (2)
    @constraint(model, Rest2[ω in Ω, t in T, b in B], sum(g[t, ω, i] for i in Ub[b]) + sum(f[t, ω, j] for j in L_plus[b]) - sum(f[t, ω, j] for j in L_minus[b]) +
                                                    δ[t, ω, b] - γ[t, ω, b] == d̂[t, ω, b])

    # Restrição (3)
    @constraint(model, Rest3[ω in Ω, t in T, b in B], sum(g[t, ω, i] +  Δ[t, ω, i] for i in Ub[b]) + 
                                                    sum(f_RT[t, ω, j] for j in L_plus[b]) - sum(f_RT[t, ω, j] for j in L_minus[b]) + 
                                                    δ_RT[t, ω, b] - γ_RT[t, ω, b] == d[t, ω, b])

    # Restrição (4)
    @constraint(model, [ω in Ω, t in T, j in Lines], f[t, ω, j] == (θ[t, ω, b_plus[j]...] - θ[t, ω, b_minus[j]...]) / x[j])

    # Restrição (5)
    @constraint(model, [ω in Ω, t in T, j in Lines], f_RT[t, ω, j] == (θ_RT[t, ω, b_plus[j]...] - θ_RT[t, ω, b_minus[j]...]) / x[j])

    # Restrição (6)
    @constraint(model, [t in T, ω in Ω, j in Lines], -F[j] <= f[t, ω, j] <= F[j])

    # Restrição (7)
    @constraint(model, [t in T, ω in Ω, j in Lines], -F[j] <= f_RT[t, ω, j] <= F[j])

    @constraint(model, [t in T, ω in Ω, i in G], G_min[i] <= g[t, ω, i] <= G_max[i])  # Restrição (8)
    @constraint(model, [t in T, ω in Ω, i in G], g[t, ω, i] + r_up[t, ω, i] <= G_max[i])  # Restrição (9)
    @constraint(model, [t in T, ω in Ω, i in G], g[t, ω, i] - r_dn[t, ω, i] >= G_min[i])  # Restrição (10)

    @constraint(model, [t in T[2:end], ω in Ω, i in G], -RD_dn[i] <= g[t, ω, i] + Δ[t, ω, i] - g[t-1, ω, i] - Δ[t-1, ω, i] <= RD_up[i])  # Restrição (11)
    @constraint(model, [t in T, ω in Ω, i in G], -r_dn[t, ω, i] <= Δ[t, ω, i])  # Restrição (12a)
    @constraint(model, [t in T, ω in Ω, i in G],  r_up[t, ω, i] >= Δ[t, ω, i])  # Restrição (12b)
    @constraint(model, [t in T, ω in Ω, i in G], r_up[t, ω, i] <= R_up[i])  # Restrição (13)
    @constraint(model, [t in T, ω in Ω, i in G], r_dn[t, ω, i] <= R_dn[i])  # Restrição (14)

    @constraint(model, [t in T, ω in Ω, b in B], δ[t, ω, b] <= d̂[t, ω, b])  # Restrição (15)
    @constraint(model, [t in T, ω in Ω, b in B], δ_RT[t, ω, b] <= d[t, ω, b])  # Restrição (16)
    @constraint(model, [t in T, ω in Ω, b in B], γ[t, ω, b] <= Gb[b])  # Restrição (17)
    @constraint(model, [t in T, ω in Ω, b in B], γ_RT[t, ω, b] <= Gb[b])  # Restrição (18)

    @constraint(model, [t in T[nL+1:end], ω in Ω, i in G], g[t, ω, i] == β0_g[t, i] + sum(sum(sum(β_g[t, i, b, l, r] * d_lift[t-l, ω, b, r] for b in B) for l in L) for r in 1:R))  # Restrição (19)
    @constraint(model, [t in T[nL+1:end], ω in Ω, i in G], r_up[t, ω, i] == β0_up[t, i] + sum(sum(sum(β_up[t, i, b, l, r] * d_lift[t-l, ω, b, r] for b in B) for l in L) for r in 1:R))  # Restrição (20)
    @constraint(model, [t in T[nL+1:end], ω in Ω, i in G], r_dn[t, ω, i] == β0_dn[t, i] + sum(sum(sum(β_dn[t, i, b, l, r] * d_lift[t-l, ω, b, r] for b in B) for l in L) for r in 1:R))  # Restrição (21)

    @constraint(model, [t in T, i in G, b in B, l in L, r in 1:R], Φ_g[t, i, b, l, r] - β_g[t, i, b, l, r] >= 0)  # Restrição (24)
    @constraint(model, [t in T, i in G, b in B, l in L, r in 1:R], Φ_g[t, i, b, l, r] + β_g[t, i, b, l, r] >= 0)  # Restrição (25)
    @constraint(model, [t in T, i in G, b in B, l in L, r in 1:R], Φ_up[t, i, b, l, r] - β_up[t, i, b, l, r] >= 0)  # Restrição (26)
    @constraint(model, [t in T, i in G, b in B, l in L, r in 1:R], Φ_up[t, i, b, l, r] + β_up[t, i, b, l, r] >= 0)  # Restrição (27)
    @constraint(model, [t in T, i in G, b in B, l in L, r in 1:R], Φ_dn[t, i, b, l, r] - β_dn[t, i, b, l, r] >= 0)  # Restrição (28)
    @constraint(model, [t in T, i in G, b in B, l in L, r in 1:R], Φ_dn[t, i, b, l, r] + β_dn[t, i, b, l, r] >= 0)  # Restrição (29)

    # Resolver o modelo
    optimize!(model)

    return value.(g), value(expected_cost)
end


R     = 5
min_d = minimum(d) 
max_d = maximum(d)
Z     = uniform_intervals(R, min_d, max_d)

d_lift = zeros(T[end], NΩ, B[end], R);
d̂_lift = zeros(T[end], NΩ, B[end], R);
for b in B, ω in Ω
    d_lift[:, ω, b, :] = create_lifted_variable(d[:, ω, b], Z)
    d̂_lift[:, ω, b, :] = create_lifted_variable(d̂[:, ω, b], Z)
end

g_opt, cost_opt = DispachLinear(0);
g_opt_pwl, cost_opt_pwl = DispachPiecewiseLinear(0);