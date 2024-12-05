using JuMP, Gurobi, CSV, DataFrames
using Plots,  Random, Distributions


mutable struct Results
    cost::Float64
    penalty::Float64
    slack_penalty::Float64
    g::Array
    r_up::Array
    r_dn::Array
    γ::Array
    γ_RT::Array
    δ::Array
    δ_RT::Array
    Δ::Array
    creation_time::Float64
    optimization_time::Float64
    termination_status
end

mutable struct Coeficients
    β0_g::Array
    β0_up::Array
    β0_dn::Array
    β_g::Array
    β_up::Array
    β_dn::Array
end

function get_demand_bus(d, d̂, B, d0)
    T, NΩ = size(d)
    new_d = zeros(T[end], NΩ, length(B))
    new_d̂ = zeros(T[end], NΩ, length(B))

    prop = d0 ./ sum(d0)

    for b in 1:length(B)    
        new_d[:, :, b] .= d .* prop[b]
        new_d̂[:, :, b] .= d̂ .* prop[b]
    end
    
    return new_d, new_d̂
    
end

function scale_demand(d_simu, initial_demand)
    
    min_val = minimum(d_simu)
    max_val = maximum(d_simu)
    max_orig = maximum(initial_demand.Initial_demand)
    min_orig = minimum(initial_demand.Initial_demand)

    return (d_simu .- min_val) ./ (max_val - min_val) .* (max_orig - min_orig) .+ min_orig
end

function simulate_demand(df_season, initial_demand, NΩ, B, seed)
    std_error = mean(df_season.val_cargaenergiahomwmed_std)
    T = size(df_season, 1)
    d̂_simu = zeros(T, NΩ)
    Random.seed!(seed)
    for t in 1:T
        d̂_simu[t, :]  = rand(Normal.(df_season[t, 2], df_season[t, 3]),NΩ)
    end

    d_simu = deepcopy(d̂_simu) .+ rand(Normal(0, std_error), T, NΩ)

    d_simu = scale_demand(d_simu, initial_demand)
    d̂_simu = scale_demand(d̂_simu, initial_demand)

    d0    = initial_demand.Initial_demand
    d, d̂  = get_demand_bus(d_simu, d̂_simu, B, d0)
    return d, d̂
end

function get_data(path_data, experiment)

    if experiment == "Dataset - 300 bus system/"
        prefix = "ieee300bus_"
        initial_demand = CSV.read(path_data * experiment * prefix * "initial_demand_data.csv", DataFrame)
    else
        prefix = ""
        initial_demand = nothing
    end
    branch     = CSV.read(path_data * experiment * prefix * "branch_data.csv", DataFrame)
    bus        = CSV.read(path_data * experiment * prefix * "bus_data.csv", DataFrame)
    demand     = CSV.read(path_data * experiment * prefix * "demand_data.csv", DataFrame)
    gen        = CSV.read(path_data * experiment * prefix * "gen_data.csv", DataFrame)
    simulation = CSV.read(path_data * experiment * prefix * "simulation_data.csv", DataFrame)
    return branch, bus, demand, gen, simulation, initial_demand
end

function DispachPerfectInformation(d,d̂)
    # Definir o modelo de otimização
    t1  = time()
    model = Model(Gurobi.Optimizer)
    set_silent(model)
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

    # Função objetivo
    @expression(model, expected_cost, sum(p[ω] *
                                        sum(sum(c_g[i] * g[t, ω, i] + c_up[i] * r_up[t, ω, i] + c_dn[i] * r_dn[t, ω, i] for i in G) +
                                            sum(c_δ[b] * (δ[t, ω, b] + δ_RT[t, ω, b]) + c_γ[b] * (γ[t, ω, b] + γ_RT[t, ω, b]) for b in B)
                                            for t in T) 
                                        for ω in Ω))
    @objective(model, Min, expected_cost)

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

    t2 = time()
    # Resolver o modelo
    optimize!(model)
    t3 = time()

    return Results(value(expected_cost), 0., 0., value.(g), value.(r_up),value.(r_dn),
                    value.(γ),value.(γ_RT),value.(δ), value.(δ_RT), value.(Δ), 
                    t2 - t1, t3 - t2, termination_status(model))
end

function DispachLinear(λ, d, d̂; fixed_β0_g = missing, fixed_β0_up = missing, fixed_β0_dn = missing,
                            fixed_β_g = missing, fixed_β_up = missing, fixed_β_dn = missing)
    NΩ = size(d, 2)
    t1 = time()
    # Definir o modelo de otimização
    model = Model(Gurobi.Optimizer)
    set_silent(model)

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

    @variable(model, s_g_plus[t in T, ω in Ω, i in G] >= 0) 
    @variable(model, s_g_minus[t in T, ω in Ω, i in G] >= 0)
    @variable(model, s_rup_plus[t in T, ω in Ω, i in G] >= 0) 
    @variable(model, s_rup_minus[t in T, ω in Ω, i in G] >= 0)
    @variable(model, s_rdn_plus[t in T, ω in Ω, i in G] >= 0) 
    @variable(model, s_rdn_minus[t in T, ω in Ω, i in G] >= 0)

    if ismissing(fixed_β0_g)
        @variable(model, β0_g[t in T, i in G])
        @variable(model, β0_up[t in T, i in G])
        @variable(model, β0_dn[t in T, i in G])
        @variable(model, β_g[t in T, i in G, b in B, l in L])
        @variable(model, β_up[t in T, i in G, b in B, l in L])
        @variable(model, β_dn[t in T, i in G, b in B, l in L])

    else
        β0_g  = fixed_β0_g
        β0_up =  fixed_β0_up
        β0_dn = fixed_β0_dn

        β_g = fixed_β_g
        β_up = fixed_β_up
        β_dn =  fixed_β_dn
    end

    @variable(model, Φ_g[t in T, i in G, b in B, l in L] >= 0)   # Variável Φ^{(g)}
    @variable(model, Φ_up[t in T, i in G, b in B, l in L] >= 0)  # Variável Φ^{(up)}
    @variable(model, Φ_dn[t in T, i in G, b in B, l in L] >= 0)  # Variável Φ^{(dn)}

    # Função objetivo
    @expression(model, penalty, λ * sum(Φ_g[t, i, b, l] + Φ_up[t, i, b, l] + Φ_dn[t, i, b, l] for l in L, b in B, i in G, t in T))
    @expression(model, slack_penalty, c_slack * (1/NΩ)*sum((s_g_plus[t, ω, i] + s_g_minus[t,ω, i] + s_rup_plus[t,ω, i] + s_rup_minus[t,ω, i] + s_rdn_plus[t,ω, i] + s_rdn_minus[t,ω, i]) for t in T, i in G, ω in Ω))
    @expression(model, expected_cost, sum(p[ω] *
                                        sum(sum(c_g[i] * g[t, ω, i] + c_up[i] * r_up[t, ω, i] + c_dn[i] * r_dn[t, ω, i] for i in G) +
                                            sum(c_δ[b] * (δ[t, ω, b] + δ_RT[t, ω, b]) + c_γ[b] * (γ[t, ω, b] + γ_RT[t, ω, b]) for b in B)
                                            for t in T) 
                                        for ω in Ω))
    @objective(model, Min, expected_cost + penalty + slack_penalty)

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

    @constraint(model, [t in T[nL+1:end], ω in Ω, i in G], g[t, ω, i] == β0_g[t, i] + sum(sum(β_g[t, i, b, l] * d[t-l, ω, b] for b in B) for l in L) + s_g_plus[t, ω, i] - s_g_minus[t, ω, i])  # Restrição (19)
    @constraint(model, [t in T[nL+1:end], ω in Ω, i in G], r_up[t, ω, i] == β0_up[t, i] + sum(sum(β_up[t, i, b, l] * d[t-l, ω, b] for b in B) for l in L) + s_rup_plus[t, ω, i] - s_rup_minus[t, ω, i])  # Restrição (20)
    @constraint(model, [t in T[nL+1:end], ω in Ω, i in G], r_dn[t, ω, i] == β0_dn[t, i] + sum(sum(β_dn[t, i, b, l] * d[t-l, ω, b] for b in B) for l in L) + s_rdn_plus[t, ω, i] - s_rdn_minus[t, ω, i])  # Restrição (21)

    @constraint(model, [t in T, i in G, b in B, l in L], Φ_g[t, i, b, l] - β_g[t, i, b, l] >= 0)  # Restrição (24)
    @constraint(model, [t in T, i in G, b in B, l in L], Φ_g[t, i, b, l] + β_g[t, i, b, l] >= 0)  # RestrAição (25)
    @constraint(model, [t in T, i in G, b in B, l in L], Φ_up[t, i, b, l] - β_up[t, i, b, l] >= 0)  # Restrição (26)
    @constraint(model, [t in T, i in G, b in B, l in L], Φ_up[t, i, b, l] + β_up[t, i, b, l] >= 0)  # Restrição (27)
    @constraint(model, [t in T, i in G, b in B, l in L], Φ_dn[t, i, b, l] - β_dn[t, i, b, l] >= 0)  # Restrição (28)
    @constraint(model, [t in T, i in G, b in B, l in L], Φ_dn[t, i, b, l] + β_dn[t, i, b, l] >= 0)  # Restrição (29)

    # Resolver o modelo
    t2 = time()
    # Resolver o modelo
    optimize!(model)
    t3 = time()
    if ismissing(fixed_β0_dn)
        return  Results(value(expected_cost), value.(penalty), value(slack_penalty), value.(g), value.(r_up),value.(r_dn),
                            value.(γ),value.(γ_RT),value.(δ), value.(δ_RT), value.(Δ), 
                            t2 - t1, t3 - t2, termination_status(model)),
                    Coeficients(value.(β0_g), value.(β0_up), value.(β0_dn), value.(β_g), value.(β_up), value.(β_dn))
    else
        return  Results(value(expected_cost), value.(penalty), value(slack_penalty), value.(g), value.(r_up),value.(r_dn),
                            value.(γ),value.(γ_RT),value.(δ), value.(δ_RT), value.(Δ), 
                            t2 - t1, t3 - t2, termination_status(model)),
                            Coeficients(β0_g, β0_up, β0_dn, β_g, β_up, β_dn)
    end

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

function lift_demand(d,  max_d, min_d, R)

    NΩ = size(d)[2]
    Z = uniform_intervals(R, min_d, max_d)
    #println(Z)
    d_lift = zeros(T[end], NΩ, R);
    for ω in Ω
        d_lift[:, ω, :] = create_lifted_variable(d[:, ω], Z)
    end
    return d_lift
end

function lift_demand_buses(d, d̂, R, B)
    
    d_lift = zeros(size(d)[1], size(d)[2], size(d)[3], R)
    d̂_lift = zeros(size(d)[1], size(d)[2], size(d)[3], R)

    for b in B
        max_d = maximum(d[:, :, b])
        min_d = minimum(d[:, :, b])

        d_lift[:, :, b, :] = lift_demand(d[:, :, b], max_d, min_d, R)
        d̂_lift[:, :, b, :] = lift_demand(d̂[:, :, b], max_d, min_d, R)
    end
    return d_lift, d̂_lift
end

function DispachPiecewiseLinear(λ, d, d̂, d_lift, d̂_lift; fixed_β0_g = missing, fixed_β0_up = missing, fixed_β0_dn = missing,
                                    fixed_β_g = missing, fixed_β_up = missing, fixed_β_dn = missing)
    NΩ = size(d, 2)
    t1 = time()
    # Definir o modelo de otimização
    model = Model(Gurobi.Optimizer)
    set_silent(model)
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

    # @variable(model, β0_g[t in T, i in G])
    # @variable(model, β0_up[t in T, i in G])
    # @variable(model, β0_dn[t in T, i in G])
    # @variable(model, β_g[t in T, i in G, b in B, l in L, r in 1:R])
    # @variable(model, β_up[t in T, i in G, b in B, l in L, r in 1:R])
    # @variable(model, β_dn[t in T, i in G, b in B, l in L, r in 1:R])

    # if !ismissing(fixed_β0_g)
    #     fix.(model[:β0_g], fixed_β0_g; force = true)
    #     fix.(model[:β0_up], fixed_β0_up; force = true)
    #     fix.(model[:β0_dn], fixed_β0_dn; force = true)

    #     fix.(model[:β_g], fixed_β_g; force = true)
    #     fix.(model[:β_up], fixed_β_up; force = true)
    #     fix.(model[:β_dn], fixed_β_dn; force = true)
    # end

    @variable(model, s_g_plus[t in T, ω in Ω, i in G] >= 0) 
    @variable(model, s_g_minus[t in T, ω in Ω, i in G] >= 0)
    @variable(model, s_rup_plus[t in T, ω in Ω, i in G] >= 0) 
    @variable(model, s_rup_minus[t in T, ω in Ω, i in G] >= 0)
    @variable(model, s_rdn_plus[t in T, ω in Ω, i in G] >= 0) 
    @variable(model, s_rdn_minus[t in T, ω in Ω, i in G] >= 0)

    if ismissing(fixed_β0_g)
        @variable(model, β0_g[t in T, i in G])
        @variable(model, β0_up[t in T, i in G])
        @variable(model, β0_dn[t in T, i in G])
        @variable(model, β_g[t in T, i in G, b in B, l in L, r in 1:R])
        @variable(model, β_up[t in T, i in G, b in B, l in L, r in 1:R])
        @variable(model, β_dn[t in T, i in G, b in B, l in L, r in 1:R])
    else
        β0_g  = fixed_β0_g
        β0_up =  fixed_β0_up
        β0_dn = fixed_β0_dn

        β_g = fixed_β_g
        β_up = fixed_β_up
        β_dn =  fixed_β_dn
    end

    @variable(model, Φ_g[t in T, i in G, b in B, l in L, r in 1:R] >= 0)  # Variável Φ^{(g)}
    @variable(model, Φ_up[t in T, i in G, b in B, l in L, r in 1:R] >= 0)  # Variável Φ^{(up)}
    @variable(model, Φ_dn[t in T, i in G, b in B, l in L, r in 1:R] >= 0)  # Variável Φ^{(dn)}

    # Função objetivo
    @expression(model, penalty, λ * sum(Φ_g[t, i, b, l, r] + Φ_up[t, i, b, l, r] + Φ_dn[t, i, b, l, r] for l in L, b in B, i in G, t in T, r in 1:R))
    @expression(model, slack_penalty, c_slack * (1/NΩ)*sum((s_g_plus[t, ω, i] + s_g_minus[t,ω, i] + s_rup_plus[t,ω, i] + s_rup_minus[t,ω, i] + s_rdn_plus[t,ω, i] + s_rdn_minus[t,ω, i]) for t in T, i in G, ω in Ω))
    @expression(model, expected_cost, sum(p[ω] *
                                        sum(sum(c_g[i] * g[t, ω, i] + c_up[i] * r_up[t, ω, i] + c_dn[i] * r_dn[t, ω, i] for i in G) +
                                            sum(c_δ[b] * (δ[t, ω, b] + δ_RT[t, ω, b]) + c_γ[b] * (γ[t, ω, b] + γ_RT[t, ω, b]) for b in B)
                                            for t in T) 
                                        for ω in Ω))
    @objective(model, Min, expected_cost + penalty + slack_penalty)

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

    @constraint(model, [t in T[nL+1:end], ω in Ω, i in G], g[t, ω, i] == β0_g[t, i] + sum(sum(sum(β_g[t, i, b, l, r] * d_lift[t-l, ω, b, r] for b in B) for l in L) for r in 1:R) + s_g_plus[t, ω, i] - s_g_minus[t, ω, i])  # Restrição (19)
    @constraint(model, [t in T[nL+1:end], ω in Ω, i in G], r_up[t, ω, i] == β0_up[t, i] + sum(sum(sum(β_up[t, i, b, l, r] * d_lift[t-l, ω, b, r] for b in B) for l in L) for r in 1:R) + s_rup_plus[t, ω, i] - s_rup_minus[t, ω, i])  # Restrição (20)
    @constraint(model, [t in T[nL+1:end], ω in Ω, i in G], r_dn[t, ω, i] == β0_dn[t, i] + sum(sum(sum(β_dn[t, i, b, l, r] * d_lift[t-l, ω, b, r] for b in B) for l in L) for r in 1:R) + s_rdn_plus[t, ω, i] - s_rdn_minus[t, ω, i])  # Restrição (21)

    @constraint(model, [t in T, i in G, b in B, l in L, r in 1:R], Φ_g[t, i, b, l, r] - β_g[t, i, b, l, r] >= 0)  # Restrição (24)
    @constraint(model, [t in T, i in G, b in B, l in L, r in 1:R], Φ_g[t, i, b, l, r] + β_g[t, i, b, l, r] >= 0)  # Restrição (25)
    @constraint(model, [t in T, i in G, b in B, l in L, r in 1:R], Φ_up[t, i, b, l, r] - β_up[t, i, b, l, r] >= 0)  # Restrição (26)
    @constraint(model, [t in T, i in G, b in B, l in L, r in 1:R], Φ_up[t, i, b, l, r] + β_up[t, i, b, l, r] >= 0)  # Restrição (27)
    @constraint(model, [t in T, i in G, b in B, l in L, r in 1:R], Φ_dn[t, i, b, l, r] - β_dn[t, i, b, l, r] >= 0)  # Restrição (28)
    @constraint(model, [t in T, i in G, b in B, l in L, r in 1:R], Φ_dn[t, i, b, l, r] + β_dn[t, i, b, l, r] >= 0)  # Restrição (29)

    # Resolver o modelo
    t2 = time()
    # Resolver o modelo
    optimize!(model)
    t3 = time()

    if ismissing(fixed_β0_dn)
        return Results(value(expected_cost), value.(penalty), value(slack_penalty), value.(g), value.(r_up),value.(r_dn),
                            value.(γ),value.(γ_RT),value.(δ), value.(δ_RT), value.(Δ), 
                            t2 - t1, t3 - t2, termination_status(model)) ,
                    Coeficients(value.(β0_g), value.(β0_up), value.(β0_dn), value.(β_g), value.(β_up), value.(β_dn))
    else
        return Results(value(expected_cost), value.(penalty), value(slack_penalty), value.(g), value.(r_up),value.(r_dn),
                            value.(γ),value.(γ_RT),value.(δ), value.(δ_RT), value.(Δ), 
                            t2 - t1, t3 - t2, termination_status(model)) ,
                    Coeficients(β0_g, β0_up, β0_dn, β_g, β_up, β_dn)
    end
end

function compute_metrics(results::Results)

    g    = mean(sum(results.g[t, :, i] for t in 1:24, i in 1:length(G)))
    r_up = mean(sum(results.r_up[t, :, i] for t in 1:24, i in 1:length(G)))
    r_dn = mean(sum(results.r_dn[t, :, i] for t in 1:24, i in 1:length(G)))
    γ    = mean(sum(results.γ[t, :, b] for t in 1:24, b in 1:length(B)))
    γ_RT = mean(sum(results.γ_RT[t, :, b] for t in 1:24, b in 1:length(B)))
    δ    = mean(sum(results.δ[t, :, b] for t in 1:24, b in 1:length(B)))
    δ_RT = mean(sum(results.δ_RT[t, :, b] for t in 1:24, b in 1:length(B)))
    Δ    = mean(sum(results.Δ[t, :, i] for t in 1:24, i in 1:length(G)))
    results.termination_status == OPTIMAL ? status = 1 : status = 0
    
    return Dict("g" => g, "r_up" => r_up, "r_dn" => r_dn , "γ" => γ, "γ_RT" => γ_RT, "δ"=> δ, "δ_RT" => δ_RT, "Δ" => Δ,
                "cost" => results.cost, "penalty" => results.penalty, "slack_penalty" => results.slack_penalty, "status" => status)
end

function plot_results(T, y_pi, y_ldr, y_pwl, title, ylabel)
    # theme(:ggplot)
    y_pi_mean  = vec(mean(sum(y_pi[:, :, i] for i in 1:size(y_pi, 3)), dims = 2))
    y_ldr_mean = vec(mean(sum(y_ldr[:, :, i] for i in 1:size(y_ldr, 3)), dims = 2))
    y_pwl_mean = vec(mean(sum(y_pwl[:, :, i] for i in 1:size(y_pwl, 3)), dims = 2))
    plot(1:T, y_pi_mean, label = "PI", linewidth = 2.0)
    plot!(1:T, y_ldr_mean, label = "LDR", linewidth = 2.0)
    plot!(1:T, y_pwl_mean, label = "PWL", linewidth = 2.0)
    title!(title)
    xlabel!("Hours")
    ylabel!(ylabel)
end

path_data = pwd() * "/data/"
path_results = pwd()*"/results/"
experiment = "Dataset - 300 bus system/"

branch, bus, demand, gen, simulation, initial_demand = get_data(path_data, experiment);

df_season = CSV.read(path_data*"typical_day_SE_fall.csv", DataFrame)

# Sets and indexes
Lines   = collect(axes(branch, 1))
B       = collect(1:length(bus.Bus_id))
conv   = Dict(enumerate(bus.Bus_id))
conv_i = Dict([(el, i) for (i, el) in enumerate(bus.Bus_id)])
L_plus  = [findall(x -> x == conv[b], branch.From_Bus) for b in B]
L_minus = [findall(x -> x == conv[b], branch.To_Bus) for b in B]
b_plus  = [[conv_i[j]] for j in branch.From_Bus]
b_minus = [[conv_i[j]] for j in branch.To_Bus]
G       = gen.Gen_id
Ub = [[] for b in 1:maximum(B)]
for b in B
    Ub[b] = findall(x -> x == conv[b], gen.Gen_bus)
end
nL      = simulation.L[1]
L       = collect(1:nL)
T       = collect(1:simulation.T[1])


Load = [i == 0 ? 0 : 1 for i in bus.Load_id]

# Parameters
d0    = initial_demand.Initial_demand
c_g   = gen.Energy_Cost
c_up  = gen.Up_reserve_cost
c_dn  = gen.Dn_reserve_cost
G_max = gen.Max_gen
G_min = gen.Min_gen
c_δ   = ones(length(B))*(2.1*maximum(c_g))
c_slack = maximum(c_δ)
c_γ   = ones(length(B))*(0.1*maximum(c_g))
F     = branch.Transmission_Capacity 
x     = branch.Reactance
RD_dn = gen.Ramp_Limit
RD_up = gen.Ramp_Limit
R_dn  = gen.Max_gen
R_up  = gen.Max_gen
Gb = zeros(maximum(B))
for b in B
    println(b)
    if !isempty(Ub[b])
        Gb[b] = sum(G_max[g] for g in Ub[b]) 
    end
end

R = 5

# Cross Validation for \lambda
λ_values = Int.(collect(400:100:1e3))

NΩ = 5
Ω  = collect(1:NΩ)
p  = ones(NΩ)/NΩ

d, d̂ = simulate_demand(df_season, initial_demand, NΩ, B, 123)
d_lift, d̂_lift = lift_demand_buses(d, d̂, R, B)

# out of sample
NΩ_out = 25
Ω_out  = collect(1:NΩ_out)
p_out  = ones(NΩ_out)/NΩ_out

d_out, d̂_out = simulate_demand(df_season, initial_demand, NΩ_out, B, 456)
d_lift_out, d̂_lift_out = lift_demand_buses(d_out, d̂_out, R, B)

cv_ldr_in = []
cv_pwl_in = []

cv_ldr_out = []
cv_pwl_out = []
# Salvar os resultados para cada lambda
for λ in Int.(λ_values)
    @info(λ)
    @info("In sample...")
    results_ldr_in, coefs_ldr = DispachLinear(λ, d, d̂);
    push!(cv_ldr_in, results_ldr_in.cost);
    metrics_ldr_cv_in = round.(DataFrame(compute_metrics(results_ldr_in)), digits = 3);
    CSV.write(path_results*"results_ldr_insample_cv_lambda_$(λ)_scenarios.csv", metrics_ldr_cv_in);
    results_pwl_in, coefs_pwl = DispachPiecewiseLinear(λ, d, d̂, d_lift, d̂_lift);
    push!(cv_pwl_in, results_pwl_in.cost)
    metrics_pwl_cv_in = round.(DataFrame(compute_metrics(results_pwl_in)), digits = 3)
    CSV.write(path_results*"results_pwl_insample_cv_lambda_$(λ)_scenarios.csv", metrics_pwl_cv_in)

    @info("Out of sample...")
    results_ldr_out, coefs_ldr = DispachLinear(λ, d_out, d̂_out; fixed_β0_g = coefs_ldr.β0_g, 
                                               fixed_β0_up = coefs_ldr.β0_up, fixed_β0_dn = coefs_ldr.β0_dn,
                                               fixed_β_g = coefs_ldr.β_g, fixed_β_up = coefs_ldr.β_up,
                                               fixed_β_dn = coefs_ldr.β_dn);
    push!(cv_ldr_out, results_ldr_out.cost)
    metrics_ldr_cv_out = round.(DataFrame(compute_metrics(results_ldr_out)), digits = 3)
    CSV.write(path_results*"results_ldr_outofsample_cv_lambda_$(λ)_scenarios.csv", metrics_ldr_cv_out)

    results_pwl_out, coefs_pwl = DispachPiecewiseLinear(λ, d_out, d̂_out, d_lift_out, d̂_lift_out; fixed_β0_g = coefs_pwl.β0_g, 
                                                fixed_β0_up = coefs_pwl.β0_up, fixed_β0_dn = coefs_pwl.β0_dn,
                                                fixed_β_g = coefs_pwl.β_g, fixed_β_up = coefs_pwl.β_up,
                                                fixed_β_dn = coefs_pwl.β_dn);
    push!(cv_pwl_out, results_pwl_out.cost);
    metrics_pwl_cv_out = round.(DataFrame(compute_metrics(results_pwl_out)), digits = 3);
    CSV.write(path_results*"results_pwl_outofsample_cv_lambda_$(λ)_scenarios.csv", metrics_pwl_cv_out);
end

λ_ldr = λ_values[findmin(cv_ldr_out)[2]]
λ_pwl = λ_values[findmin(cv_pwl_out)[2]]

# In sample models
NΩ = 40
Ω  = collect(1:NΩ)
p  = ones(NΩ)/NΩ

d, d̂ = simulate_demand(df_season, initial_demand, NΩ, B, 123)
d_lift, d̂_lift = lift_demand_buses(d, d̂, R, B)

results_pi             = DispachPerfectInformation();
results_ldr, coefs_ldr = DispachLinear(λ_ldr, d, d̂);
results_pwl, coefs_pwl = DispachPiecewiseLinear(λ_pwl, d, d̂, d_lift, d̂_lift);


metrics_pi_train  = round.(DataFrame(compute_metrics(results_pi)), digits = 3)
CSV.write(path_results*"results_pi_lambda_$(λ)_scenarios_$(NΩ).csv", metrics_pi_train)

metrics_ldr_train = round.(DataFrame(compute_metrics(results_ldr)), digits = 3)
CSV.write(path_results*"results_ldr_lambda_$(λ)_scenarios_$(NΩ).csv", metrics_ldr_train)

metrics_pwl_train = round.(DataFrame(compute_metrics(results_pwl)), digits = 3)
CSV.write(path_results*"results_pwl_lambda_$(λ)_scenarios_$(NΩ).csv", metrics_pwl_train)

# OUt of sample models
NΩ = simulation.S[1]
d_new, d̂_new = simulate_demand(df_season, initial_demand, NΩ, B, 0)
d_lift_new, d̂_lift_new = lift_demand_buses(d_new, d̂_new, R, B)

results_pi_out = DispachPerfectInformation();
results_ldr_out, coefs_ldr_out = DispachLinear(λ_ldr, d_new, d̂_new; fixed_β0_g = coefs_ldr.β0_g, fixed_β0_up = coefs_ldr.β0_up,
                                            fixed_β0_dn = coefs_ldr.β0_dn, fixed_β_g = coefs_ldr.β_g, 
                                            fixed_β_up = coefs_ldr.β_up, fixed_β_dn = coefs_ldr.β_dn);

results_pwl_out, coefs_pwl_out = DispachPiecewiseLinear(λ_pwl, d_new, d̂_new, d_lift_new, d̂_lift_new; 
                                                        fixed_β0_g = coefs_pwl.β0_g, fixed_β0_up = coefs_pwl.β0_up, 
                                                        fixed_β0_dn = coefs_pwl.β0_dn, fixed_β_g = coefs_pwl.β_g, 
                                                        fixed_β_up = coefs_pwl.β_up, fixed_β_dn = coefs_pwl.β_dn);

metrics_pi_test  = round.(DataFrame(compute_metrics(results_pi_out)), digits = 3)
CSV.write(path_results*"results_pi_out_lambda_$(λ)_scenarios_$(NΩ).csv", metrics_pi_test)

metrics_ldr_test = round.(DataFrame(compute_metrics(results_ldr_out)), digits = 3)
CSV.write(path_results*"results_ldr_out_lambda_$(λ)_scenarios_$(NΩ).csv", metrics_ldr_test)

metrics_pwl_test = round.(DataFrame(compute_metrics(results_pwl_out)), digits = 3)
CSV.write(path_results*"results_pwl_out_lambda_$(λ)_scenarios_$(NΩ).csv", metrics_pwl_test)



plot_results(24, results_pi.g, results_ldr.g, results_pi.g, "title", "label")