using CSV, DataFrames

path_results = pwd()*"/results/"
path_presentation = pwd()*"/presentation"

df_pi_in = CSV.read(path_results*"/results_insample_pi_scenarios_20.csv", DataFrame)
df_ldr_in = CSV.read(path_results*"/results_insample_ldr_lambda_20_scenarios_20.csv", DataFrame)
df_pwl_in = CSV.read(path_results*"/results_insample_pwl_lambda_20_scenarios_20.csv", DataFrame)

df_pi_out = CSV.read(path_results*"/results_pi_out_scenarios_1000.csv", DataFrame)
df_ldr_out = CSV.read(path_results*"/results_ldr_out_lambda_20_scenarios_1000.csv", DataFrame)
df_pwl_out = CSV.read(path_results*"/results_pwl_out_lambda_20_scenarios_1000.csv", DataFrame)

header1 = ["", "In Sample", "In Sample", "In Sample", "Out of Sample", "Out of Sample", "Out of Sample"]
header2 = ["", "PI", "LDR", "PWL", "PI", "LDR", "PWL"]

variables = ["Total Generation", "Redispach", 
            "Reserve Up", "Reserve Down", 
            "Curtailment", "Real Time Curtailment",
            "Load Shedding", "Real Time Load Shedding",
            "Model Creation Time", "Optimization Time"]

g    = vcat(df_pi_in.g, df_ldr_in.g, df_pwl_in.g, df_pi_out.g, df_ldr_out.g, df_pwl_out.g)
Δ    = vcat(df_pi_in.Δ, df_ldr_in.Δ, df_pwl_in.Δ, df_pi_out.Δ, df_ldr_out.Δ, df_pwl_out.Δ)
r_up = vcat(df_pi_in.r_up, df_ldr_in.r_up, df_pwl_in.r_up, df_pi_out.r_up, df_ldr_out.r_up, df_pwl_out.r_up)
r_dn = vcat(df_pi_in.r_dn, df_ldr_in.r_dn, df_pwl_in.r_dn, df_pi_out.r_dn, df_ldr_out.r_dn, df_pwl_out.r_dn)
γ    = vcat(df_pi_in.γ, df_ldr_in.γ, df_pwl_in.γ, df_pi_out.γ, df_ldr_out.γ, df_pwl_out.γ)
γ_RT = vcat(df_pi_in.γ_RT, df_ldr_in.γ_RT, df_pwl_in.γ_RT, df_pi_out.γ_RT, df_ldr_out.γ_RT, df_pwl_out.γ_RT)
δ    = vcat(df_pi_in.δ, df_ldr_in.δ, df_pwl_in.δ, df_pi_out.δ, df_ldr_out.δ, df_pwl_out.δ)
δ_RT = vcat(df_pi_in.δ_RT, df_ldr_in.δ_RT, df_pwl_in.δ_RT, df_pi_out.δ_RT, df_ldr_out.δ_RT, df_pwl_out.δ_RT)
creation_time = vcat(df_pi_in.creation_time, df_ldr_in.creation_time, df_pwl_in.creation_time, df_pi_out.creation_time, df_ldr_out.creation_time, df_pwl_out.creation_time)
optimization_time = vcat(df_pi_in.optimization_time, df_ldr_in.optimization_time, df_pwl_in.optimization_time, df_pi_out.optimization_time, df_ldr_out.optimization_time, df_pwl_out.optimization_time)

values = Matrix(hcat(g , Δ, r_up, r_dn, γ, γ_RT, δ, δ_RT, creation_time, optimization_time)')
results = hcat(variables, values)
results = vcat(reshape(header1, :, 7), reshape(header2, :, 7), results)
results = DataFrame(results,:auto)

CSV.write(path_presentation*"/comparative_results.csv", results)


costs_names = ["Total Generation",
            "Reserve Up", "Reserve Down", 
            "Curtailment", "Real Time Curtailment",
            "Load Shedding", "Real Time Load Shedding"]

total_cost = vcat(df_pi_in.cost, df_ldr_in.cost, df_pwl_in.cost, df_pi_out.cost, df_ldr_out.cost, df_pwl_out.cost)
r_up_cost = vcat(df_pi_in.rup_cost, df_ldr_in.rup_cost, df_pwl_in.rup_cost, df_pi_out.rup_cost, df_ldr_out.rup_cost, df_pwl_out.rup_cost)
r_dn_cost = vcat(df_pi_in.rdn_cost, df_ldr_in.rdn_cost, df_pwl_in.rdn_cost, df_pi_out.rdn_cost, df_ldr_out.rdn_cost, df_pwl_out.rdn_cost)
γ_cost    = vcat(df_pi_in.γ_cost, df_ldr_in.γ_cost, df_pwl_in.γ_cost, df_pi_out.γ_cost, df_ldr_out.γ_cost, df_pwl_out.γ_cost)
γ_RT_cost = vcat(df_pi_in.γ_RT_cost, df_ldr_in.γ_RT_cost, df_pwl_in.γ_RT_cost, df_pi_out.γ_RT_cost, df_ldr_out.γ_RT_cost, df_pwl_out.γ_RT_cost)
δ_cost    = vcat(df_pi_in.δ_cost, df_ldr_in.δ_cost, df_pwl_in.δ_cost, df_pi_out.δ_cost, df_ldr_out.δ_cost, df_pwl_out.δ_cost)
δ_RT_cost = vcat(df_pi_in.δ_RT_cost, df_ldr_in.δ_RT_cost, df_pwl_in.δ_RT_cost, df_pi_out.δ_RT_cost, df_ldr_out.δ_RT_cost, df_pwl_out.δ_RT_cost)


costs_values = Matrix(hcat(total_cost, r_up_cost, r_dn_cost, γ_cost, γ_RT_cost, δ_cost, δ_RT_cost)')
results_costs = hcat(costs_names, costs_values)
results_costs = vcat(reshape(header1, :, 7), reshape(header2, :, 7), results_costs)
results_costs = DataFrame(results_costs,:auto)

CSV.write(path_presentation*"/comparative_costs.csv", results_costs)