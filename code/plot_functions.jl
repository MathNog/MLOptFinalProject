using Plots, CSV, DataFrames

function plot_λ_cv(λ_values, total_cost)
    idx_opt, cost_opt = findmin(total_cost)
    λ_opt = λ_values[idx_opt]

    plot(λ_values, total_cost./1e7, legend = false)
    plot!([λ_opt], [cost_opt], st=:scatter, marker=:star)
    title!("In sample λ cross-validation")
    xlabel!("λ")
    ylabel!("total cost (\$ 10^7)")
end