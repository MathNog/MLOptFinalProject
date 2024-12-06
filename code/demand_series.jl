using CSV, DataFrames, Plots, Statistics, Dates, Distributions, Random

path_data = pwd()*"/data/"
path_presentation = pwd()*"/presentation/"

function create_datetime_columns!(df, col_name)
    df[!, col_name] = DateTime.(df[!, col_name], "yyyy-mm-dd HH:MM:SS")

    df.day   = day.(df[!, col_name])
    df.month = month.(df[!, col_name])
    df.hour  = Dates.Time.(df[!, col_name]);
end

function filter_subsystem(df, col_name, subsystem)
    return filter(row -> row[col_name] == subsystem, df)
end


function get_season_indexes(df)
    fall_init   = DateTime("2023-03-20 00:00:00", "yyyy-mm-dd HH:MM:SS")
    winter_init = DateTime("2023-06-21 00:00:00", "yyyy-mm-dd HH:MM:SS")
    spring_init = DateTime("2023-09-23 00:00:00", "yyyy-mm-dd HH:MM:SS")
    summer_init = DateTime("2023-12-22 00:00:00", "yyyy-mm-dd HH:MM:SS")

    i_fall   = df[!, "din_instante"] .>= fall_init .&& df[!, "din_instante"] .< winter_init
    i_winter = df[!, "din_instante"] .>= winter_init .&& df[!, "din_instante"] .< spring_init
    i_spring = df[!, "din_instante"] .>= spring_init .&& df[!, "din_instante"] .< summer_init
    i_summer = df[!, "din_instante"] .< fall_init .|| df[!, "din_instante"] .>= summer_init

    return i_fall, i_winter, i_spring, i_summer
end

function get_typical_day(df, col_name)
    return combine(groupby(df, ["hour"]), col_name => mean, col_name => std)
end

function compute_typical_day_demand(df, system)
    df_sist = filter_subsystem(df, "id_subsistema", system)
    create_datetime_columns!(df_sist, "din_instante");
    i_fall, i_winter, i_spring, i_summer = get_season_indexes(df_sist)
    df_fall = get_typical_day(df_sist[i_fall, :], "val_cargaenergiahomwmed")
    df_winter = get_typical_day(df_sist[i_winter, :], "val_cargaenergiahomwmed")
    df_spring = get_typical_day(df_sist[i_spring, :], "val_cargaenergiahomwmed")
    df_summer = get_typical_day(df_sist[i_summer, :], "val_cargaenergiahomwmed") 

    return df_fall, df_winter, df_spring, df_summer
end

df = CSV.read(path_data*"CURVA_CARGA_2023.csv", DataFrame)

system = "SE"

df_fall, df_winter, df_spring, df_summer = compute_typical_day_demand(df, system)


p_mean = plot(df_fall.hour, df_fall.val_cargaenergiahomwmed_mean,
    color = :chocolate, label = "", xlabel = "Hour", ylabel = "Mean Demand (MWh)", title = "South-east hourly mean demand Fall 2023")
p_std = plot(df_fall.hour, df_fall.val_cargaenergiahomwmed_std,
    color = :chocolate, label = "", xlabel = "Hour", ylabel = "Std Demand (MWh)", title = "South-east hourly std demand Fall 2023")
plot(p_mean, p_std, layout = (2,1), size=(800,600))
savefig(path_presentation*"/fall_typical_day.png")


df_sist = filter_subsystem(df, "id_subsistema", system)
create_datetime_columns!(df_sist, "din_instante");

CSV.write(path_data * "typical_day_$(system)_fall.csv", df_fall)
CSV.write(path_data * "typical_day_$(system)_spring.csv", df_spring)
CSV.write(path_data * "typical_day_$(system)_winter.csv", df_winter)
CSV.write(path_data * "typical_day_$(system)_summer.csv", df_summer)

i_fall, i_winter, i_spring, i_summer = get_season_indexes(df_sist)


first_fall   = findfirst(x -> x==1, i_fall)
first_winter = findfirst(x -> x==1, i_winter)
first_spring = findfirst(x -> x==1, i_spring)
first_summer = findlast(x -> x==0, i_summer)

plot(df_sist.din_instante, df_sist.val_cargaenergiahomwmed, label= "")
vline!(df_sist.din_instante[first_fall:first_fall], label = "summer/fall", color = :darkorange)
vline!(df_sist.din_instante[first_winter:first_winter], label = "fall/winter", color = :navyblue)
vline!(df_sist.din_instante[first_spring:first_spring], label = "winter/spring", color = :mediumspringgreen)
vline!(df_sist.din_instante[first_summer:first_summer], label = "spring/summer", color = :orangered)
plot!(title = "South-east Hourly demand 2023",
    xlabel = "Hour", ylabel = "Demand (MWh)",
    size = (850,500))
savefig(path_presentation*"/SE_2023_demand.png")


plot(df_sist.din_instante[1:first_fall], df_sist.val_cargaenergiahomwmed[1:first_fall], label = "summer", color = :orangered)
plot!(df_sist.din_instante[first_fall:first_winter], df_sist.val_cargaenergiahomwmed[first_fall:first_winter],label = "fall", color = :darkorange)
plot!(df_sist.din_instante[first_winter:first_spring], df_sist.val_cargaenergiahomwmed[first_winter:first_spring],label = "winter", color = :navyblue)
plot!(df_sist.din_instante[first_spring:first_summer], df_sist.val_cargaenergiahomwmed[first_spring:first_summer],label = "spring", color = :mediumspringgreen)
plot!(df_sist.din_instante[first_summer:end], df_sist.val_cargaenergiahomwmed[first_summer:end],label = "", color = :orangered)
plot!(title = "South-east Hourly demand 2023",
    xlabel = "Hour", ylabel = "Demand (MWh)",
    size = (850,500))
savefig(path_presentation*"/SE_2023_demand_colors.png")

df_fall_full = df_sist[i_fall, :]
plot(df_fall_full.din_instante, df_fall_full.val_cargaenergiahomwmed, 
        color = :chocolate, label = "", xlabel = "Hour", ylabel = "Demand (MWh)", title = "South-east Hourly demand Fall 2023")
savefig(path_presentation*"/SE_2023_fall_demand.png")
