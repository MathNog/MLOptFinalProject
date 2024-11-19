using CSV, DataFrames, Plots, Statistics, Dates, Distributions, Random

path_data = pwd()*"/data/"

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

plot(collect(1:24), df_day.demand_mw)

CSV.write(path_data * "typical_day_$(system)_fall.csv", df_fall)
CSV.write(path_data * "typical_day_$(system)_spring.csv", df_spring)
CSV.write(path_data * "typical_day_$(system)_winter.csv", df_winter)
CSV.write(path_data * "typical_day_$(system)_summer.csv", df_summer)