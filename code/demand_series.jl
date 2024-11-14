using CSV, DataFrames, Plots, Statistics, Dates

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

function get_typical_day(df, col_name) 
    return combine(groupby(df, ["hour"]), col_name => mean)
end

function compute_typical_day_demand(df, system)
    df_sist = filter_subsystem(df, "id_subsistema", system)
    create_datetime_columns!(df_sist, "din_instante");
    df_day = get_typical_day(df_sist, "val_cargaenergiahomwmed") 
    rename!(df_day, "val_cargaenergiahomwmed_mean" => "demand_mw")
    return df_day
end

df = CSV.read(path_data*"CURVA_CARGA_2023.csv", DataFrame)

system = "SE"

df_day = compute_typical_day_demand(df, system)

plot(df_day.demand_mw)

for system in ["SE", "S", "N", "NE"]
    df_day = compute_typical_day_demand(df, system)
    CSV.write(path_data * "typical_day_$(system).csv", df_day)
end