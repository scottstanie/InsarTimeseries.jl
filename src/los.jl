using SQLite

struct Ellp
    a::Float64
    e2::Float64
end

const SOL = 299792458
const ELLP = Ellp(6378137.0, 0.0066943799901499996)


function load_table_params(dbfile, param_list; param_types=nothing, tablename="file")
    db = SQLite.DB(dbfile)

    # TODO: How do you really do string interpolation for this library :\
    param_list_str = join(param_list, "\",\"")
    querystr = """SELECT * FROM $tablename WHERE name IN ("$param_list_str") """

    table = SQLite.Query(db, querystr)
    if isnothing(param_types)
        return Dict(row.name => row.value for row in table)
    else
        return Dict(row.name => parse(param_types[row.name], row.value)
                    for row in table)
    end
end

function load_all_params(dbfile; tablename="file")
    param_dict = Dict{String, Any}()
    param_list = [
        "azimuthBursts",
        "linesPerBurst",
        "samplesPerBurst",
        "prf",
        "slantRangeTime",
        "rangeSamplingRate",
    ]
    # To parse into the correct types, since all are strings
    param_types = Dict(
        "azimuthBursts" => Int,
        "linesPerBurst" => Int,
        "samplesPerBurst" => Int,
        "prf" => Float64,
        "rangeSamplingRate" => Float64,
        "slantRangeTime" => Float64,
    )
    d = load_table_params(dbfile, param_list, param_types, tablename=tablename)
    merge!(param_dict, d)

    # One of each of these per burst number up to "azimuthBursts"
    for idx = 1:param_dict["azimuthBursts"]
        burst_params = [
            "azimuthTimeSeconds$idx",
            "firstValidLine$idx",
            "lastValidLine$idx"
        ]
        types = Dict(
            "azimuthTimeSeconds$idx" => Float64,
            "firstValidLine$idx" => Float64,
            "lastValidLine$idx" => Float64,
        )
        d = load_table_params(dbfile, burst_params, types, tablename=tablename)
        merge!(param_dict, d)
    end

    # Finally, the orbinfo file is a string
    merge!(param_dict, 
           load_table_params(dbfile, ["orbinfo"], tablename=tablename))
    return param_dict
end
