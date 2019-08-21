using SQLite
using Glob
using LinearAlgebra

const EARTH_SMA = 6378137.0
const EARTH_E2 = 0.0066943799901499996

const SOL = 299792458


"""Calculate the LOS vector in ENU from a lat/lon Array
Can also take in the pre-computed xyz los vectors

Args:
    lat_lons (Array[float]): list of (lat, lon) coordinates to compute LOS vecs
    Optional: xyz_los_vecs (list[tuple[float]]): list of xyz LOS vectors
    Optional: dbfile to read orbit parameters

"""
function los_to_enu(lat_lons::Array{<:AbstractFloat, 2}; xyz_los_vecs=nothing, dbfile=nothing, param_dict=nothing)
    if isnothing(xyz_los_vecs)
        xyz_los_vecs = calculate_los_xyz(lat_lons, dbfile=dbfile, param_dict=param_dict)
    end
    # In projections file:
    return convert_xyz_latlon_to_enu(reshape(lat_lons, 2, :),
                                     reshape(xyz_los_vecs, 3, :))
end

function los_to_enu(lat_lons::Array{<:AbstractFloat, 1}; xyz_los_vecs=nothing, dbfile=nothing, param_dict=nothing)
    return los_to_enu(reshape(lat_lons, :, 1), xyz_los_vecs=xyz_los_vecs, dbfile=dbfile, param_dict=param_dict)
end


filepath(full_path) = joinpath(splitpath(full_path)[1:end-1]...)
#
# Start of fortran converted functions:
#
function calculate_los_xyz(lat::T, lon::T; dbfile::Union{String, Nothing}=nothing, param_dict=nothing) where {T<:AbstractFloat}
    if isnothing(param_dict)
        if isnothing(dbfile)
            dbfile_list = Glob.glob("*.db*")
            dbfile = dbfile_list[1]
        end
        param_dict = load_all_params(dbfile)
    end
    xyz = _compute_xyz(lat, lon)

    orbinfo_filename = param_dict["orbinfo"]  # The .db file doesn't save path
    dbpath = filepath(dbfile)
    orbinfo_file = joinpath(dbpath, orbinfo_filename)

    timeorbit, xx, vv = read_orbit_vector(orbinfo_file)

    for idx=1:param_dict["azimuthBursts"]
        idx = 9
        vecr = compute_burst_vec(xyz, param_dict, idx, timeorbit, xx, vv)
        if !isnothing(vecr)
            return vecr
        end
    end
    return nothing
end


function calculate_los_xyz(lat::T, lon::T, dem, dem_rsc, param_dict, timeorbit, xx, vv) where {T<:AbstractFloat}
    xyz = _compute_xyz(lat, lon, dem, dem_rsc)
    for idx=1:param_dict["azimuthBursts"]
        idx = 9
        vecr = compute_burst_vec(xyz, param_dict, idx, timeorbit, xx, vv)
        if !isnothing(vecr)
            return vecr
        end
    end
    return nothing
end

function calculate_los_xyz(lat_lon::Array{<:AbstractFloat, 1}; dbfile=nothing, param_dict=nothing)
    lat, lon = lat_lon
    return calculate_los_xyz(lat, lon, dbfile=dbfile, param_dict=param_dict)
end

function calculate_los_xyz(lat_lon_vecs::Array{<:AbstractFloat, 2}; dbfile=nothing, param_dict=nothing)
    num_vecs = size(lat_lon_vecs, 2)
    xyz_vecs = Array{eltype(lat_lon_vecs)}(undef, (3, ))
    for i=1:num_vecs
        lat, lon = lat_lon_vecs[:, i]
        xyz_vecs[:, i] .= calculate_los_xyz(lat, lon, dbfile=dbfile, param_dict=param_dict)
    end
    return xyz_vecs
end

function _compute_xyz(lat, lon)
    # TODO: do i wanna get rid of this "params" file?
    # dem_file, dem_rsc_file = readlines("params")

    # TODO: what directory should this be??
    dem_rsc_file = sario.find_rsc_file(directory="../")
    dem_file = replace(dem_rsc_file, ".rsc" => "")
    dem_rsc = load(dem_rsc_file)
    
    row, col = latlon_to_rowcol(dem_rsc, lat, lon)
    # println("($lat, $lon) is at ($row, $col)")

    # Load just the 1 value from the DEM
    dem_height = load(dem_file, (row, col))

    llh = [ deg2rad(lat), deg2rad(lon), dem_height ]
    xyz = llh_to_xyz(llh)
end

function _compute_xyz(lat, lon, dem, dem_rsc)
    row, col = latlon_to_rowcol(dem_rsc, lat, lon)

    # Load just the 1 value from the DEM
    dem_height = dem[row, col]

    llh = [ deg2rad(lat), deg2rad(lon), dem_height ]
    xyz = llh_to_xyz(llh)
end

function latlon_to_rowcol(dem_rsc, lat, lon)
    firstlon, firstlat = dem_rsc["x_first"], dem_rsc["y_first"]
    deltalon, deltalat = dem_rsc["x_step"], dem_rsc["y_step"]

    col = round(Int, (lon-firstlon)/deltalon) + 1
    row = round(Int, (lat-firstlat)/deltalat) + 1
    return row, col
end


function compute_burst_vec(xyz, param_dict, idx, timeorbit, xx, vv)
    dtaz = 1.0 / param_dict["prf"]  # Nazlooks / prf
    tstart = param_dict["azimuthTimeSeconds$idx"]
    tend  = tstart + (param_dict["linesPerBurst"] - 1) * dtaz
    tmid = (tstart + tend) / 2.0

    # println("Burst $idx, Start, stop Acquisition time: $tstart,$tend")
  
    rngstart = param_dict["slantRangeTime"] * SOL / 2.0
    dmrg = SOL / 2.0 / param_dict["rangeSamplingRate"] # Nnum rnglooks * drho
    rngend = rngstart + (param_dict["samplesPerBurst"]-1)*dmrg
    rngmid = 0.50*(rngstart+rngend)

    # TODO: Does it matter if the lat/lon is within this specific burst?
    # If we interpolate to the middle it seems to always be the same
    # 
    # latlons = bounds(tstart, tend, rngstart, rngend, timeorbit, xx, vv)
    # if lat > latlons[1] || lat < latlons[2]
    #     print("$lat not within bounds $latlons")
    #     return nothing
    # end


    xyz_mid, vel_mid = intp_orbit(timeorbit, xx, vv, tmid)

    # println("Satellite midpoint time,position,velocity: $tmid $xyz_mid $vel_mid")

    tline, range = orbitrangetime(xyz,timeorbit, xx, vv, tmid, xyz_mid, vel_mid)
    if isnothing(tline) || isnothing(range)
        println("Failed on burst $idx")
        return nothing
    end


    satx, satv = intp_orbit(timeorbit, xx, vv, tline)
    # TESTING: see how much the East and North change at beginning/end of orbit
    # satx = xx[:, end]
    # satv = vv[:, end]

    # @show satx, satv
    dr = xyz - satx
    vecr = dr / range

    # println("los vector (away from satellite):, $vecr")
    return vecr

end

function load_table_params(dbfile, param_list, param_types=nothing; tablename="file")
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

    # For this, one of each of these per burst up to "azimuthBursts"
    for idx = 1:param_dict["azimuthBursts"]
        burst_params = [
            "azimuthTimeSeconds$idx",
        ]
        types = Dict(
            "azimuthTimeSeconds$idx" => Float64,
        )
        d = load_table_params(dbfile, burst_params, types, tablename=tablename)
        merge!(param_dict, d)
    end

    # Finally, the orbinfo file is a string
    merge!(param_dict, 
           load_table_params(dbfile, ["orbinfo"], tablename=tablename))
    return param_dict
end


function read_orbit_vector(orbtiming_file)
    lines = readlines(orbtiming_file)
    timefirst, timeend, nlines, num_state_vec = map(line -> parse(Int, line), lines[1:4])
    timeorbit = zeros(num_state_vec)
    xx = Array{Float64, 2}(undef, (3, num_state_vec))
    vv = similar(xx)
    # Note: we don't use the accel
    for i = 1:num_state_vec

        floats = map(num -> parse(Float64, num), split(lines[4 + i]))
        timeorbit[i] = floats[1]
        xx[:, i] = floats[2:4]
        vv[:, i] = floats[5:7]
    end
    return timeorbit, xx, vv
end


function llh_to_xyz(llh::Array{<:AbstractFloat, 1})
    lat, lon, h = llh

    rad_earth = EARTH_SMA/sqrt(1.0 - EARTH_E2*sin(lat)^2)

    xyz = similar(llh)
    xyz[1] = (rad_earth + h)*cos(lat)*cos(lon)
    xyz[2] = (rad_earth + h)*cos(lat)*sin(lon)
    xyz[3] = (rad_earth*(1.0 - EARTH_E2) + h) * sin(lat)
    return xyz
end

function orbitrangetime(xyz, timeorbit, xx, vv, tline0, satx0, satv0, tol=1e-6)
    # starting state
    tline = tline0
    satx = satx0
    satv = satv0
    
    idx, max_iter = 1, 100
    tprev = tline + 1  # Need starting guess
    while abs(tline - tprev) > tol && idx < max_iter
        tprev = tline

        dr = xyz - satx
        range = norm(dr, 2)

        fn = dr' * satv
        fnprime = -norm(satv, 2)^2

        tline = tline - fn / fnprime

        satx, satv = intp_orbit(timeorbit, xx, vv, tline)
        if isnothing(satx) || isnothing(satv)
            return nothing, nothing
        end
        idx += 1
    end
    if abs(tline - tprev) > tol
        println("Warning: orbitrangetime didn't converge within $tol (residual $(abs(tline - tprev))")
    end

    dr = xyz - satx
    range = norm(dr, 2)

    return tline, range
end


function intp_orbit(timeorbit, xx, vv, time)
    if isinf(time) || isnan(time)
        return nothing, nothing
    end
    ilocation = (time - timeorbit[1]) / (timeorbit[2]-timeorbit[1])
    # @show time
    # @show ilocation
    ilocation = round(Int, clamp(ilocation, 2, length(timeorbit) - 2))

    xyz_mid, vel_mid = orbithermite(xx[:,ilocation-1:end], vv[:,ilocation-1:end],
                                    timeorbit[ilocation-1:end], time)
    # xyz_mid, vel_mid = orbithermite(xx, vv,timeorbit, time)

    return xyz_mid, vel_mid
end


"""orbithermite - hermite polynomial interpolation of orbits"""
function orbithermite(x,v,t,time)
# inputs
#  x - 3x4 matrix of positions at four times
#  v - 3x4 matrix of velocities
#  t - 4-vector of times for each of the above data points
#  time - time to interpolate orbit to
# 
# outputs
#  xx - position at time time
#  vv - velocity at time time
    n = 4

    h, hdot = zeros(n), zeros(n)
    f0, f1 = zeros(n), zeros(n)
    g0, g1 = zeros(n), zeros(n)
    sum, product = 0, 0


    for i=1:n
        f1[i]=time-t[i]
        sum=0.
        for j=1:n
            if j != i
                sum += 1.0/(t[i]-t[j])
            end
        end
        f0[i]=1.0 - 2.0 * (time-t[i])*sum
    end

    for i=1:n
        product=1.0
        for k=1:n
            if k != i
                product *= (time-t[k])/(t[i]-t[k])
            end
        end
        h[i]=product

        sum=0.0
        for j=1:n
            product=1.0
            for k=1:n
               if k != i && k != j
                   product *= (time-t[k])/(t[i]-t[k])
               end
            end
            if j != i
                sum += 1.0 / (t[i]-t[j])*product
            end
        end
        hdot[i]=sum
    end

    for i=1:n
        g1[i]=h[i] + 2.0 * (time-t[i]) * hdot[i]
        sum=0.
        for j=1:n
            if i != j
                sum += 1.0 / (t[i]-t[j])
            end
        end
        g0[i]=2.0 *(f0[i]*hdot[i]-h[i]*sum)
    end
    
    xout = zeros(3)
    vout = zeros(3)
    for k=1:3
        sum=0.0
        for i=1:n
            sum += (x[k,i]*f0[i]+v[k,i]*f1[i])*h[i]*h[i]
        end
        xout[k]=sum

        sum=0.
        for i=1:n
            sum += (x[k,i]*g0[i]+v[k,i]*g1[i])*h[i]  #*h[i] extra in pdf
        end
        vout[k] = sum
    end

    return xout, vout
end


function _find_db_path(geo_path)
    extra_path = joinpath(geo_path, "extra_files")
    if isdir(extra_path) && length(Glob.glob(extra_path * "/*.db*")) > 0
        return extra_path
    else
        return geo_path
    end
end


function grid(dem_rsc)
    x = range(dem_rsc["x_first"], step=dem_rsc["x_step"], length=dem_rsc["width"])
    y = range(dem_rsc["y_first"], step=dem_rsc["y_step"], length=dem_rsc["file_length"])
    return x, y
end

function los_map(dem_rsc, dbfile, outfile=nothing)
    param_dict = load_all_params(dbfile)
    orbinfo_filename = param_dict["orbinfo"]  # The .db file doesn't save path
    dbpath = filepath(dbfile)
    orbinfo_file = joinpath(dbpath, orbinfo_filename)
    timeorbit, xorbit, vorbit = read_orbit_vector(orbinfo_file)

    dem_rsc_file = sario.find_rsc_file(directory=".")
    @show dem_rsc_file
    dem_rsc = load(dem_rsc_file)
    dem_file = replace(sario.find_rsc_file(directory=".."), ".rsc" => "")
    @show dem_file
    dem = load(dem_file)

    xx, yy = InsarTimeseries.grid(dem_rsc)
    out = Array{Float64, 3}(undef, (length(yy), length(xx), 3))
    Threads.@threads for j in 1:length(xx)
        for i in 1:length(yy)
    # for (j, x) in enumerate(xx)
        # for (i, y) in enumerate(yy)
            # @show i, j, y, x
            y = yy[i]
            x = xx[j]
            xyz_los_vecs = calculate_los_xyz(y, x, dem, dem_rsc, param_dict, timeorbit, xorbit, vorbit)
            # println("$y $x is at ", InsarTimeseries.latlon_to_rowcol(dem_rsc, y, x))
            out[i, j, :] = InsarTimeseries.los_to_enu([y, x], xyz_los_vecs=xyz_los_vecs)
        end
    end
    if !isnothing(outfile)
        println("Writing to $outfile dset 'stack'")
        h5write(outfile, "stack", permutedims(out, (2, 1, 3)))
    end

    return out
end
