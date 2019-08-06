using SQLite
using Glob
using LinearAlgebra

struct Ellp
    a::Float64
    e2::Float64
end

const SOL = 299792458
const ELLP = Ellp(6378137.0, 0.0066943799901499996)


function calculate_los(lat, lon, dbfile=nothing)
    if isnothing(dbfile)
        dbfile_list = Glob.glob("*.db*")
        dbfile = dbfile_list[1]
    end
    param_dict = load_all_params(dbfile)

    xyz = _compute_xyz(lat, lon)

    timeorbit, xx, vv = read_orbit_vector(param_dict["orbinfo"])

    for idx=1:param_dict["azimuthBursts"]
        vecr = loop_burst(xyz, param_dict, idx, timeorbit, xx, vv)
        if !isnothing(vecr)
            return vecr
        end
    end
    return nothing
end

function _compute_xyz(lat, lon)
    # TODO: do i wanna get rid of this "params" file?
    dem_file, dem_rsc_file = readlines("params")
    dem_rsc = sario.load(dem_rsc_file)
    dem = sario.load(dem_file)

    firstlon, firstlat = dem_rsc["x_first"], dem_rsc["y_first"]
    deltalon, deltalat = dem_rsc["x_step"], dem_rsc["y_step"]


    col = round(Int, (lon-firstlon)/deltalon)
    row = round(Int, (lat-firstlat)/deltalat)
    
    if row < 1 || col < 1 || row > size(dem, 1) || col > size(dem, 2)
        println("($lat, $lon) out of bounds at ($row, $col) for DEM size $(size(dem))")
        return
    else
        println("($lat, $lon) is at ($row, $col)")
    end

    llh = [ deg2rad(lat), deg2rad(lon), dem[row, col] ]
    xyz = llh_to_xyz(llh)
end


function loop_burst(xyz, param_dict, idx, timeorbit, xx, vv)
    dtaz = 1.0 / param_dict["prf"]  # Nazlooks / prf
    tstart = param_dict["azimuthTimeSeconds$idx"]
    tend  = tstart + (param_dict["linesPerBurst"] - 1) * dtaz
    tmid = (tstart+tend) / 2.0

    println("Burst $idx, Start, stop Acquisition time: $tstart,$tend")
  
    rngstart = param_dict["slantRangeTime"] * SOL / 2.0
    dmrg = SOL / 2.0 / param_dict["rangeSamplingRate"] # Nnum rnglooks * drho
    rngend = rngstart + (param_dict["samplesPerBurst"]-1)*dmrg
    rngmid = 0.50*(rngstart+rngend)

    # latlons = bounds(tstart, tend, rngstart, rngend, timeorbit, xx, vv)

    # if lat > latlons[1] || lat < latlons[2]
    #     print("$lat not within bounds $latlons")
    #     return nothing
    # end


    xyz_mid, vel_mid = intp_orbit(timeorbit, xx, vv, tmid)
    println("Satellite midpoint time,position,velocity: $tmid $xyz_mid $vel_mid")


    satx = xyz_mid
    satv = vel_mid
    tline, rngpix = orbitrangetime(xyz,timeorbit, xx, vv, tmid, xyz_mid, vel_mid)
    if isnothing(tline) || isnothing(rngpix)
        println("Failed on burst $idx")
        return nothing
    end


    satx, satv = intp_orbit(timeorbit, xx, vv, tline)
    dr = xyz - satx
    vecr = dr / rngpix

    println("los vector (away from satellite):, $vecr")
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
    a = ELLP.a
    e2 = ELLP.e2
    lat, lon, h = llh

    rad_earth = a/sqrt(1.0 - e2*sin(lat)^2)

    xyz = similar(llh)
    xyz[1] = (rad_earth + h)*cos(lat)*cos(lon)
    xyz[2] = (rad_earth + h)*cos(lat)*sin(lon)
    xyz[3] = (rad_earth*(1.0 - e2) + h) * sin(lat)
    return xyz
end

function orbitrangetime(xyz, timeorbit, xx, vv, tline0, satx0, satv0)
    # starting state
    tline = tline0
    satx = satx0
    satv = satv0
    
    tprev = tline + 1
    while abs(tline - tprev) > 5e-9
        tprev = tline

        dr = xyz - satx
        rngpix = norm(dr, 2)

        fn = dr' * satv
        fnprime = -norm(satv, 2)^2

        tline = tline - fn / fnprime

        satx, satv = intp_orbit(timeorbit, xx, vv, tline)
        if isnothing(satx) || isnothing(satv)
            return nothing, nothing
        end
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
function  orbithermite(x,v,t,time)
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

