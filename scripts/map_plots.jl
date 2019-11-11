using DataFrames
import Interpolations: interpolate, extrapolate, BSpline, Constant, Flat
import ImageFiltering: imfilter, Kernel

import MapImages: nearest_pixel, nearest, grid, km_to_deg


### Gridding Functions for summing CSVs into bins ###
#
oob(row, col, out) = row < 1 || row > size(out, 1) || col < 1 || col > size(out, 2)

function bin_vals(demrsc, df; valcol=:sum15_17, loncol=:LongNAD27, latcol=:LatNAD27,
                  smoothing::Union{Nothing, <:Real}=nothing)
    out = zeros(size(demrsc))
    for (lon, lat, val) in eachrow(df[:, [loncol, latcol, valcol]])
        row, col = nearest_pixel(demrsc, lat, lon)
        oob(row, col, out) && continue
        out[row, col] += coalesce(val, 0)
    end
    
    return (!isnothing(smoothing) && smoothing > 0) ? imfilter(out, Kernel.gaussian(smoothing)) : out
end

function coarse_bin_vals(demrsc, df; step_km=nothing, digits=nothing, 
                         valcol=:sum15_17, loncol=:LongNAD27, latcol=:LatNAD27,
                         samesize=true, sum_vals=true)
    if !isnothing(step_km)
        lons, lats = coarse_grid_km(demrsc, step_km)
    elseif !isnothing(digits)
        lons, lats = coarse_grid_deg(demrsc, digits)
    else
        error("Need step_km or digits")
    end
    out = zeros(size(lats, 1), size(lons, 1))
    for (lon, lat, val) in eachrow(df[:, [loncol, latcol, valcol]])
        row = nearest(lats, lat)
        col = nearest(lons, lon)
        oob(row, col, out) && continue
        out[row, col] += (sum_vals ? val : 1)
    end
    return samesize ? resize_nearest(out, size(demrsc)...) : out
end

function resize_nearest(img, nrows, ncols)
    intp = interpolate(img, BSpline(Constant()))
    etpf = extrapolate(intp, Flat())
    return etpf(range(.5, stop=size(img, 1)+.5, length=nrows), range(.5, stop=size(img, 2)+.5, length=ncols))
end


"""Create a grid in DemRsc area with spacing `step_km`"""
function coarse_grid_km(demrsc, step_km::Union{AbstractFloat, Nothing}=(10/0.62))
    lons, lats = grid(demrsc, sparse=true)
    step = km_to_deg(step_km)
    return lons[1]:step:lons[end], lats[1]:(-step):lats[end]
end

"""Create a grid from a DemRsc on the 10^(-`digits`) grid lines"""
function coarse_grid_deg(demrsc, digits::Union{Int, Nothing}=2)
    step = 10.0^(-digits)
    lons, lats = grid(demrsc, sparse=true)
    lat1, lat2 = round.(extrema(lats), digits=digits)
    lon1, lon2 = round.(extrema(lons), digits=digits)
    return lon1:step:lon2, lat2:(-step):lat1
end

function plot_loc_per_mile(df, demrsc)
    grid_occurs = coarse_bin_vals(demrsc, df, sum_vals=false, step_km=1/0.62);
end
