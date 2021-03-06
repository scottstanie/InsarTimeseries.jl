import Base: in, view, size, eltype, names

using Dates
using HDF5

# TODO: remove duplication with Sario constants
const MASK_FILENAME = "masks.h5"
const UNW_FILENAME = "unw_stack.h5"
const CC_FILENAME = "cc_stack.h5"

# dataset names for general 3D stacks
const STACK_DSET = "stack"
const STACK_MEAN_DSET = "stack_mean"
const STACK_FLAT_DSET = "stack_flat"
const STACK_FLAT_SHIFTED_DSET = "stack_flat_shifted"

# Mask file datasets
const GEO_MASK_DSET = "geo"
const GEO_MASK_SUM_DSET = "geo_sum"
const IGRAM_MASK_DSET = "igram"

const DEM_RSC_DSET = "dem_rsc"
const GEOLIST_DSET = "geo_dates"
const INTLIST_DSET = "int_dates"

const REFERENCE_ATTR = "reference"
const REFERENCE_STATION_ATTR = "reference_station"

# Type alias for commonly used compositite type
const Igram = Tuple{Date,Date}

# To check if either date of an igram is contained with a geo date list
Base.in(igram::Igram, geo_date_list::AbstractArray{Date}) =
    ((igram[1] in geo_date_list) || (igram[2] in geo_date_list))

Base.in(date::Date, igram_list::AbstractArray{Igram}) =
    [date in igram for igram in igram_list]

temporal_baseline(igram::Igram) = (igram[2] - igram[1]).value
temporal_baseline(igram_list::Array{Igram}) =
    [temporal_baseline(igram) for igram in igram_list]
span(dates::AbstractArray{Date,1}) = (dates[end] - dates[1]).value
span(igrams::AbstractArray{Igram,1}) = (sort(igrams)[end][2] - sort(igrams)[1][1]).value

_get_day_nums(dts::AbstractArray{Date,1}) = [(d - dts[1]).value for d in dts]


# Extend views to work for datasets
function Base.view(dset::HDF5Dataset, I::Vararg{Any})
    idxs = to_indices(dset, I)
    return dset[idxs...]
end

function Base.view(dset::HDF5Dataset, i::Colon, j::Colon, k)
    # Trick to get the slice view of a 3D dset into 2D instead of (m, n, 1)
    return view(dset[i, j, k], :, :)
end

function Base.view(dset::HDF5Dataset, i::Int, j::Int, k::Colon)
    # Hack for 1D depth slice to have 1D instead of (1, 1, N)
    return view(dset[i, j, k], 1, 1, :)
end


function eltype(h5file::String, dset::String)
    h5open(h5file) do f
        return eltype(f[dset])
    end
end

# Allow the HDF5 `names` to work on just filenam
function names(fname::AbstractString)
    h5open(fname) do f
        return names(f)
    end
end
function names(fname::AbstractString, group::AbstractString)
    h5open(fname) do f
        return names(f[group])
    end
end

function get_chunk(fname::AbstractString, dset::AbstractString)
    h5open(fname) do f
        return HDF5.get_chunk(f[dset])
    end
end

# From https://github.com/JuliaLang/julia/blob/master/base/broadcast.jl#L662
# broadcastable(x::Union{AbstractArray,Number,Ref,Tuple,Broadcasted}) = x
# # Default to collecting iterables — which will error for non-iterables
# broadcastable(x) = collect(x)
Base.broadcastable(x::HDF5Dataset) = read(x)
Base.collect(x::HDF5Dataset) = collect(read(x))
Base.iterate(x::HDF5Dataset) = iterate(read(x))
Base.iterate(x::HDF5Dataset, state) = iterate(read(x), state)

"""Get the number of megabytes of RAM available
Note this is higher than "Sys.free_memory()"
"""
function getmemavail()
    out = Pipe()
    lines = read(`free -m`, String)
    memline = split(lines, '\n')[2]
    return parse(Float64, split(memline)[end])
end
