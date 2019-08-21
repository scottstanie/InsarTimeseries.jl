import Base.view
import Base.lastindex
using Dates
using HDF5

const DATE_FMT = "yyyymmdd"

# dataset names for general 3D stacks
const STACK_DSET = "stack"
const STACK_MEAN_DSET = "stack_mean"
const STACK_FLAT_DSET1 = "stack_deramped_1"
const STACK_FLAT_DSET2 = "stack_deramped_2"
const STACK_FLAT_SHIFTED_DSET1 = "stack_deramped_1_shifted"
const STACK_FLAT_SHIFTED_DSET2 = "stack_deramped_2_shifted"

# Mask file datasets
const GEO_MASK_DSET = "geo"
const IGRAM_MASK_DSET = "igram"
const IGRAM_MASK_SUM_DSET = "igram_sum"

const DEM_RSC_DSET = "dem_rsc"
const GEOLIST_DSET = "geo_dates"
const INTLIST_DSET = "int_dates"

const REFERENCE_ATTR = "reference"
const REFERENCE_STATION_ATTR = "reference_station"

# Type alias for commonly used compositite type
const Igram = Tuple{Date, Date}

import Base.in
# To check if either date of an igram is contained with a geo date list
Base.in(igram::Igram, geo_date_list::AbstractArray{Date}) = ((igram[1] in geo_date_list) || (igram[2] in geo_date_list))

Base.in(date::Date, igram_list::AbstractArray{Igram}) = [date in igram for igram in igram_list]

temporal_baseline(igram::Igram) = (igram[2] - igram[1]).value
temporal_baseline(igram_list::Array{Igram}) = [temporal_baseline(igram) for igram in igram_list]

# Allow the slicing of c[1:end, :]
# TODO: why isn't this defined in HDF5?? prob should do a PR to them
function lastindex(dset::HDF5Dataset, d::Int)
    return size(dset, d)
end

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
