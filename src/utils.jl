import Base.view
import Base.lastindex
using HDF5

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
