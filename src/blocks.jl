struct BlockIterator
    rows::Cint
    cols::Cint
    ni::Cint
    nj::Cint
    n::Cint
    xbsize::Cint
    ybsize::Cint
end

function Base.iterate(obj::BlockIterator, iter::Int = 0)
    iter == obj.n && return nothing
    j = floor(Int, iter / obj.ni)
    i = iter % obj.ni
    nrows = if (i + 1) * obj.ybsize < obj.rows
        obj.ybsize
    else
        obj.rows - i * obj.ybsize
    end
    ncols = if (j + 1) * obj.xbsize < obj.cols
        obj.xbsize
    else
        obj.cols - j * obj.xbsize
    end
    (((i, j), (nrows, ncols)), iter + 1)
end

function blocks(fname::AbstractString, dset_name::AbstractString)
    (chunkrows, chunkcols, _) = get_chunk(fname, dset_name)
    rows, cols, _ = size(fname, dset_name)
    ni = ceil(Cint, rows / chunkrows)
    nj = ceil(Cint, cols / chunkcols)
    BlockIterator(rows, cols, ni, nj, ni * nj, chunkcols, chunkrows)
end

Base.length(b::BlockIterator) = b.n

struct WindowIterator
    blockiter::BlockIterator
end

function Base.iterate(obj::WindowIterator, iter::Int = 0)
    handle = obj.blockiter
    next = Base.iterate(handle, iter)
    next == nothing && return nothing
    (((i, j), (nrows, ncols)), iter) = next
    # (((1:ncols) .+ j * handle.xbsize, (1:nrows) .+ i * handle.ybsize), iter)
    (((1:nrows) .+ i * handle.ybsize, (1:ncols) .+ j * handle.xbsize), iter)
end

windows(fname::AbstractString, dset_name::AbstractString) = WindowIterator(blocks(fname, dset_name))

Base.length(w::WindowIterator) = w.blockiter.n
Base.size(w::WindowIterator) = (length(w),)