"""Function to calculate a simple stack averaging solution"""

const STACK_DSET = "stack"
const STACK_FLAT_DSET = "deramped_stack"
const STACK_FLAT_SHIFTED_DSET = "deramped_shifted_stack"


function run_stackavg(unw_stack_file::String, geolist::Array{Date, 1}, valid_igram_list::Array{Tuple{Date, Date}, 1})

    # Figure out which of all the igrams we want to use
    chosen_igrams = pick_igrams(geolist)

    # Get the full list including invalid/ignored to figure out indices to pick
    # from the .h5 file (we dont know this just from valid_igram_list
    full_igram_list = load_intlist_from_h5(unw_stack_file)

    picked_igram_indices = indices_from_full(chosen_igrams, full_igram_list)
    println("Out of $(length(picked_igram_indices)) desired igrams, "*
            "$(sum(isnothing.(picked_igram_indices))) are not available")
    filter!(.!isnothing, picked_igram_indices)

    # Load only these 
    println("Reading $(length(picked_igram_indices)) igrams out of unw stack")


    stack_dset = STACK_DSET
    println("Using '$stack_dset' dataset from file for averaging")
    @time unw_stack = load_hdf5_stack(unw_stack_file, stack_dset, picked_igram_indices)

    # Now with proper igrams picked, just divide the total phase by time diff sum
    timediffs = [(later - early).value for (early, later) in chosen_igrams]
    phase_sum = sum(unw_stack, dims=3)
    avg_velo = phase_sum ./ sum(timediffs)
    return avg_velo
end

"""Given the set of acquisition dates, find which pairs should be used 
as interferograms to average"""
function stack_indices(geolist::Array{Date, 1})
    num_geos = length(geolist)
    # We round up here so that if there is an odd number, 
    # we just ignore the middle date 
    # e.g. dates [1,2,3,4,5] would be [1,4], [2,5]
    idx_span = ceil(num_geos / 2)
    num_igrams = floor(num_geos / 2)
    return [Int.((i, i + idx_span)) for i in 1:num_igrams]
end


function pick_igrams(geolist::Array{Date, 1})::Array{Tuple{Date, Date}, 1}
    geo_pairs = stack_indices(geolist)
    return [(geolist[i1], geolist[i2]) for (i1, i2) in geo_pairs]
end


"""Find the indices of the stack-picked igrams within the full stack
of all compute igrams (that are stored in the unw_stack file)"""
indices_from_full(chosen_igrams, full_igram_list) = indexin(chosen_igrams, full_igram_list)
