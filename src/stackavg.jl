"""Function to calculate a simple stack averaging solution"""


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
igram_indices(stack_igrams, full_igram_list) = indexin(stack_igrams, full_igram_list)

function run_stackavg(geolist::Array{Date, 1}, intlist::Array{Tuple{Date, Date}, 1})
    return
end
