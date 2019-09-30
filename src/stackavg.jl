"""Function to calculate a simple stack averaging solution"""

"""Averages all igrams in a stack

saves a 2D array of velocities: mm per year
"""
function run_stackavg(unw_stack_file::String, input_dset::String, outfile::String, outdset::String,
                      geolist::Array{Date, 1}, valid_igram_list::Array{Igram, 1}; stack_all=true)

    # Figure out which of all the igrams we want to use
    # chosen_igrams = pick_igrams(geolist)
    if stack_all
        chosen_igrams = valid_igram_list
    else
        chosen_igrams = pick_igrams(geolist, valid_igram_list)
    end

    # Get the full list including invalid/ignored to figure out indices to pick
    # from the .h5 file (we dont know this just from valid_igram_list, some are ignored)
    full_igram_list = load_intlist_from_h5(unw_stack_file)

    picked_igram_indices = indices_from_full(chosen_igrams, full_igram_list)
    println("Out of $(length(picked_igram_indices)) desired igrams, "*
            "$(sum(isnothing.(picked_igram_indices))) are missing")
    filter!(.!isnothing, picked_igram_indices)

    # Load only these 
    # input_dset = STACK_DSET
    # input_dset = STACK_FLAT_DSET
    println("Reading $(length(picked_igram_indices)) igrams out of '$input_dset' dataset from file ")
    @time unw_stack = load_hdf5_stack(unw_stack_file, input_dset, picked_igram_indices)

    # Now with proper igrams picked, just divide the total phase by time diff sum
    timediffs = temporal_baseline(chosen_igrams)
    phase_sum = sum(unw_stack, dims=3)[:, :, 1]
    avg_velo = phase_sum ./ sum(timediffs)

    # Finally, save as a mm/year velocity
    out = Float32.(P2MM * avg_velo)

    println("Writing solution into $outfile : $outdset")
    h5open(outfile, "cw") do f
        f[outdset] = permutedims(out)
        # TODO: Do I care to add this for stack when it's all the same?
        # f[_count_dset(cur_outdset)] = countstack
    end

    return outfile, outdset
end


"""Given the set of acquisition dates, find which pairs should be used 
as interferograms to average

"""
function stack_indices(geolist::Array{Date, 1})
    num_geos = length(geolist)
    # We round up here so that if there is an odd number, 
    # we just ignore the middle date 
    # e.g. dates [1,2,3,4,5] would be [1,4], [2,5]
    idx_span = ceil(num_geos / 2)
    num_igrams = floor(num_geos / 2)
    return [Int.((i, i + idx_span)) for i in 1:num_igrams]
end



"""Remove the first instance of `value` from `array`
If it is not containing, will return array
"""
function removefirst!(array, value)
    idx = findfirst(isequal(value), array)
    return isnothing(idx) ? array : deleteat!(array, idx)
end

function remove_all_igrams!(array, geo_dates) 
    for date in geo_dates
        for idx = 1:2
            deleteat!(array, findall(x -> x[idx] == date, array))
        end
    end
end


"""Get all igrams that start with `geo_date`"""
find_igrams(igram_list::Array{Igram, 1}, geo_date::Date) = filter(ig -> ig[1] == geo_date, igram_list)

max_igram(igram_list, geo_date) = find_igrams(igram_list, geo_date)[end]

function pick_igrams(full_geolist::Array{Date, 1}, full_igram_list::Array{Igram, 1})
    geolist = copy(full_geolist)
    igram_list = copy(full_igram_list)

    total_timespan = (geolist[end] - geolist[1]).value  
    # Find this value to make sure no igram goes over half
    half_span = total_timespan / 2

    chosen_igrams = Array{Igram, 1}()

    while length(geolist) > 1
        cur_geo = popfirst!(geolist)
        longest_igram = max_igram(igram_list, cur_geo)
        push!(chosen_igrams, longest_igram)

        early, late = longest_igram
        remove_all_igrams!(igram_list, [early, late])
        removefirst!(geolist, late)  # Only late, since first was popped off
    end

    return chosen_igrams

end



""" Ideal way to pick stack igrams: evenly spaced"""
function pick_igrams_default(geolist::Array{Date, 1})::Array{Igram, 1}
    geo_pairs = stack_indices(geolist)
    return [(geolist[i1], geolist[i2]) for (i1, i2) in geo_pairs]
end


"""Find the indices of the stack-picked igrams within the full stack
of all compute igrams (that are stored in the unw_stack file)"""
indices_from_full(chosen_igrams, full_igram_list) = indexin(chosen_igrams, full_igram_list)
