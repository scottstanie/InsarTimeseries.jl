prunes = [false, false, false, true, true]
L1arr = [false, true, false, true, false]
constant_vel = [false, true, true, true, true]  # First is just stack avg
stackavgarr = .!constant_vel


for (doprune, isl1, constvel, stackavg) in zip(prunes, L1arr, constant_vel, stackavgarr)
    soltype = isl1 ? "l1" : (stackavg ? "stackavg" : "l2")
    prunestr = doprune ? "prune" : "noprune"
    outfile = "velocities_$(prunestr)_$soltype.h5"
    @show outfile

    try
        @time InsarTimeseries.run_inversion(
            constant_velocity = constvel,
            stack_average = stackavg,
            max_temporal_baseline = 500,
            ignore_geo_file = "geolist_ignore.txt",
            L1 = isl1,
            prune = doprune,
            outfile = outfile,
        )
    catch
        println("Skipping $outfile")
        continue
    end

    # stackavg && continue
    # try
    #     InsarTimeseries._save_std_attrs(outfile, "velos/1")
    # catch
    #     print()
    # end
    # for ds in ["stddev_raw/1", "stddev/1"]
    #     tmp = h5read(outfile, ds)
    #     h5open(outfile, "cw") do f
    #         tmp = abs(InsarTimeseries.PHASE_TO_CM) .* read(f, ds)
    #         o_delete(f[ds])
    #         f[ds] = tmp
    #     end
    # end
end
