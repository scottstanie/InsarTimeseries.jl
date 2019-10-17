import Glob

function proc_subset(lats, lons, out_dir; ref_station=nothing, ref_row=nothing, ref_col=nothing)
    intf = Glob.glob("*.int")[1]  # Dummy file
    i1 = MapImages.MapImage(intf)

    # demrsc = Sario.load("dem.rsc")
    cropped = i1[lats, lons]
    Sario.save(joinpath(out_dir, "dem.rsc"), cropped.demrsc)

    # Now get the range for cropping all the other files
    rowrange, colrange = MapImages.latlon_to_rowcol(i1.demrsc, lats, lons)
    
    @show rowrange, colrange, lats, lons
    println("New dem:")
    @show cropped.demrsc
    println("Ints:")
    @time InsarTimeseries.loop_over_files(x->x[rowrange, colrange], ".", ".int", ".int", out_dir=out_dir, overwrite=false, do_permute=true);
    println("ccs:")
    @time InsarTimeseries.loop_over_files(x->x[rowrange, colrange], ".", ".cc", ".cc", out_dir=out_dir, overwrite=false, do_permute=true);
    println("geo masks:")
    @time InsarTimeseries.loop_over_files(x->x[rowrange, colrange], ".", ".geo.mask", ".geo.mask", out_dir=out_dir, overwrite=false, do_permute=true);

    cd(out_dir)
    run(`insar process --step 8`)
    @time InsarTimeseries.prepare_stacks(".", ref_station=ref_station, ref_row=ref_row, ref_col=ref_col, overwrite=false, zero_masked=false)
end
