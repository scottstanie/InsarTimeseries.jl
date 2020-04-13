import ArchGDAL
function geolocate(
    slc::AbstractArray,
    lons::AbstractArray,
    lats::AbstractArray,
    step_deg = 0.1,
)
    geolons, geolats = coarse_grid_deg(lons, lats, step_deg)
    geoslc = zeros(eltype(slc), size(geolats, 1), size(geolons, 1))
    # Add counts to take avg when multiple pixels fall in same output bucket
    counts = zeros(size(geoslc))

    println("Using step size $step_deg")
    println("Input/output sizes:")
    @show size(lons), size(lats), size(slc)
    @show size(geolons), size(geolats), size(geoslc)
    for idx in eachindex(lons, lats, slc)
        (lon, lat, val) = lons[idx], lats[idx], slc[idx]
        row = nearest(geolats, lat)
        col = nearest(geolons, lon)
        oob(row, col, geoslc) && continue
        geoslc[row, col] += val
        counts[row, col] += 1
    end

    return (geoslc ./ counts), geolons, geolats
end

function coarse_grid_deg(lons, lats, step)
    lon1, lon2 = extrema(lons)
    lat1, lat2 = extrema(lats)
    return lon1:step:lon2, lat2:(-step):lat1
end

# if we have just an array (with even step size) and point:
nearest(arr::AbstractArray, point) = Int(round((point - first(arr)) / (arr[2] - arr[1])))

# Out of bounds check (if we decide to make output area smaller than input)
oob(row, col, out) = row < 1 || row > size(out, 1) || col < 1 || col > size(out, 2)



function loadfiles()
    rows, cols = 7319, 4720
    slc = Array{ComplexF32,2}(undef, (cols, rows))
    lons = Array{Float32,2}(undef, (cols, rows))
    lats = similar(lons)

    read!("10913_11584.int", slc)
    slc = permutedims(slc)

    read!("lon", lons)
    lons = permutedims(lons)

    read!("lat", lats)
    lats = permutedims(lats)
    return (slc, lats, lons)
end

function write_geotiff(geoslc, geolons, geolats, outname = "geoslc.tif")
    height, width = size(geoslc)

    ArchGDAL.create(
        "geoslc.tif",
        driver = ArchGDAL.getdriver("GTiff"),
        width = width,
        height = height,
        nbands = 1,
        dtype = Float32,
    ) do dataset
        # ArchGDAL expects the dimensions of the buffer to be (cols, rows, bands) or (cols, rows).
        ArchGDAL.write!(dataset, permutedims(Float32.(log10.(abs.(geoslc)))), 1)
    end

    lry, uly = extrema(geolats)
    ulx, lrx = extrema(geolons)
    run(`gdal_edit.py geoslc.tif -a_srs EPSG:4326 -a_ullr $ulx $uly $lrx $lry`)
    # xres = geolons[2] - geolons[1]
    # yres = geolats[2] - geolats[1]
    # run(`gdal_edit.py geoslc.tif -a_srs EPSG:4236 -tr $xres $yres -a_ullr $ulx $uly $lrx $lry
end

function testgeocode(step_deg = 0.0002)
    # step_deg = 0.001  # Leads to 935x892
    slc, lats, lons = loadfiles()
    geoslc, geolons, geolats = geolocate(slc, lons, lats, step_deg)
    write_geotiff(geoslc, geolons, geolats, outname = "geoslc.tif")
end
