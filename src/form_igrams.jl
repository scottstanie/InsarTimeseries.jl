import Sario
import Sario: take_looks, take_looks!
import Glob: glob

function make_igam(slc1::AbstractArray, slc2::AbstractArray, rowlooks::Int, collooks::Int)
    return take_looks(slc1 .* conj.(slc2), rowlooks, collooks)
end

function make_igam!(
    igram::AbstractArray,
    slc1::AbstractArray,
    slc2::AbstractArray,
    rowlooks::Int,
    collooks::Int,
)
    # We are assuming it's okay to overwrite slc2
    slc2 .*= conj.(slc1)
    slc2 .= conj.(slc2)
    return take_looks!(igram, slc2, rowlooks, collooks)
end


function powlooks(image::AbstractArray, rowlooks::Int, collooks::Int)
    return take_looks(abs2.(image), rowlooks, collooks)
end

function powlooks!(out::AbstractArray, image::AbstractArray, rowlooks::Int, collooks::Int)
    return take_looks!(out, abs2.(image), rowlooks, collooks)
end

function make_int_cor(
    slc1::AbstractArray,
    slc2::AbstractArray,
    rowlooks::Int,
    collooks::Int,
)
    igram = make_igam(slc1, slc2, rowlooks, collooks)
    return _make_cor(igram, slc1, slc2, rowlooks, collooks)
end

function _make_cor(
    igram::AbstractArray,
    slc1::AbstractArray,
    slc2::AbstractArray,
    rowlooks::Int,
    collooks::Int,
)
    ampslc1 = sqrt.(powlooks(slc1, rowlooks, collooks))
    ampslc2 = sqrt.(powlooks(slc2, rowlooks, collooks))
    amp = @. real(abs(igram))
    # @show extrema(ampslc1), extrema(ampslc2), extrema(amp)
    cor = real.(amp ./ (eps(Float32) .+ (ampslc1 .* ampslc2)))
    return cor, amp, igram
end

# julia> lines = readlines("sbas_list")
# 4-element Array{String,1}:
#  "../S1A_20141104.geo ../S1A_20141128.geo 24.0    29.539676548307892     " ...
function form_igram_names()
    sbas_lines = readlines("sbas_list")
    # TODO: use the parsers to get the dates...
    out = []
    for line in sbas_lines
        early_file, late_file, temp, spatial = split(line)
        # "../S1A_20141104.geo"
        igram_name = join(map(_get_date, [early_file, late_file]), '_') * ".int"
        push!(out, (igram_name, String(early_file), String(late_file)))
    end
    # Note: Sorting so that ALL igrams with `early_file` are formed in a for
    return sort(out)
end

function create_igrams(rowlooks = 1, collooks = 1)
    current_ints = glob("*.int")
    current_cors = glob("*.cc")

    fulldemrsc = Sario.load("../elevation.dem.rsc")
    # New version
    demrsc = take_looks(fulldemrsc, rowlooks, collooks)
    Sario.save("dem.rsc", demrsc)

    println("Preallocating arrays")
    # Pre-allocate arrays for speed
    # Note transpose of cols/rows, since we'll read in
    # without transposing the data
    fullrows, fullcols = size(fulldemrsc)
    early = zeros(ComplexF32, fullcols, fullrows)
    late = zeros(ComplexF32, fullcols, fullrows)
    rows, cols = size(demrsc)
    igram = zeros(ComplexF32, cols, rows)

    ampslc1 = zeros(Float32, cols, rows)
    # ampslc2 = similar(ampslc1)
    ampslc2 = zeros(Float32, cols, rows)
    amp = similar(ampslc1)
    cor = similar(ampslc1)
    outcor = zeros(Float32, (cols, rows, 2))

    cur_early_file = ""
    for (igram_name, early_file, late_file) in form_igram_names()
        cor_name = replace(igram_name, ".int" => ".cc")
        if (igram_name in current_ints) && (cor_name in current_cors)
            println("Skipping $igram_name and $cor_name: exists")
            continue
        else
            println("Forming $igram_name and $cor_name")
        end

        # Keep early in memory for all pairs: only load for new set
        if cur_early_file != early_file
            println("Loading $early_file")
            @time read!(early_file, early)
            cur_early_file = early_file
        end
        # But we load late every time
        println("Loading $late_file")
        @time read!(late_file, late)

        # @show extrema(abs.(late)), extrema(abs.(early))
        # transposed, so flip row/col
        # println("Forming int, cor")
        # @time cor, amp, igram = make_int_cor(early, late, collooks, rowlooks)

        # TODO: see about ArchGdal operating on blocks for low memory footprint

        # Note again: flipped collooks/rowlooks since we read in as transposed
        println("Forming amps")
        powlooks!(ampslc1, early, collooks, rowlooks)
        ampslc1 .= sqrt.(ampslc1)

        powlooks!(ampslc2, late, collooks, rowlooks)
        ampslc2 .= sqrt.(ampslc2)
        println("Forming igram")
        @time igram .= make_igam!(igram, early, late, collooks, rowlooks)

        println("Forming cor")
        amp .= real.(abs.(igram))
        # @show extrema(ampslc1), extrema(ampslc2), extrema(amp)
        cor .= real.(amp ./ (eps(Float32) .+ (ampslc1 .* ampslc2)))

        println("Saving $cor_name, $igram_name")
        outcor[:, :, 1] .= amp
        outcor[:, :, 2] .= cor
        Sario.save(cor_name, outcor, do_permute = false)
        Sario.save(igram_name, igram, do_permute = false)
    end
    GC.gc()
end

# julia> "../S1A_20141128.geo" |> x-> split(x, '_')[2] |> x -> split(x, '.')[1]
#     "20141128"
_get_date(geo_name) = geo_name |> x -> split(x, '_')[2] |> x -> split(x, '.')[1]
