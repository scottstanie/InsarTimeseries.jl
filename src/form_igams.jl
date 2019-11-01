take_looks = Sario.take_looks

# TODO: check if preallocating and doing ! version is worth speed
function make_igam(slc1::AbstractArray, slc2::AbstractArray, rowlooks::Int, collooks::Int)
    return take_looks(slc1 .* conj.(slc2), rowlooks, collooks)
end

function powlooks(image::AbstractArray, rowlooks::Int, collooks::Int)
    return sqrt.(take_looks(abs2.(image), rowlooks, collooks))
end

function make_int_cor(slc1::AbstractArray, slc2::AbstractArray, rowlooks::Int, collooks::Int)
    igram = make_igam(slc1, slc2, rowlooks, collooks)
    ampslc1 = powlooks(slc1, rowlooks, collooks)
    ampslc2 = powlooks(slc2, rowlooks, collooks)
    amp = @. real(abs(igram))
    @show typeof(ampslc1), typeof(amp)
    cor = real.(amp ./ (eps(eltype(ampslc1)) .+ ampslc1 .* ampslc2))
    return cor, amp, igram
    # c = amp/(np.finfo(float).eps+slclooks1*slclooks2)
    # cc = amp+1j*c
    # tosave[:,0:nrg] = np.real(c)
    # tosave[:,nrg:] = np.imag(c)
end

# julia> lines = readlines("sbas_list");
# 
# julia> lines[1:2]
# 4-element Array{String,1}:
#  "../S1A_20141104.geo ../S1A_20141128.geo 24.0    29.539676548307892     "
#  "../S1A_20141104.geo ../S1A_20141222.geo 48.0    11.520530983588465     "
function form_igram_names()
    sbas_lines = readlines("sbas_list")
    # TODO: use the parsers to get the dates...
    out = []
    for line in sbas_lines
        early_file, late_file, temp, spatial = split(line)
        # "../S1A_20141104.geo"
        igram_name = join(map(_get_date, [early_file, late_file]), '_') *".int"
        push!(out, (igram_name, String(early_file), String(late_file)))
    end
    # Note: Sorting so that ALL igrams with `early_file` are formed in a for
    return sort(out)
end
function create_igrams(rowlooks=1, collooks=1)
    current_ints = Glob.glob("*.int")
    current_cors = Glob.glob("*.cc")
    cur_early_file = ""
    early = zeros(ComplexF32, 10, 10)
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
            @time early = Sario.load(early_file, do_permute=false)
            cur_early_file = early_file
        end
        # But we load late every time
        println("Loading $late_file")
        @time late = Sario.load(late_file, do_permute=false)

        # transposed, so flip row/col
        println("Forming int, cor")
        @time cor, amp, igram = make_int_cor(early, late, collooks, rowlooks)

        println("Saving $cor_name, $igram_name")
        outcor = cat(cor, amp, dims=3)
        Sario.save(cor_name, outcor, do_permute=false)
        Sario.save(igram_name, igram, do_permute=false)
    end
end

# julia> "../S1A_20141128.geo" |> x-> split(x, '_')[2] |> x -> split(x, '.')[1]
#     "20141128"
_get_date(geo_name) = geo_name |> x-> split(x, '_')[2] |> x -> split(x, '.')[1]

