take_looks = Sario.take_looks

function make_igam(slc1::AbstractArray, slc2::AbstractArray, rowlooks::Int, collooks::Int)
    return take_looks(slc1 .* conj.(slc2), rowlooks, collooks)
end

function powlooks(image::AbstractArray, rowlooks::Int, collooks::Int)
    return sqrt.(take_looks(abs2.(image), rowlooks, collooks))
end

function correlation(slc1::AbstractArray, slc2::AbstractArray, rowlooks::Int, collooks::Int)
    return correlation(make_igam(slc1y, slc2, rowlooks, collooks))
end

function igram_correlation(igram::AbstractArray, rowlooks::Int, collooks::Int)
    ampslc1 = powlooks(slc1,rowlooks,collooks)
    ampslc2 = powlooks(slc2,rowlooks,collooks)
    amp = @. real(abs(igram))
    cor = real.(amp ./ eps(eltype(igram) .+ ampslc1 .* ampslc2))
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
function create_igrams()
    current_ints = Glob.glob("*.int")

    sbas_lines = readlines("sbas_list")
    # TODO: use the parsers to get the dates...
    for line in sbas_lines[1:3]
        early, late, temp, spatial = split(line)
        # "../S1A_20141104.geo"
        igram_name = join(map(_get_date, early, late), '_') *".int"
    end
    return
end

# julia> "../S1A_20141128.geo" |> x-> split(x, '_')[2] |> x -> split(x, '.')[1]
#     "20141128"
_get_date(geo_name) = geo_name |> x-> split(x, '_')[2] |> x -> split(x, '.')[1]

