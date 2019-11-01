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

