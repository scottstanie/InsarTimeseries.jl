try
    rm("test1.h5")
catch
end
try
    rm("test2.h5")
catch
end
function write_in_layers_repack(filename = "test1.h5", dset_name = "test")
    file_list = Glob.glob("*.unw")
    stack_size = (720, 720, length(file_list))

    h5open(filename, "cw") do f
        d_create(
            f,
            dset_name,
            datatype(Float32),
            dataspace(stack_size),
            # "chunk", (1, 1, size(stack_in, 3)),  # Seems better to repack?
        )
    end


    println("Reading layers and writing")
    arr_buf = zeros(Float32, (720, 720))
    h5open(filename, "cw") do hf
        dset = hf[dset_name]
        for (idx, f) in enumerate(file_list)
            arr_buf .= Sario.load(f, do_permute = true)
            dset[:, :, idx] = arr_buf
        end
    end
    InsarTimeseries.repack(filename, dset_name)

    return
end

function write_in_pixel_chunks(filename = "test2.h5", dset_name = "test")
    file_list = Glob.glob("*.unw")
    stack_size = (720, 720, length(file_list))
    m = Array{Float32,3}(undef, stack_size)
    println("Reading all layers into memory")
    for (idx, f) in enumerate(file_list)
        m[:, :, idx] .= Sario.load(f)
    end

    h5open(filename, "cw") do f
        d_create(
            f,
            dset_name,
            datatype(Float32),
            dataspace(stack_size),
            "chunk",
            (1, 1, size(m, 3)),  # Seems better to repack?
        )
    end

    println("writing pixelwise")
    h5open(filename, "r+") do f
        # NOTE: this seems to take 10x slower than repack
        # for col in 1:size(m, 1)
        #     for row in 1:size(m, 2)
        #         f[dset_name][row, col, :] = m[row, col, :]
        #     end
        # end
        f[dset_name][:, :, :] = m

    end
    return
end

@time write_in_layers_repack();
@time write_in_pixel_chunks();
