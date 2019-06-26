import Glob

function read_stack(file_ext::String, rows::Int, cols::Int)
	allfiles = Glob.glob(file_ext, ".")
	buffer = Array{Float32, 2}(undef, (rows, cols))
    stack = Array{Float32, 3}(undef, (rows, cols, length(allfiles)))
    for (idx, f) in enumerate(allfiles)
 		read!(f, buffer)
 	   	stack[:, :, idx] = buffer
    end
    return stack
end
