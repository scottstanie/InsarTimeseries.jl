using BenchmarkTools

function mycopy(x::Vector{T}) where T
   inds = axes(x, 1)
   out = similar(Array{T}, inds, inds, inds)
   for j = inds
       for i = inds
           for k = inds
               out[i,j,k] = x[k]
           end
       end
   end
   return out
end


function copy0(x::Vector{T}) where T
   inds = axes(x, 1)
   out = similar(Array{T}, inds, inds, inds)
   for k = inds
       for j = inds
           for i = inds
               out[i,j,k] = x[i]
           end
       end
   end
   return out
end
function copy01(x::Vector{T}) where T
   inds = axes(x, 1)
   out = similar(Array{T}, inds, inds, inds)
   for k = inds
       for j = inds
           out[:,j,k] = x
       end
   end
   return out
end

function copy01dot(x::Vector{T}) where T
   inds = axes(x, 1)
   out = similar(Array{T}, inds, inds, inds)
   for k = inds
       for j = inds
           out[:,j,k] .= x
       end
   end
   return out
end

function copy1(x::Vector{T}) where T
   inds = axes(x, 1)
   out = similar(Array{T}, inds, inds, inds)
   for i = inds
       for j = inds
           out[i,j,:] .= x
       end
   end
   return out
end

function copy2(x::Vector{T}) where T
   inds = axes(x, 1)
   out = similar(Array{T}, inds, inds, inds)
   for j = inds
       for i = inds
           out[i,j,:] .= x
       end
   end
   return out
end

function copy3(x::Vector{T}) where T
   lx = length(x)
   out = similar(Array{T}, lx, lx*lx)
   for j = 1:lx
       out[:,j] .= x
   end
   return reshape(out, (length(x), length(x), length(x)))
end

function copy4(x::Vector{T}) where T
   lx = length(x)
   out = similar(Array{T}, lx*lx, lx)
   for j = 1:lx
       out[j,:] = x
   end
   return reshape(out, (length(x), length(x), length(x)))
end


function copy5(x::Vector{T}) where T
   lx = length(x)
   out = similar(Array{T}, lx,lx,lx)
   for idx in eachindex(out)
       out[idx] = x[idx % lx]
   end
   return out
end

x = randn(500);
println("mycopy"); @btime mycopy(x);
println("copy0"); @btime copy0(x);
println("copy01"); @btime copy01(x);
println("copy01dot"); @btime copy01dot(x);
# println("copy1"); @btime copy1(x);
println("copy2"); @btime copy2(x);
println("copy3"); @btime copy3(x);
# @btime copy4(x);
# @println("copy5"); btime copy5(x);
