
# I hope this is correct!
# This was generated by Claude 2 based on the Octave C++ implementation at:
# https://github.com/shsajjadi/OctaveUnimplemented/blob/master/graydist.cc 
"""
    graydist(I, seeds; metric="quasi-euclidean")

Compute gray weighted distance transform of image I using quasi-euclidean metric.
seeds is a boolean mask indicating seed points. 
"""
function graydist(I::AbstractMatrix, seeds::AbstractMatrix; metric="quasi-euclidean")

    # Check inputs
    size(I) == size(seeds) || error("I and seeds must have same size")

    metrics = ["quasi-euclidean"]
    metric ∈ metrics || error("Invalid metric $metric")

    T = similar(I, Float32) 
    fill!(T, Inf)
    
    idx = similar(seeds, Int)
    
    root2 = sqrt(2)
    h = 0.5
    
    shifts = ((-1,0), (1,0), (0,-1), (0,1), (-1,-1), (-1,1), (1,-1), (1,1))
    
    fringe = CartesianIndex{2}[]
    for i in CartesianIndices(seeds)
        if seeds[i]
            push!(fringe, i)
            T[i] = 0
            idx[i] = 0
        end 
    end
    
    while !isempty(fringe)
        u = pop!(fringe)
        
        fu = I[u]
        for shift in shifts
            v = u + CartesianIndex(shift)
            if checkbounds(Bool, I, v) && T[v] > T[u] + cost(fu, I[v], u, v)
                T[v] = T[u] + cost(fu, I[v], u, v)
                idx[v] = idx[u]
                push!(fringe, v)
            end
        end
    end
    
    return T, idx
end

function cost(fu, fv, u, v)
    x1, y1 = Tuple(u)
    x2, y2 = Tuple(v)
    
    if abs(x1 - x2) > abs(y1 - y2)
        return h*(fu + fv) + (root2 - 1)*abs(y1 - y2)
    else
        return (root2 - 1)*abs(x1 - x2) + h*(fu + fv) 
    end
end





####### Old Stuff


# function graydist(I::AbstractMatrix, seeds::AbstractMatrix; metric="quasi-euclidean")
    
#     # Check inputs
#     size(I) == size(seeds) || error("I and seeds must have same size")
    
#     metrics = ["quasi-euclidean"]
#     metric ∈ metrics || error("Invalid metric $metric")
    
#     T = similar(I, Float32)
#     fill!(T, Inf)
#     idx = similar(seeds, Int)
    
#     h = 0.5
#     root2 = sqrt(2)
    
#     shifts = ((-1, 0), (1, 0), (0, -1), (0, 1))
        
#     fringe = CartesianIndex{2}[]
#     for i in CartesianIndices(seeds)
#         if seeds[i]
#             push!(fringe, i)
#             T[i] = 0
#             idx[i] = 0
#         end
#     end
    
#     while !isempty(fringe)
#         u = pop!(fringe)
        
#         fu = I[u]
#         for shift in shifts
#             v = u + CartesianIndex(shift)
#             if checkbounds(Bool, I, v) && (T[v] > T[u] + metric_cost(fu, I[v], metric))
#                 T[v] = T[u] + metric_cost(fu, I[v], metric)
#                 idx[v] = idx[u]
#                 push!(fringe, v) 
#             end
#         end
#     end
    
#     return T, idx
# end

# function metric_cost(fu, fv, metric)
#     if metric == "quasi-euclidean"
#         return h*(fu + fv) + root2*abs(fu - fv)
#     else
#         error("Invalid metric $metric") 
#     end
# end




# function create_quasi_euclidean_chamfer_weights()
#     h = 0.5
#     sq2 = √0.5
#     weight_quasi = Vector{Vector{Float64}}([
#         Float64[],
#         [h, sq2, h],
#         [h, h, sq2, h, sq2],
#         [h, sq2, h],
#         [h, sq2, h, h, sq2],
#         [sq2, h, sq2, h, h, sq2, h, sq2],
#         [sq2, h, h, sq2, h],
#         [h, sq2, h],
#         [sq2, h, sq2, h, h],
#         [sq2, h, h]
#     ])
#     return weight_quasi
# end


# using Images

# function graydist_cpp(img::Array{T,N}, seeds::Array{Bool,N}; metric::String="quasi-euclidean") where {T,N}

#     ResultType = promote_type(Float32, T)
#     IndexType = UInt32

#     dist = fill(typemax(ResultType), size(img))
#     idx_segment = similar(dist, IndexType)
#     idx_predecessor = similar(dist, IndexType)

#     offset = get_neighbor_offsets(N)
#     weights = create_quasi_euclidean_chamfer_weights(N)

#     Q = PriorityQueue{Tuple{CartesianIndex{N}, ResultType}}()

#     for i in CartesianIndices(seeds)
#         if seeds[i]
#             dist[i] = zero(ResultType)
#             idx_segment[i] = LinearIndices(img)[i]
#             idx_predecessor[i] = 0
#             enqueue!(Q, (i, zero(ResultType)))
#         end
#     end

#     while !isempty(Q)
#         u, udist = dequeue!(Q)

#         for (i,offset) in enumerate(offset)
#             v = u + offset
#             if checkbounds(Bool, img, v) && !isassigned(dist, v)
#                 vdist = udist + weights[i] * (img[u] + img[v])
#                 if vdist < dist[v]
#                     dist[v] = vdist
#                     idx_segment[v] = idx_segment[u]
#                     idx_predecessor[v] = u
#                     enqueue!(Q, (v, vdist))
#                 end
#             end
#         end
#     end

#     return dist, idx_segment, idx_predecessor
# end

# function get_quasi_euclidean_weights(N)
#     weights = [sqrt(0.5) for i in 1:2N]
#     return weights
# end

# function get_neighbor_offsets(N)
#     offsets = [CartesianIndex(ntuple(d -> d==i ? [-1,1] : [0], Val(N))) for i in 1:N]
#     return offsets
# end



