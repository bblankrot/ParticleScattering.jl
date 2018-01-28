function binarySearch(f, key, low, high)
    #given monotonically increasing function f, find lowest index i such that
    #f(i) >= key. For monotonically decreasing simply negate f and key, and this
    #will return lowest i such that f(i) <= key
    low > high && error("binarySearch: low > high")
    val = 0
    mid = 0
    while low != high
        #compute function
        mid = div(low + high, 2);
        val = f(mid)
        if val <= key
            low = mid + 1
        else
            high = mid
        end
    end
    if low > mid
        return low,-f(low) #have to recompute in case we're done after "low=mid+1"
    else
        return low,-val
    end
end

function pLeftOfLine(p, l1, l2)
    #tests if a point p is left (>0), on (0) or right (<0) of line between l1
    #and l2
    return ((l2[1] - l1[1])*(p[2] - l1[2])
            - (p[1] - l1[1])*(l2[2] - l1[2]))
end

function pInPolygon(p, ft)
    #use winding number algorithm to test if a point p is in polygon a.
    #returns 1 if it is in, -1 if out, 0 if on.
    #TODO: Currently returns false if it is on a corner, fix this...
    wn = Int64(0)
    n = size(ft,1)
    for i = 1:n
        i_1 = i == n ? 1 : i + 1
        if (ft[i,2] <= p[2])
            #if point is between points of line, and to the left of it, +1
            if (ft[i_1,2] > p[2])
                ret = pLeftOfLine(p, ft[i,:], ft[i_1,:])
                (ret > 0) && (wn += 1)
            end
        else
            #if point is between points of line, and to the right of it, -1
            if (ft[i_1,2] <= p[2])
                ret = pLeftOfLine(p, ft[i,:], ft[i_1,:])
                (ret < 0) && (wn -= 1)
            end
        end
    end
    return (wn != 0)
end

function my_uniqueind(v::Vector{T}) where T <: Number
    #returns inds,vals such that vals[inds[i]] == v[i]
    unique_inds = Array{Int64}(0)
    unique_vals = Array{T}(0)
    k = 0
    for val in v
        ind = findfirst(unique_vals,val)
        if ind == 0
            k += 1
            push!(unique_vals, val)
            push!(unique_inds, k)
        else
            push!(unique_inds, ind)
        end
    end
    return unique_inds,unique_vals
end