module JankoUtils

export read_dict, write_dict, load_float64_array, pprint, diff_dict, idxmod, inspect, scale_cm

import SparseArrays
const SA = SparseArrays

"""
    diff_dict(old_dict, new_dict)

highlight the differences between two Dict{String, Float64} objects
- added keys are shown in green
- removed keys are shown in red
- modified values are shown in bold
"""
function diff_dict(old_dict::Dict{String, Float64}, new_dict::Dict{String, Float64})
    keys1 = collect(keys(old_dict))
    keys2 = collect(keys(new_dict))

    plus = setdiff(keys2, keys1)
    minus = setdiff(keys1, keys2)
    common = intersect(keys1, keys2)
    changed = Array{String,1}()
    allkeys = sort(union(keys1, keys2), lt=(x,y)->lowercase(x) < lowercase(y))

    for k in common
        if old_dict[k] != new_dict[k]
            push!(changed, k)
        end
    end

    for k in allkeys
        color = :normal

        d = old_dict
        if !(k in keys1)
            d = new_dict
        end

        if k in plus
            color = :green
        elseif k in minus
            color = :red
        elseif k in changed
            color = :bold
        end
        if color == :normal
            s = string(lpad(k, 15, " "), ": ", d[k], "\n")
        else
            s = string(lpad(k, 15, " "), ": ", old_dict[k], " -> ", new_dict[k], "\n")
        end
        print_with_color(color, s)
    end
end

"""
    pprint(dict)

pretty print a Dict object
"""
function pprint(d::Dict{String, Float64})
    for k in sort(collect(keys(d)), lt=(x,y)->lowercase(x) < lowercase(y))
        println(lpad(k, 15, " "), ": ", d[k])
    end
end

"""
    read_dict(filename)

read a Dict{String, Float64} object from a file
each line of the file should have one key/value pair per line,
separated by a single space
"""
function read_dict(filename::String)
    if !isfile(filename)
        throw(string(filename, " is not a file"))
    end
    d = Dict{String, Float64}()
    f = open(filename)
    l = readline(f)
    while !isempty(l)
        (k, v) = split(l)
        d[k] = parse(Float64, v)
        l = readline(f)
    end
    close(f)
    return d
end

"""
    write_dict(filename, dict)

write a Dict{String, Float64} object to a file, in the format
readable by the read_dict function above
"""
function write_dict(filename::String, d::Dict{String, Float64})
    f = open(filename, "w+")
    for (k, v) in d
        write(f, string(k, " ", v, "\n"))
    end
    close(f)
end

"""
    load_float64_array(dirname, filename, cols)

wrapper around load_float64_array
"""
function load_float64_array(dirname::String, filename::String, cols::Int=1)
    return load_float64_array(string(dirname, "/", filename), cols)
end

"""
    load_float64_array(filename, cols)

wrapper around `read` to read an Array{Float64} object from a file
with the specified number of columns
"""
function load_float64_array(filename::String, cols::Int=1)::Array{Float64}
    s = div(filesize(filename), 8)
    array = read(filename, Float64, s)
    if cols > 1
        return reshape(array, (div(size(array,1), cols), cols))
    end
    return array
end

"""
    bin_to_txt_array(src, dst, cols)

convert an array stored in binary form to text form
"""
function bin_to_txt_array(source::String, dest::String, cols::Int=1)
    error("not checked")
    a = load_float64_array(source, cols)
    f = open(dest, "w")
    writedlm(f, a)
    close(f)
end

"""
    flatten_tuple_array(array, t)

convert an array of tuples to a Array{t,2} object
"""
function flatten_tuple_array(a::Array, t::Type)
    if isempty(a)
        return Array{t}()
    end
    n = size(a,1)
    m = length(a[1])
    return reinterpret(t, a, (m, n))'
end

"""
    idxmod(i, n)

helper function to wrap indices around the bounds of an array of length n
"""
function idxmod(i, n::Int)
    if 1<=i<=n
        return i
    else
        return mod(i-1, n)+1
    end
end

macro var_name(x)
    string(x)
end

"""
    inspect(x)

print the name and the value of variable x
"""
function inspect(x)
    println(@var_name(x), ": ", x)
    return
end

"""
    spdiagm_const(coeffs, diag_idx, N)

Constructs a sparse N-by-N matrix, where the diagonal number diag_idx[i] is filled with value coeffs[i].
The diagonal is wrapped to extend over the diagonal -sign(diag_idx[i])*(N-abs(diag_idx[i])).
This is useful for the approximation of differential operators on functions defined on the torus.
"""
function spdiagm_const(coeffs::Array{Float64,1}, diag_idx::Array{Int64,1}, N::Int64)
    if length(coeffs) != length(diag_idx)
        throw("the number of diagonals does not match the number of indices")
    end
    pairs = ()
    for n in 1:length(coeffs)
        i, c = diag_idx[n], coeffs[n]
        conjug = -sign(i)*(N-abs(i))
        if i == 0
            pairs = (pairs..., i => c*ones(N))
            # d_idx = (d_idx..., i)
            # d = (d..., c*ones(N))
            continue
        end

        pairs = (pairs..., i => c*ones(N-abs(i)), conjug => c*ones(N-abs(conjug)))
        # d_idx = (d_idx..., i, conjug)
        # d = (d..., c*ones(N-abs(i)), c*ones(N-abs(conjug)))
    end
    return SA.spdiagm(pairs...)
end

"""
    scale_cm(x, cmap; range_min=r_min, range_max=r_max)

    Returns the values of x mapped through cmap after rescaling to the [0, 256] interval
    If range_min or range_max are specified, the scaling is performed with these values
    instead of minimum(x) and maximum(x), respectively.
"""

function scale_cm(x::Array{Float64,1}, cm;
                  range_min::Union{Float64, Missing}=missing,
                  range_max::Union{Float64, Missing}=missing)
    x_min = ismissing(range_min) ? minimum(x) : range_min
    x_max = ismissing(range_max) ? maximum(x) : range_max

    x_mapped = @. 1.0/(x_max-x_min)*(x-x_min)

    return cm(x_mapped)
end

"""
check if the nodes of a polygon are ordered
in a counter clock wise orientation
"""
function check_ccw_polygon(x::Array{Float64,2})
    n = size(x, 1)
    cm = sum(x; dims=1)/n
    α = angle.(x[:,1] .- cm[1] .+ (x[:,2] .- cm[2])*im)

    total_angle = 0.0

    for i = 1:n
        if i < n
            d = α[i+1] - α[i]
        else
            d = α[1] - α[i]
        end
        if d >= π
            total_angle += d - 2π
        elseif d <= -π
            total_angle += d + 2π
        else
            total_angle += d
        end
    end

    return total_angle >= 0
end

end # module JankoUtils
