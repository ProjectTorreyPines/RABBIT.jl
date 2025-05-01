using Printf

function cropdata_f(data::Vector{Float64}, nw::Int)
    text = ""
    if length(data) < nw
        for j in 1:length(data)
            if data[j] < 0
                formatted_value = @sprintf("      %8.7f", data[j])[1:16]  # One less space for negatives
            else
                formatted_value = @sprintf("       %8.7f", data[j])[1:16]
            end
            text *= formatted_value
        end
        text *= "\n"
    else
        for i in 0:floor(Int, (length(data) - 1) / nw)
            start_index = i * nw + 1
            end_index = min((i + 1) * nw, length(data))
            for j in start_index:end_index
                if data[j] < 0
                    formatted_value = @sprintf("      %8.7f", data[j])[1:16]  # One less space for negatives
                else
                    formatted_value = @sprintf("       %8.7f", data[j])[1:16]
                end
                text *= formatted_value
            end
            text *= "\n"
        end
    end
    return text
end

function cropdata_e(data::Vector{Float64}, nw::Int)
    text = ""
    if length(data) < nw
        for j in 1:length(data)
            text *= @sprintf("   %.7e", data[j])
        end
        text *= "\n"
    else
        for i in 0:floor(Int, (length(data) - 1) / nw)
            start_index = i * nw + 1
            end_index = min((i + 1) * nw, length(data))
            for j in start_index:end_index
                text *= @sprintf("   %.7e", data[j])
            end
            text *= "\n"
        end
    end
    return text
end

function cropdata_f(data::Vector{Vector{Float64}}, nw::Int)
    text = ""
    for i in 1:length(data)
        text *= cropdata_f(data[i], nw)
    end 
    return text 
end 

function cropdata_e(data::Vector{Vector{Float64}}, nw::Int)
    text = ""
    for i in 1:length(data)
        text *= cropdata_e(data[i], nw)
    end 
    return text 
end 

function cropdata_f(data::Matrix{Float64}, nw::Int)
    dataf = data[:]
    cropdata_f(dataf, nw)
end 

function print6(d)
    text = ""
    i0 = 1
    n = length(d)
    for i = 1:div(n, 6)
        row = [ @sprintf("%13.6f", d[j]) for j = i0:i0+5 ]
        text *= join(row, "") * "\n"
        i0 += 6
    end
    for i = i0:n
        text *= @sprintf("%13.6f", d[i])
    end
    text *= "\n"
    return text
end

function pretty_print_vector(data::Vector{Float64})
    text = "    "
    for i in 1:length(data)
        text *= string(data[i])
        text *= "       "
    end
    return text
end

function struct_unpack(data::Vector{UInt8})
    return reinterpret(Float32, data)[1]
end