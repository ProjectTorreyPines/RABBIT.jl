using Printf

function cropdata_f(data::Vector{Float64}, nw::Int)
    text = ""
    if length(data) < nw
        for j in 1:length(data)
            text *= @sprintf("%16.7f", data[j])
        end
        text *= "\n"
    else
        for i in 0:floor(Int, (length(data) - 1) / nw)
            start_index = i * nw + 1
            end_index = min((i + 1) * nw, length(data))
            for j in start_index:end_index
                text *= @sprintf("%16.7f", data[j])
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
            text *= @sprintf("%16.7e", data[j])
        end
        text *= "\n"
    else
        for i in 0:floor(Int, (length(data) - 1) / nw)
            start_index = i * nw + 1
            end_index = min((i + 1) * nw, length(data))
            for j in start_index:end_index
                text *= @sprintf("%16.7e", data[j])
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