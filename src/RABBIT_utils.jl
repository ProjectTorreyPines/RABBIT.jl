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
            text *= "\n" # need to remove the new line when it's the last value i.e. don't end with whitespace
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

function get_cp1d_time_slice(dd::IMAS.dd, field::Symbol, time_slice::Int)
    cp1d = dd.core_profiles.profiles_1d
    res = getfield(IMAS.freeze(cp1d[time_slice]), field)
    return res
end