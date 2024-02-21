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
