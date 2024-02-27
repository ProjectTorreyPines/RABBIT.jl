Base.@kwdef mutable struct RABBITtimetraces 
    n_time::Union{Real,Missing} = missing
    n_rho::Union{Real,Missing} = missing
    time::Union{Vector{Float64},Missing} = missing 
    rho::Union{Vector{Float64},Missing} = missing 
    te::Union{Vector{Vector{Float64}},Missing} = missing
    ti::Union{Vector{Vector{Float64}},Missing} = missing
    dene::Union{Vector{Vector{Float64}},Missing} = missing
    rot_freq_tor::Union{Vector{Vector{Float64}},Missing} = missing
    zeff::Union{Vector{Vector{Float64}},Missing} = missing
    pnbi::Union{Array{Float64},Missing} = missing
end 

# Base.@kwdef mutable struct RABBITbeams

# end

# function FUSEtoRABBITbeams()
    
# end

# function FUSEtoRABBITequilibria()

# end

function FUSEtoRABBITtimetraces(dd::IMAS.dd)
    cp1d = dd.core_profiles.profiles_1d
    eV_to_keV = 1e-3
    cm3_to_m3 = 1e-6

    timetraces = RABBITtimetraces()

    timetraces.n_time = length(cp1d)
    timetraces.n_rho = length(cp1d[1].grid.rho_tor_norm)
    timetraces.rho = cp1d[1].grid.rho_tor_norm

    timetraces.time = [cp1d[i].time for i in eachindex(cp1d)] .* 1e3
    timetraces.te = [cp1d[i].electrons.temperature .* eV_to_keV for i in eachindex(cp1d)]
    timetraces.dene = [get_cp1d_time_slice(dd, :electrons, i).density .* cm3_to_m3 for i in eachindex(cp1d)]
    timetraces.rot_freq_tor = [cp1d[i].rotation_frequency_tor_sonic for i in eachindex(cp1d)]
    timetraces.zeff = [cp1d[i].zeff for i in eachindex(cp1d)]

    ion_temps = []

    for j in 1:length(cp1d)
        for i in 1:length(cp1d[1].ion)
            push!(ion_temps, cp1d[j].ion[i].temperature .* eV_to_keV)
        end
    end
    
    timetraces.ti = ion_temps

    pnbis = []
    for i in 1:length(dd.nbi.unit)
        for j in 1:length(dd.nbi.unit[1].power_launched.data)
            push!(pnbis, dd.nbi.unit[i].power_launched.data[j])
        end
    end
    timetraces.pnbi = pnbis

    return timetraces
end

function write_timetraces(timetraces::RABBITtimetraces, filename::AbstractString)
    nw = 5
    
    open(filename, "w") do io 
        println(io, "         ", timetraces.n_time)
        println(io, "         ", timetraces.n_rho)
        println(io, "rho_tor")
        print(io, cropdata_f(timetraces.time, nw))
        print(io, cropdata_f(timetraces.rho, nw))
        print(io, cropdata_f(timetraces.te, nw))
        print(io, cropdata_f(timetraces.ti, nw))
        print(io, cropdata_e(timetraces.dene, nw))
        print(io, cropdata_f(timetraces.rot_freq_tor, nw))
        print(io, cropdata_f(timetraces.zeff, nw))
        print(io, cropdata_f(timetraces.pnbi, nw))
    end
end


