# ActorRABBIT whose act params write to beams 

Base.@kwdef mutable struct RABBITtimetraces 
    n_time::Union{Real,Missing} = missing
    n_rho::Union{Real,Missing} = missing
    time::Union{Vector{Float64},Missing} = missing 
    rho::Union{Vector{Float64},Missing} = missing 
    te::Union{Vector{Float64},Missing} = missing

    ti1::Union{Vector{Float64},Missing} = missing
    ti2::Union{Vector{Float64},Missing} = missing
    ti3::Union{Vector{Float64},Missing} = missing
    ti4::Union{Vector{Float64},Missing} = missing
    ti5::Union{Vector{Float64},Missing} = missing
    ti6::Union{Vector{Float64},Missing} = missing # how many ion species does RABBIT support? 

    dene::Union{Vector{Float64},Missing} = missing 
    rot_freq_tor::Union{Vector{Float64},Missing} = missing 
    zeff::Union{Vector{Float64},Missing} = missing
    pnbi::Union{Real,Missing} = missing
 end

# Base.@kwdef mutable struct RABBITbeams

# end

# function FUSEtoRABBITbeams()
    
# end

# function FUSEtoRABBITequilibria()

# end

function FUSEtoRABBITtimetraces(dd::IMAS.dd)
    cp1d = dd.core_profiles.profiles_1d
    ions = cp1d[1].ion
    eV_to_keV = 1e-3

    timetraces = RABBITtimetraces()

    timetraces.n_time = length(cp1d)
    timetraces.n_rho = length(cp1d[1].grid.rho_tor_norm)
    timetraces.time = [cp1d[1].time] # need to think about adding appropriate time slice 
    timetraces.rho = cp1d[1].grid.rho_tor_norm

    timetraces.te = cp1d[1].electrons.temperature * eV_to_keV

    for idx in eachindex(ions)
        species = idx + 1
        setfield!(timetraces, Symbol("ti$species"), ions[idx].temperature * eV_to_keV)
    end

    timetraces.dene = cp1d[1].electrons.density 

    timetraces.rot_freq_tor = cp1d[1].rotation_frequency_tor_sonic
    timetraces.zeff = cp1d[1].zeff
    timetraces.pnbi = dd.nbi.unit[1].power_launched.data[1] # this should probably be a full profile rather than just a single value 
    
    return timetraces
end

function write_timetraces(timetraces::RABBITtimetraces, filename::AbstractString)

    nw = 5
    
    open(filename, "w") do io 
        println(io, "         ", timetraces.n_time)
        println(io, "         ", timetraces.n_rho)
        println(io, "rho_tor")
        println(io, cropdata_f(timetraces.time, nw)) # this should be cropped but here we have only a single time slice rather than many 
        println(io, cropdata_f(timetraces.rho, nw)) # is 5 the correct length of the thing
        println(io, cropdata_f(timetraces.te, nw))

        ### then all the ion temperature profiles ###

        println(io, cropdata_e(timetraces.dene, nw))
        println(io, cropdata_f(timetraces.rot_freq_tor, nw))
        println(io, cropdata_f(timetraces.zeff, nw))
        # println(io, cropdata_f(timetraces.pnbi, nw)) # this should be an array with multiple time slices
    end
end


