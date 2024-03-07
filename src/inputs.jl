using IMAS 
using Printf

Base.@kwdef mutable struct RABBITinput
    time::Union{Float64,Missing} = missing 
    rho::Union{Vector{Float64},Missing} = missing 

    # timetraces 
    n_time::Union{Int,Missing} = missing
    n_rho::Union{Int,Missing} = missing
    te::Union{Vector{Float64},Missing} = missing
    ti::Union{Vector{Float64},Missing} = missing
    dene::Union{Vector{Float64},Missing} = missing
    rot_freq_tor::Union{Vector{Float64},Missing} = missing
    zeff::Union{Vector{Float64},Missing} = missing
    pnbi::Union{Array{Float64},Missing} = missing

    # equilibria 
    nw::Union{Int,Missing} = missing
    nh::Union{Int,Missing} = missing
    psirz::Union{Matrix{Float64},Missing} = missing
    npsi1d::Union{Int,Missing} = missing
    qpsi::Union{Vector{Float64},Missing} = missing
    fpol::Union{Vector{Float64},Missing} = missing
    sibry::Union{Float64,Missing} = missing
    simag::Union{Float64,Missing} = missing
    signip::Union{Float64,Missing} = missing
    rmaxis::Union{Float64,Missing} = missing
    zmaxis::Union{Float64,Missing} = missing

    # AuxQuantities
    r::Union{Vector{Float64},Missing} = missing
    z::Union{Vector{Float64},Missing} = missing
    rhorz::Union{Matrix{Float64},Missing} = missing
    psi::Union{Vector{Float64},Missing} = missing
    vol::Union{Vector{Float64},Missing} = missing
    area::Union{Vector{Float64},Missing} = missing

    # beams 
end

function FUSEtoRABBITinput(dd::IMAS.dd)
    eV_to_keV = 1e-3
    cm3_to_m3 = 1e-6

    eq = dd.equilibrium

    all_inputs = RABBITinput[]

    for (i,eqt) in enumerate(eq.time_slice)
        time = eqt.time

        eqt2d = findfirst(:rectangular, eqt.profiles_2d)
        if eqt2d === nothing
            continue
        end

        inp = RABBITinput()
        inp.time = time * 1e3

        inp.nw = length(eqt2d.grid.dim1)
        inp.nh = length(eqt2d.grid.dim2)

        inp.psirz = eqt2d.psi
        inp.npsi1d = length(eqt.profiles_1d.psi)
        inp.qpsi = abs.(eqt.profiles_1d.q) # RABBIT assumes positive q 
        inp.fpol = eqt.profiles_1d.f
        inp.sibry = eqt.global_quantities.psi_boundary 
        inp.simag = eqt.global_quantities.psi_axis
        inp.signip = sign(eqt.global_quantities.ip)
        inp.rmaxis = eqt.global_quantities.magnetic_axis.r
        inp.zmaxis = eqt.global_quantities.magnetic_axis.z 

        inp.r = eqt2d.grid.dim1
        inp.z = eqt2d.grid.dim2

        inp.rhorz = eqt2d.phi ./ (pi .* eq.vacuum_toroidal_field.b0[i])
        inp.psi = eqt.profiles_1d.psi
        inp.vol = eqt.profiles_1d.volume
        inp.area = eqt.profiles_1d.area

        cp1d = dd.core_profiles.profiles_1d[time]

        inp.rho = range(0.0, stop=1.0, length=101)
        inp.n_rho = length(inp.rho)
        inp.te = IMAS.interp1d(cp1d.grid.rho_tor_norm,cp1d.electrons.temperature).(inp.rho) .* eV_to_keV
        inp.dene = IMAS.interp1d(cp1d.grid.rho_tor_norm,cp1d.electrons.density).(inp.rho) .* cm3_to_m3
        inp.rot_freq_tor = inp.rho .* 0.0
        inp.zeff = IMAS.interp1d(cp1d.grid.rho_tor_norm,cp1d.zeff).(inp.rho)
    
        for i in 2:length(cp1d.ion)
            @assert cp1d.ion[1].temperature == cp1d.ion[i].temperature "All ion temperatures should be the same"
        end
        inp.ti = IMAS.interp1d(cp1d.grid.rho_tor_norm,cp1d.ion[1].temperature).(inp.rho) .* eV_to_keV

        pnbis = []
        for i in 1:length(dd.nbi.unit)
            for j in 1:length(dd.nbi.unit[1].power_launched.data)
                push!(pnbis, dd.nbi.unit[i].power_launched.data[j]) # take this from pulse schedule instead 
            end
        end
        inp.pnbi = pnbis


        push!(all_inputs, inp)
    end

    if length(all_inputs) == 1
        inp = deepcopy(all_inputs[1])
        inp.time -= 1E6
        push!(all_inputs, inp)
    end

    return all_inputs

end

function write_timetraces(all_inputs::Vector{RABBITinput})
    nw = 5
    
    open("timetraces.dat", "w") do io 
        println(io, "         ", length(all_inputs))
        println(io, "         ", all_inputs[1].n_rho)
        println(io, "rho_tor")
        print(io, cropdata_f([all_inputs[i].time for i in eachindex(all_inputs)], nw))
        print(io, cropdata_f(all_inputs[1].rho, nw))
        print(io, cropdata_f([all_inputs[i].te for i in eachindex(all_inputs)], nw))
        print(io, cropdata_f([all_inputs[i].ti for i in eachindex(all_inputs)], nw))
        print(io, cropdata_e([all_inputs[i].dene for i in eachindex(all_inputs)], nw))
        print(io, cropdata_f([all_inputs[i].rot_freq_tor for i in eachindex(all_inputs)], nw))
        print(io, cropdata_f([all_inputs[i].zeff for i in eachindex(all_inputs)], nw))
        print(io, cropdata_f(all_inputs[1].pnbi, nw))
    end
end

function write_equilibria(input::RABBITinput, filename::AbstractString)
    open(filename, "w") do io 
        println(io, "          ", input.nw)
        println(io, "          ", input.nh)
        print(io, print6(input.r))
        print(io, print6(input.z))

        for i in range(1,length(input.psirz), step = input.nw)
            print(io, print6(input.psirz[i:i+(input.nw - 1)]))
        end

        for i in range(1,length(input.rhorz), step = input.nw)
            print(io, print6(input.rhorz[i:i+(input.nw - 1)]))
        end
        println(io, "          ", input.npsi1d)
        print(io, print6(input.psi))
        print(io, print6(input.vol)) 
        print(io, print6(input.area))
        print(io, print6(input.rho))
        print(io, print6(input.qpsi)) 
        print(io, print6(input.fpol))
        print(io, @sprintf("%12.6f%12.6f%12.6f%12.6f%12.6f\n", input.sibry, input.simag, input.signip, input.rmaxis, input.zmaxis))

    end
end

function write_equilibria(all_inputs::Vector{RABBITinput})
    mkdir("equ")
    for i in eachindex(all_inputs)
        filename = "equ/equ_$i.dat"
        write_equilibria(all_inputs[i], filename)
    end
end