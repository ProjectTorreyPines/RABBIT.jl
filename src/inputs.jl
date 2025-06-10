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
    pnbi::Union{Vector{Vector{Float64}},Missing} = missing

    # equilibria 
    nw::Union{Int,Missing} = missing
    nh::Union{Int,Missing} = missing
    psirz::Union{Array{Float64},Missing} = missing
    npsi1d::Union{Int,Missing} = missing
    qpsi::Union{Vector{Float64},Missing} = missing
    fpol::Union{Vector{Float64},Missing} = missing
    sibry::Union{Float64,Missing} = missing
    simag::Union{Float64,Missing} = missing
    signip::Union{Float64,Missing} = missing
    rmaxis::Union{Float64,Missing} = missing
    zmaxis::Union{Float64,Missing} = missing
    eq_rho::Union{Vector{Float64},Missing} = missing

    # AuxQuantities
    r::Union{Vector{Float64},Missing} = missing
    z::Union{Vector{Float64},Missing} = missing
    rhorz::Union{Array{Float64},Missing} = missing
    psi::Union{Vector{Float64},Missing} = missing
    vol::Union{Vector{Float64},Missing} = missing
    area::Union{Vector{Float64},Missing} = missing

    # beams 
    n_sources::Union{Int,Missing} = missing
    nv::Union{Int,Missing} = missing
    start_pos::Union{Array{Float64},Missing} = missing
    beam_unit_vector::Union{Array{Float64},Missing} = missing
    beam_width_polynomial_coefficients::Union{Array{Float64},Missing} = missing
    injection_energy::Union{Array{Float64},Missing} = missing
    particle_fraction::Union{Array{Float64},Missing} = missing # particle fraction of full/half/third energy
    a_beam::Union{Array{Float64},Missing} = missing

end

Base.@kwdef mutable struct Timetraces
    te::Union{Vector{Vector{Float64}},Missing}
    ti::Union{Vector{Vector{Float64}},Missing}
    dene::Union{Vector{Vector{Float64}},Missing}
    rot_freq_tor::Union{Vector{Vector{Float64}},Missing}
    zeff::Union{Vector{Vector{Float64}},Missing}
end

function extract_timetraces(all_inputs, varname)
    n_timeslices = length(all_inputs)
    n_rho = length(all_inputs[1].te)
    return [[getfield(all_inputs[itime], varname)[irho] for itime in 1:n_timeslices] for irho in 1:n_rho]
end


function write_timetraces(all_inputs::Vector{RABBITinput})
    timetraces = Timetraces(
        extract_timetraces(all_inputs, :te),
        extract_timetraces(all_inputs, :ti),
        extract_timetraces(all_inputs, :dene),
        extract_timetraces(all_inputs, :rot_freq_tor),
        extract_timetraces(all_inputs, :zeff))

    nw = 5
    open("timetraces.dat", "w") do io
        println(io, "         ", length(all_inputs))
        println(io, "         ", all_inputs[1].n_rho)
        println(io, "rho_tor")
        print(io, cropdata_f([all_inputs[i].time for i in eachindex(all_inputs)], nw))
        print(io, cropdata_f(all_inputs[1].rho, nw))
        print(io, cropdata_f([timetraces.te[i] for i in 1:length(timetraces.te)], nw))
        print(io, cropdata_f([timetraces.ti[i] for i in 1:length(timetraces.ti)], nw))
        print(io, cropdata_e([timetraces.dene[i] for i in 1:length(timetraces.dene)], nw))
        print(io, cropdata_f([timetraces.rot_freq_tor[i] for i in 1:length(timetraces.rot_freq_tor)], nw))
        print(io, cropdata_f([timetraces.zeff[i] for i in 1:length(timetraces.zeff)], nw))
        return print(io, cropdata_f(all_inputs[1].pnbi, nw))
    end
end

function write_equilibria(input::RABBITinput, filename::AbstractString)
    open(filename, "w") do io
        println(io, "          ", input.nw)
        println(io, "          ", input.nh)
        print(io, print6(input.r))
        print(io, print6(input.z))

        for i in range(1, length(input.psirz); step=input.nw)
            print(io, print6(input.psirz[i:i+(input.nw-1)]))
        end

        for i in range(1, length(input.rhorz); step=input.nw)
            print(io, print6(input.rhorz[i:i+(input.nw-1)]))
        end
        println(io, "          ", input.npsi1d)
        print(io, print6(input.psi))
        print(io, print6(input.vol))
        print(io, print6(input.area))
        print(io, print6(input.eq_rho))
        print(io, print6(input.qpsi))
        print(io, print6(input.fpol))
        return print(io, @sprintf("%12.6f%12.6f%12.6f%12.6f%12.6f\n", input.sibry, input.simag, input.signip, input.rmaxis, input.zmaxis))

    end
end

function write_equilibria(all_inputs::Vector{RABBITinput})
    mkdir("equ")
    for i in eachindex(all_inputs)
        filename = "equ/equ_$i.dat"
        write_equilibria(all_inputs[i], filename)
    end
end

function write_options(vessel_hfs::Float64, vessel_lfs::Float64)
    table_path = abspath(joinpath(dirname(@__DIR__), "tables_highRes"))
    open("options.nml", "w") do io
        println(io, "&species
 Aimp=12.00 ! mass of impurity species
 Zimp=6.00 ! charge of impurity species
 /
 &output_settings
 it_orbout(:) = 71,-1, 115,115,115,115, -1,-1,-1,-1, -1, 200,200, -1
 nrhoout = 20 ! length of output rho grid
 writedistfun = .false. ! toggle writing of distribution function
 writedepo2d = .false. ! toggle writing of 2D deposition profile
 /
 &physics
 jumpcor=2 ! smoothed number of gridpoints in plasma center
 flast=10. ! fudge factor in colrad model to increase loss of highest state
 Rmax = ", vessel_lfs)
        println(io, " ! start of NBI rays (m) -> set similar to vacuum vessel
  Rmin = ", vessel_hfs)
        println(io, " ! end of NBI rays (m) -> set similar to vacuum vessel
  Rlim = 2.26 ! flux surface passing through Rlim and zlim is considered as limiter
  zlim = -0.05
  torqjxb_model = 3 ! 0 = no orbit averaging (oav) deposition, 1 = calc jxb, 2 = oav deposition, 3 = calc jxb, but rescale it to match deposited torque
  beamlosswall=.false. ! option to simulate partial beam scraping on the vessel wall - switch off for DIII-D
  table_path= '", table_path, "'")
        return println(
            io,
            "/
 &numerics
 distfun_nv=20 ! number of velocity divisions in distribution function; if 0, no output of distribution function + neutron rate
 distfun_vmax=3.3e6 ! maximum velocity in distribution function (m/s)
 norbits=20 ! number of  orbits (per beam energy component)
 rescalecor=.false.
 vchangecor=.false.
 /
 &fpsol
 /


        "
        )
    end
end

function write_beams(all_inputs::Vector{RABBITinput})
    open("beams.dat", "w") do io
        println(io, "# no. of sources:")
        println(io, "        ", all_inputs[1].n_sources)
        println(io, "# nv:")
        println(io, "        ", all_inputs[1].nv)
        println(io, "# start_pos: [m]")
        print(io, cropdata_f(all_inputs[1].start_pos, 3))
        println(io, "# beam unit vector:")
        print(io, cropdata_f(all_inputs[1].beam_unit_vector, 3))
        println(io, "# beam-width-polynomial coefficients:")
        print(io, cropdata_f(all_inputs[1].beam_width_polynomial_coefficients, 3))
        println(io, "# Injection energy [eV]:")
        print(io, cropdata_f(all_inputs[1].injection_energy, 5))
        println(io, "# Particle fraction of full/half/third energy:")
        print(io, cropdata_f(all_inputs[1].particle_fraction, 3))
        println(io, "# A beam [u]")
        return print(io, cropdata_f(all_inputs[1].a_beam, 5))
    end

end

"""
    run_RABBIT(all_inputs::Vector{RABBITinput}; remove_inputs::Bool=true)

Writes RABBIT input files (equ_X.dat, timetraces.dat, beams.dat, options.nml) to a run directory and executes RABBIT on that directory.

Set remove_inputs=false to keep run directory containing full input and output files.
"""
function run_RABBIT(all_inputs::Vector{RABBITinput}, vessel_hfs::Float64, vessel_lfs::Float64; remove_inputs::Bool=true)
    folder = mktempdir()

    exec_path = abspath(joinpath(dirname(@__DIR__), "rabbit"))

    old_dir = pwd()
    try
        cd(folder)

        write_equilibria(all_inputs)

        write_timetraces(all_inputs)

        write_beams(all_inputs)

        write_options(vessel_hfs, vessel_lfs)

    finally
        cd("../")
    end

    open("command.sh", "w") do io
        return write(io, string(exec_path), " $filename &> command.log")
    end

    output = try
        run(Cmd(`bash command.sh`))
        read_outputs(pwd(); filename)
    catch e
        txt = open("command.log", "r") do io
            return split(read(io, String), "\n")
        end
        @error "Error running RABBIT" * join(txt[max(1, length(txt) - 50):end], "\n")
        rethrow(e)
    finally
        cd(old_dir)
    end

    if remove_inputs
        rm(folder; force=true, recursive=true)
    end

    return output
end