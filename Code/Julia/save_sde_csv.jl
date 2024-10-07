### Save Out Models to CSV  ###
using RData
using Interpolations
using DifferentialEquations
using ProgressBars
using Symbolics

################################################
## Extracts the simulated communities from SDE solution for user defined period, and adds supporting info prior to saving
################################################
function ts_harvest_save(solution::Union{RODESolution,ODESolution,EnsembleSolution},tspan::Tuple)
    out_dict = Dict()

if isa(solution,EnsembleSolution)
    for i in eachindex(solution)
        ens_dict = Dict()
        temp_sol = solution[i]

        if length(unique(floor.(solution[i].t))) == length(tspan[1]:tspan[2]) 
        
        ens_dict[:sol] = transpose(convert(Array,VectorOfArray(temp_sol.u[indexin(unique(floor.(solution[i].t)),floor.(solution[i].t))]))) #convert from array of vectors to array
        #due to some duplication of timepoints, only keep unique timepoints

        #ens_dict[:sol] = temp_sol.u
        #ens_dict[:sol] = temp_sol
        ens_dict[:par] = Dict(:A => temp_sol.prob.p[:A], 
        :r => temp_sol.prob.p[:r], 
        :spp => temp_sol.prob.p[:spp],
        :ctrl => temp_sol.prob.p[:signal][end-1]) #extract community info

        #ens_dict[:par] = temp_sol.prob.p
        ens_dict[:u0] = temp_sol.prob.u0 
        out_dict[string(i)] = ens_dict
        else
            ens_dict[:sol] = nothing  

            ens_dict[:par] = Dict(:A => temp_sol.prob.p[:A], 
            :r => temp_sol.prob.p[:r], 
            :spp => temp_sol.prob.p[:spp],
            :ctrl => temp_sol.prob.p[:signal][end-1])

            ens_dict[:u0] = temp_sol.prob.u0
            out_dict[string(i)] = ens_dict
        end
    end
else

    if length(unique(floor.(solution.t))) == length(tspan[1]:tspan[2]) 

    out_dict[:sol] = transpose(convert(Array,VectorOfArray(solution.u[indexin(unique(floor.(solution.t)),floor.(solution.t))])))
    #out_dict[:sol] = solution.u
    out_dict[:par] = Dict(:A => solution.prob.p[:A], 
    :r => solution.prob.p[:r], 
    :spp => solution.prob.p[:spp],
    :ctrl => solution.prob.p[:signal][end-1])

    #out_dict[:par] = solution.prob.p
    out_dict[:u0] = solution.prob.u0
    else
        out_dict[:sol] = nothing   
        out_dict[:par] = Dict(:A => solution.prob.p[:A], 
    :r => solution.prob.p[:r], 
    :spp => solution.prob.p[:spp],
    :ctrl => solution.prob.p[:signal][end-1])

    out_dict[:u0] = solution.prob.u0
    end
end
return sort(out_dict,by = x -> parse(Int, x)) #order dictionary by number not alphabetically
end

# motif = motif_dat
# spp = [3]
# tspan = (0,500)
# signal = signal_motif
# sims = 100
# saveat=1.0
# noise_col = "white"
# solution = stoc_gen_sol

# function motif_sim(motif::Vector{DictoVec{Any}}, spp::Union{Vector, Matrix}, noise::Union{Vector, Matrix}, tspan::Tuple, signal::Interpolations.Extrapolation, sims::Int64; noise_col = "white",saveat::Float64 = 0.1)

#     motif_dict = Dict()

#     for i in ProgressBar(1:length(motif))

#     parameters = Dict(:A => motif[i]["A_matrix"], :r => motif[i]["R"],:ex_ctrl => repeat([1.0],length(noise)), :spp => spp,:W => noise,:signal => signal)

#     extinction_event = VectorContinuousCallback(extinction!,extinction_affect!,length(noise))
#     du0 = copy(motif[i]["Neq"])
#     jac_sparsity = Symbolics.jacobian_sparsity((du,u)->gen_lotka_volterra_harvest!(du,u,parameters,0.0),du0,motif[i]["Neq"])
#     f = ODEFunction(gen_lotka_volterra_harvest!;jac_prototype=float.(jac_sparsity))

#     if noise_col == "white"
#             stoc_gen_prob = SDEProblem(f, stoc_gen_lv!,motif[i]["Neq"], tspan, parameters,noise = WienerProcess(0.0,repeat([0.0],length(noise))),dt=0.1,callback = extinction_event,progress = true)
#         end
#         if noise_col == "coloured"
#             stoc_gen_prob = SDEProblem(f, stoc_gen_lv!,motif[i]["Neq"], tspan, parameters,noise = OrnsteinUhlenbeckProcess(0.5,0.0,0.5,0.0,repeat([0.0],length(noise))),dt=0.1,callback = extinction_event,progress = true)
#         end
    
#     ens_prob = EnsembleProblem(stoc_gen_prob)
#     stoc_gen_sol = solve(ens_prob,EM(),EnsembleThreads(),saveat=saveat,trajectories=sims)
#     #out = ts_save(stoc_gen_sol,tspan)
#     #motif_dict[string(i)] = out
#     motif_dict[string(i)] = ts_harvest_save(stoc_gen_sol,tspan)

#     end

# return sort(motif_dict,by = x -> parse(Int, x))
# end

################################################
## Further tidying of simualtions prior to save out
################################################
function extract_ts_for_csv(container::AbstractDict, label::String, trim::UnitRange,round_digits::Int64)
    
    #trim = (trim[1] + 1):(trim[end] + 3) #convert from time point to index position (i.e t0 = index 1)
    trim = trim .+ 1 #convert from time point to index position (i.e t0 = index 1)

    out_df = DataFrames.DataFrame(Tables.table(container["1"]["1"][:sol][trim,:])) # if first iteration, create Dataframe
    transform!(out_df, names(out_df) .=> ByRow(x -> round(x,digits=round_digits)) .=>  names(out_df)) #round values to digits specified in argument 
    out_df.max_stress .= container["1"]["1"][:par][:ctrl]
    #out_df.target_spp .= container["1"]["1"][:par][:spp]
    out_df.target_spp .= join(container["1"]["1"][:par][:spp],"_")

    out_df.community .= "1"
    out_df.sim .= "1"

    for i in keys(container) #loop through top layers
        for j in keys(container[i]) #loop through nested layers
            if i != "1" || j != "1" #if not first iteration, append successive layers to above Dataframe
                temp = DataFrames.DataFrame(Tables.table(container[i][j][:sol][trim,:]))
                transform!(temp, names(temp) .=> ByRow(x -> round(x,digits=round_digits)) .=>  names(temp))
                temp.max_stress .= container[i][j][:par][:ctrl]
                #temp.target_spp .= container[i][j][:par][:spp]
                temp.target_spp .= join(container[i][j][:par][:spp],"_")

                temp.community .= i
                temp.sim .= j
                append!(out_df,temp)
            end
        end
    end
    out_df.motif .= label
    return out_df
end

################################################
## Simulates the harvest form of GLV ##
################################################
function motif_sim(motif::Vector{DictoVec{Any}}, u0,spp::Union{Vector, Matrix,String}, noise::Union{Vector, Matrix}, tspan::Tuple, signal::Interpolations.Extrapolation, sims::Int64; noise_col = "white",saveat::Float64 = 0.1)

    motif_dict = Dict()

    for i in ProgressBar(1:length(motif)) #loop across all communities

    parameters = Dict(:A => motif[i]["A_matrix"], :r => motif[i]["R"],:ex_ctrl => repeat([1.0],length(noise)), :spp => Vector{Int64}[],:W => noise,:signal => signal)
        # prepare community data
    if spp == "auto" #automatically define which species to stress
        #push!(parameters,:spp => [findmax(motif[i]["Neq"])[2];])
        push!(parameters,:spp => [findmax(x -> x == motif[i]["Neq"][motif[i]["tlvl"][:,1] .== 3.0][1], motif[i]["Neq"])[2];])
    else
        push!(parameters,:spp => spp) #or accept user provide species(s)
    end

    parameters = (; parameters...) #convert Dict to NamedTuple

    extinction_event = VectorContinuousCallback(extinction!,extinction_affect!,length(noise)) #define extinction function depending on number of species
    
    if ismissing(u0) #if initial abundances not explicitly provided, copy from initial data
        du0 = copy(motif[i]["Neq"])
    else
        du0 = copy(u0)
    end

    jac_sparsity = Symbolics.jacobian_sparsity((du,u)->gen_lotka_volterra_harvest!(du,u,parameters,0.0),du0,motif[i]["Neq"]) #speed up solver by estimating Jacobian sparsity
    f = ODEFunction(gen_lotka_volterra_harvest!;jac_prototype=float.(jac_sparsity)) #define ODE problem

    if noise_col == "white" #simulate using multiplicative white noise
            stoc_gen_prob = SDEProblem(f, stoc_gen_lv!,motif[i]["Neq"], tspan, parameters,noise = WienerProcess(0.0,repeat([0.0],length(noise))),dt=0.1,callback = extinction_event,progress = true)
        end
        if noise_col == "coloured" #simulate using multiplicative correlated noise (autocorrelation = 0.5)
            stoc_gen_prob = SDEProblem(f, stoc_gen_lv!,motif[i]["Neq"], tspan, parameters,noise = OrnsteinUhlenbeckProcess(0.5,0.0,0.5,0.0,repeat([0.0],length(noise))),dt=0.1,callback = extinction_event,progress = true)
        end
    
    ens_prob = EnsembleProblem(stoc_gen_prob) #indicate multiple simulations of a community are to be run
    stoc_gen_sol = solve(ens_prob,EM(),EnsembleThreads(),saveat=saveat,trajectories=sims) #simulate community multiple times

    motif_dict[string(i)] = ts_harvest_save(stoc_gen_sol,tspan) #extract solution info

    end

return sort(motif_dict,by = x -> parse(Int, x)) #order by simulation name
end
################################################################################################
############ Dakos Invasive functions ########################
################################################################################################
function ts_invasive_save(solution::Union{RODESolution,ODESolution,EnsembleSolution},tspan::Tuple)
    out_dict = Dict()

if isa(solution,EnsembleSolution)
    for i in eachindex(solution)
        ens_dict = Dict()
        temp_sol = solution[i]

        if length(unique(floor.(solution[i].t))) == length(tspan[1]:tspan[2]) 

        ens_dict[:sol] = transpose(convert(Array,VectorOfArray(temp_sol.u[indexin(unique(floor.(solution[i].t)),floor.(solution[i].t))]))) #convert from array of vectors to array
        #due to some duplication of timepoints, only keep unique timepoints

        #ens_dict[:sol] = temp_sol.u
        #ens_dict[:sol] = temp_sol
        ens_dict[:par] = Dict(:A => temp_sol.prob.p[:A], 
        :r => temp_sol.prob.p[:r], 
        :k => temp_sol.prob.p[:k], 
        :eta => temp_sol.prob.p[:eta], 
        :spp => temp_sol.prob.p[:spp],
        :ctrl => temp_sol.prob.p[:signal][end-1])

        #ens_dict[:par] = temp_sol.prob.p
        ens_dict[:u0] = temp_sol.prob.u0
        out_dict[string(i)] = ens_dict
        else
            ens_dict[:sol] = nothing  

            ens_dict[:par] = Dict(:A => temp_sol.prob.p[:A], 
            :r => temp_sol.prob.p[:r], 
            :k => temp_sol.prob.p[:k], 
            :eta => temp_sol.prob.p[:eta], 
            :spp => temp_sol.prob.p[:spp],
            :ctrl => temp_sol.prob.p[:signal][end-1])

            ens_dict[:u0] = temp_sol.prob.u0
            out_dict[string(i)] = ens_dict
        end
    end
else

    if length(unique(floor.(solution.t))) == length(tspan[1]:tspan[2]) 

    out_dict[:sol] = transpose(convert(Array,VectorOfArray(solution.u[indexin(unique(floor.(solution.t)),floor.(solution.t))])))
    #out_dict[:sol] = solution.u
    out_dict[:par] = Dict(:A => solution.prob.p[:A], 
    :r => solution.prob.p[:r], 
    :k => solution.prob.p[:k], 
    :eta => solution.prob.p[:eta], 
    :spp => solution.prob.p[:spp],
    :ctrl => solution.prob.p[:signal][end-1])

    #out_dict[:par] = solution.prob.p
    out_dict[:u0] = solution.prob.u0
    else
        out_dict[:sol] = nothing   
        out_dict[:par] =  Dict(:A => solution.prob.p[:A], 
        :r => solution.prob.p[:r], 
        :k => solution.prob.p[:k], 
        :eta => solution.prob.p[:eta], 
        :spp => solution.prob.p[:spp],
        :ctrl => solution.prob.p[:signal][end-1])

    out_dict[:u0] = solution.prob.u0
    end
end
return sort(out_dict,by = x -> parse(Int, x)) #order dictionary by number not alphabetically
end

# function invasive_sim(motif::Vector{DictoVec{}}, spp::Union{Vector, Matrix}, noise,tspan::Tuple, signal::Interpolations.Extrapolation, sims::Int64; noise_col = "white",saveat::Float64 = 0.1)

#     motif_dict = Dict()

#     for i in ProgressBar(1:length(motif))

#     parameters = Dict(:A => motif[i]["A_matrix"], 
#     :r => motif[i]["R"],
#     :k => motif[i]["K"],
#     :eta => motif[i]["eta"],
#     :imm => motif[i]["imm"],
#     :spp => spp,
#     #:W => [0.1,0.1,0.1,0.1],
#     :W => noise,
#     :ex_ctrl => repeat([1.0],length(noise)),
#     :signal => signal)

#     extinction_event = VectorContinuousCallback(extinction!,extinction_affect!,length(noise))
    
#     du0 = copy(motif[i]["Neq"])
#     jac_sparsity = Symbolics.jacobian_sparsity((du,u)->gen_lotka_volterra_invasive_harvest!(du,u,parameters,0.0),du0,motif[i]["Neq"])
#     f = ODEFunction(gen_lotka_volterra_invasive_harvest!;jac_prototype=float.(jac_sparsity))

#     if noise_col == "white"
#         stoc_gen_prob = SDEProblem(f, stoc_gen_lv!,motif[i]["Neq"], tspan, parameters,noise = WienerProcess(0.0,repeat([0.0],length(noise))),dt=0.1,callback = extinction_event,progress = true)
#     end
#     if noise_col == "coloured"
#         stoc_gen_prob = SDEProblem(f, stoc_gen_lv!,motif[i]["Neq"], tspan, parameters,noise = OrnsteinUhlenbeckProcess(0.5,0.0,0.5,0.0,repeat([0.0],length(noise))),dt=0.1,callback = extinction_event,progress = true)
#     end

#     ens_prob = EnsembleProblem(stoc_gen_prob)
#     stoc_gen_sol = solve(ens_prob,EM(),EnsembleThreads(),saveat=saveat,trajectories=sims)
#     #out = ts_save(stoc_gen_sol,tspan)
#     #motif_dict[string(i)] = out
#     motif_dict[string(i)] = ts_invasive_save(stoc_gen_sol,tspan)

#     end

# return sort(motif_dict,by = x -> parse(Int, x))
# end

function extract_ts_for_csv_invasive(container::AbstractDict, label::String, trim::UnitRange,round_digits::Int64)
    
    #trim = (trim[1] + 1):(trim[end] + 3) #convert from time point to index position (i.e t0 = index 1)
    trim = trim .+ 1 #convert from time point to index position (i.e t0 = index 1)

    out_df = DataFrames.DataFrame(Tables.table(container["1"][:sol][trim,:])) # if first iteration, create Dataframe
    transform!(out_df, names(out_df) .=> ByRow(x -> round(x,digits=round_digits)) .=>  names(out_df)) #round values to digits specified in argument 
    out_df.max_stress .= container["1"][:par][:ctrl]
    out_df.community .= label
    out_df.sim .= "1"

    for i in keys(container) #loop through top layers
            if i != "1"  #if not first iteration, append successive layers to above Dataframe
                temp = DataFrames.DataFrame(Tables.table(container[i][:sol][trim,:]))
                transform!(temp, names(temp) .=> ByRow(x -> round(x,digits=round_digits)) .=>  names(temp))
                temp.max_stress .= container[i][:par][:ctrl]
                temp.community .= label
                temp.sim .= i
                append!(out_df,temp)
        end
    end
    return out_df
end

function invasive_sim(motif::Vector{DictoVec{}}, u0,spp::Union{Vector, Matrix,String}, noise,tspan::Tuple, signal::Interpolations.Extrapolation, sims::Int64; noise_col = "white",saveat::Float64 = 0.1)

    motif_dict = Dict()

    for i in ProgressBar(1:length(motif))

    parameters = Dict(:A => motif[i]["A_matrix"], 
    :r => motif[i]["R"],
    :k => motif[i]["K"],
    :eta => motif[i]["eta"],
    :imm => motif[i]["imm"],
    :spp => Vector{Int64}[],
    #:W => [0.1,0.1,0.1,0.1],
    :W => noise,
    :ex_ctrl => repeat([1.0],length(noise)),
    :signal => signal)

#=     if spp == "auto"
        #push!(parameters,:spp => [findmin(motif[i]["Neq"])[2];])
        push!(parameters,:spp => [findmax(x -> x == motif[i]["R"][motif[i]["tlvl"][:,1] .== 3.0][1], motif[i]["R"])[2];])

    else
        push!(parameters,:spp => spp)
    end =#
    if spp == "auto"
        #push!(parameters,:spp => [findmin(motif[i]["Neq"])[2];])

        #if sum([motif[i]["tlvl"][:,1] .== 3.0][1]) == 0|1
        #push!(parameters,:spp => [findmax(x -> x == motif[i]["R"][motif[i]["tlvl"][:,1] .== 2.0][1], motif[i]["R"])[2];])
        #else
        push!(parameters,:spp => [findmax(x -> x == motif[i]["R"][motif[i]["tlvl"][:,1] .== 3.0][1], motif[i]["R"])[2];])
        #end
    elseif spp == "spp"
        push!(parameters,:spp => motif[i]["spp"])
    else 
        push!(parameters,:spp => spp)
    end
   
    parameters = (; parameters...) #convert Dict to NamedTuple

    extinction_event = VectorContinuousCallback(extinction!,extinction_affect!,length(noise))
    
    if ismissing(u0)
        du0 = copy(motif[i]["Neq"])
    else
        du0 = copy(u0)
    end

    jac_sparsity = Symbolics.jacobian_sparsity((du,u)->gen_lotka_volterra_invasive_harvest!(du,u,parameters,0.0),du0,motif[i]["Neq"])
    f = ODEFunction(gen_lotka_volterra_invasive_harvest!;jac_prototype=float.(jac_sparsity))

    if noise_col == "white"
        stoc_gen_prob = SDEProblem(f, stoc_gen_lv!,du0, tspan, parameters,noise = WienerProcess(0.0,repeat([0.0],length(noise))),dt=0.1,callback = extinction_event,progress = true)
    end
    if noise_col == "coloured"
        stoc_gen_prob = SDEProblem(f, stoc_gen_lv!,du0, tspan, parameters,noise = OrnsteinUhlenbeckProcess(0.5,0.0,0.5,0.0,repeat([0.0],length(noise))),dt=0.1,callback = extinction_event,progress = true)
    end

    ens_prob = EnsembleProblem(stoc_gen_prob)
    stoc_gen_sol = solve(ens_prob,EM(),EnsembleThreads(),saveat=saveat,trajectories=sims)

    #out = ts_save(stoc_gen_sol,tspan)
    #motif_dict[string(i)] = out
    motif_dict[string(i)] = ts_invasive_save(stoc_gen_sol,tspan)

    end

return sort(motif_dict,by = x -> parse(Int, x))
end
