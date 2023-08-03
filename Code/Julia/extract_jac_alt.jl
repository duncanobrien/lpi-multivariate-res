## Alternative extraction Jacobian from ODEProblem ##
#PREFERRED
#Numerically estimates the Jacobian of the final community using forward differencing

using DifferentialEquations
using LinearAlgebra
using ForwardDiff
#using SparsityDetection
#using FiniteDiff

function critical_sys_harvest(control_range::StepRangeLen,odeproblem::ODEProblem,control_p::Dict;abs_max = false)

    jac_vector = hcat(control_range,zeros(length(control_range))) #create array to be populated

    temp_parameters = control_p

    for j in 1:(size(jac_vector)[1])
        #setindex!(temp_parameters,jac_vector[j,1] , :signal) #change stress parameter
        temp_parameters = (temp_parameters..., signal=jac_vector[j,1])

        prob_temp = remake(odeproblem; p=temp_parameters)
        #prob_temp = remake(odeproblem; p= NamedTuple([pair for pair in temp_parameters])) #alter ODE with new parameter values
        sol_temp = solve(prob_temp)

        jac = ForwardDiff.jacobian((du,u)->gen_lotka_volterra_harvest_remake!(du,u,temp_parameters,0.0),sol_temp.u[size(sol_temp.u)[1]], sol_temp.u[size(sol_temp.u)[1]])

         #jac = zeros(size(sol_temp.u[1])[1],size(sol_temp.u[1])[1])
        #FiniteDiff.finite_difference_jacobian!(jac,(du,u)->gen_lotka_volterra_invasive_remake!(du,u,temp_parameters,0.2),sol_temp.u[size(sol_temp.u)[1]])

        #int = init(prob_temp,Rosenbrock23(autodiff=true,concrete_jac = concrete_jac))
        #try  jac = ForwardDiff.jacobian((du,u)->gen_lotka_volterra_harvest_remake!(du,u,temp_parameters,0.0),sol_temp.u[size(sol_temp.u)[1]], sol_temp.u[size(sol_temp.u)[1]]) 
        
        #if sol_temp.retcode != :Default #if integrator fails, set value to -Inf
        #    jac_vector[j,2] = -Inf
            #break
        #else
            temp = real(LinearAlgebra.eigvals(jac)) #extract eigenvalues of Jacobian
            if abs_max == true
                jac_vector[j,2] = temp[abs.(temp) .== maximum(abs.(temp))][1] #return absolute maximum eigenvalue
            else
                jac_vector[j,2] = temp[temp .== maximum(temp)][1] #return only maximum eigenvalue
            end
        #end
        #catch
        #    jac_vector[j,2] = -Inf
        #end

    end
    return jac_vector
end

function critical_sys_invasive(control_range::StepRangeLen,odeproblem::ODEProblem,control_p::Dict;abs_max = false)

    jac_vector = hcat(control_range,zeros(length(control_range))) #create array to be populated

    temp_parameters = control_p

    for j in 1:(size(jac_vector)[1])
        setindex!(temp_parameters,jac_vector[j,1] , :signal) #change stress parameter
        temp_parameters = (temp_parameters..., signal=jac_vector[j,1])

        prob_temp = remake(odeproblem; p=temp_parameters)
        #prob_temp = remake(odeproblem; p= NamedTuple([pair for pair in temp_parameters])) #alter ODE with new parameter values
        sol_temp = solve(prob_temp)

        jac = ForwardDiff.jacobian((du,u)->gen_lotka_volterra_invasive_remake!(du,u,temp_parameters,0.0),sol_temp.u[size(sol_temp.u)[1]], sol_temp.u[size(sol_temp.u)[1]])

        #jac = zeros(size(sol_temp.u[1])[1],size(sol_temp.u[1])[1])
        #FiniteDiff.finite_difference_jacobian!(jac,(du,u)->gen_lotka_volterra_invasive_remake!(du,u,temp_parameters,0.2),sol_temp.u[size(sol_temp.u)[1]])

        #int = init(prob_temp,Rosenbrock23(autodiff=true,concrete_jac = concrete_jac))
        #try  jac = ForwardDiff.jacobian((du,u)->gen_lotka_volterra_harvest_remake!(du,u,temp_parameters,0.0),sol_temp.u[size(sol_temp.u)[1]], sol_temp.u[size(sol_temp.u)[1]]) 
        
        #if sol_temp.retcode != :Default #if integrator fails, set value to -Inf
        #    jac_vector[j,2] = -Inf
            #break
        #else
            temp = real(LinearAlgebra.eigvals(jac)) #extract eigenvalues of Jacobian
            if abs_max == true
                jac_vector[j,2] = temp[abs.(temp) .== maximum(abs.(temp))][1] #return absolute maximum eigenvalue
            else
                jac_vector[j,2] = temp[temp .== maximum(temp)][1] #return only maximum eigenvalue
            end
        #end
        #catch
        #    jac_vector[j,2] = -Inf
        #end

    end
    return jac_vector
end

function extract_max_eigvec_alt(container,control_range::StepRangeLen, t::Int64, spp; model = "harvest",abs_max = false)

    jac_marker =  hcat(1:length(container),zeros(1:length(container),3))

    for j in ProgressBar(1:length(container))

    if model == "invasive"
        params_remake = Dict(:A => container[j]["A_matrix"], 
        :r => container[j]["R"],
        :k => container[j]["K"],
        :eta => container[j]["eta"],
        :imm => container[j]["imm"],
        :spp => Vector{Int64}[],
        :W => repeat([0.1],length(container[j]["R"])),
        :ex_ctrl =>repeat([1.0],length(container[j]["R"])),
        :signal => 0.2)
#= 
        if spp == "auto"
            #push!(params_remake,:spp => [findmin(container[j]["Neq"])[2];])
            push!(params_remake,:spp => [findmax(x -> x == container[j]["R"][container[j]["tlvl"][:,1] .== 3.0][1], container[j]["R"])[2];])

        else
            push!(params_remake,:spp => spp)
        end =#

        if spp == "auto"
            #if sum([container[j]["tlvl"][:,1] .== 3.0][1]) == 0|1
            #    push!(params_remake,:spp => [findmax(x -> x == container[j]["R"][container[j]["tlvl"][:,1] .== 2.0][1], container[j]["R"])[2];])
            #else
                push!(params_remake,:spp => [findmax(x -> x == container[j]["R"][container[j]["tlvl"][:,1] .== 3.0][1], container[j]["R"])[2];])
            #end
        else
            push!(params_remake,:spp => spp)
        end
        params_remake = (; params_remake...) #convert Dict to NamedTuple

        extinction_event = VectorContinuousCallback(extinction!,extinction_affect!,length(params_remake[:W]))
        gen_prob_remake = ODEProblem(gen_lotka_volterra_invasive_remake!, container[j]["Neq"], (0,t+200), params_remake,callback =  extinction_event)
        tmp = critical_sys_invasive(control_range,gen_prob_remake,params_remake,abs_max = abs_max)
    
        #if any(tmp[2:end,2] .> -0.15) && any(tmp[2:end,2] .> tmp[1,2]) #only collapsng if close to zero and stress parameter = 0 is not the largest value
        if any(tmp[2:end,2] .>= -0.25) && tmp[1,2] <= -0.25 #only collapsng if close to zero and 0 stress is not the largest value
            jac_marker[j,2] = true
        else
            jac_marker[j,2] = false
        end
        try 
            jac_marker[j,3] = tmp[tmp[:,2] .==  tmp[:,2][findfirst(tmp[:,2] .>= -0.25)],:][1]
            jac_marker[j,4] = tmp[:,2][findfirst(tmp[:,2] .>= -0.25)]
        catch 
            jac_marker[j,3] = NaN
            jac_marker[j,4] = NaN
        end
    end

    if model == "harvest"
        params_remake = Dict(:A => container[j]["A_matrix"], 
        :r => container[j]["R"],
        :spp => spp,
        :W => repeat([0.1],length(container[j]["R"])),
        :ex_ctrl => repeat([1.0],length(container[j]["R"])),
        :signal => 0.01)
 
        if spp == "auto"
            #push!(params_remake,:spp => [findmax(container[j]["Neq"])[2];])
            push!(params_remake,:spp => [findmax(x -> x == container[j]["Neq"][container[j]["tlvl"][:,1] .== 3.0][1], container[j]["Neq"])[2];])
            #push!(params_remake,:spp => [findmax(x -> x == container[j]["Neq"][container[j]["tlvl"][:,1] .== 4.0][1], container[j]["Neq"])[2];])

        else
            push!(params_remake,:spp => spp)
        end

        extinction_event = VectorContinuousCallback(extinction!,extinction_affect!,length(params_remake[:W]))
        gen_prob_remake = ODEProblem(gen_lotka_volterra_harvest_remake!, container[j]["Neq"], (0,t+200), params_remake,callback =  extinction_event)
        tmp = critical_sys_harvest(control_range,gen_prob_remake,params_remake,abs_max = abs_max)
        
        #if any(tmp[2:end,2] .>= -0.05) && any(tmp[2:end,2] .> tmp[1,2]) #only collapsng if close to zero and 0 stress is not the largest value
        if any(tmp[2:end,2] .>= -0.0001) && any(tmp[2:end,2] .> tmp[1,2]) #only collapsng if close to zero and 0 stress is not the largest value
            jac_marker[j,2] = true
        else
            jac_marker[j,2] = false
        end
        try 
            jac_marker[j,3] = tmp[tmp[:,2] .==  tmp[:,2][findfirst(tmp[:,2] .>= -0.001)],:][1]
            jac_marker[j,4] = tmp[:,2][findfirst(tmp[:,2] .>= -0.001)]
        catch 
            jac_marker[j,3] = NaN
            jac_marker[j,4] = NaN
        end
    end
    #jac_marker[j,3] = tmp[tmp[:,2] .==  maximum(tmp[:,2]),:][1]
    #jac_marker[j,4] = maximum(tmp[:,2])

    end
    jac_marker = DataFrames.DataFrame(Tables.table(jac_marker))
    rename!(jac_marker,[:community,:collapse,:stress_param,:max_eigval])
    transform!(jac_marker, :max_eigval .=> ByRow(x -> round(x,digits=5)) .=>  :max_eigval) #round values to digits specified in argument 
    transform!(jac_marker, :community => ByRow(x -> string(trunc(Int,x))) => :community)
    return jac_marker
end
