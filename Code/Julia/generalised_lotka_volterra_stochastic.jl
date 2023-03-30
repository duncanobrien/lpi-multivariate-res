#### Generalised Lotka Volterra - Stochastic ####
using DifferentialEquations
using Plots
using LinearAlgebra
using Interpolations #interpolate for harvesting/control parameters
using DataFrames #create dataframes
using CSV #export CSV files
using Distributions #uniform function

################################################
## Generic interaction function ##
################################################
function gen_interaction(u, i, p)
    return dot(u, p[:A][i, :])
end

################################################
## Generalised Lotka Volterra formula
################################################
function gen_lotka_volterra!(du, u, par, t)
    for i in eachindex(u)
        if t == 0 #due to issues with extinction_affect! modifying input ex, resetting all ex_ctrl to 1.0 at beginning of run required
            par[:ex_ctrl][i] = 1.0
        end
        du[i] = u[i] * (par[:r][i] + gen_interaction(u, i, par))
    end
end

################################################
## Harvesting form of Generalised Lotka Volterra formula
################################################
function gen_lotka_volterra_harvest!(du, u, par, t)
    C = par[:signal](t)
        for i in eachindex(u)
        if t == 0 #due to issues with extinction_affect! modifying input ex, resetting all ex_ctrl to 1.0 at beginning of run required
             par[:ex_ctrl][i] = 1.0
        end
        if i in par[:spp] #invoke control parameter on target species
            if t < 200
                C = 0.0
            else 
                C = C
            end
            du[i] = (u[i] * (par[:r][i] + gen_interaction(u, i, par))) - (C * (u[i]^2 / (u[i]^2 + 1^2)))
            else 
            du[i] = (u[i] * (par[:r][i] + gen_interaction(u, i, par)))
            end
    end
   end

function gen_lotka_volterra_harvest_remake!(du, u, par, t)
    C = par[:signal]
        for i in eachindex(u)
        if t == 0 #due to issues with extinction_affect! modifying input ex, resetting all ex_ctrl to 1.0 at beginning of run required
             par[:ex_ctrl][i] = 1.0
        end
        #if any(in.(par[:spp], Ref(i))) #invoke control parameter C if i equals species intended to be affected
        if i in par[:spp]
            if t < 200
                C = 0.0
            else 
                C = C
            end

            du[i] = (u[i] * (par[:r][i] + gen_interaction(u, i, par))) - (C * (u[i]^2 / (u[i]^2 + 1^2)))
            else 
            du[i] = (u[i] * (par[:r][i] + gen_interaction(u, i, par)))
            end
    end
   end


   function gen_lotka_volterra_harvest_v2!(du, u, par, t)
    C = par[:signal](t)
        for i in eachindex(u)
        if t == 0 #due to issues with extinction_affect! modifying input ex, resetting all ex_ctrl to 1.0 at beginning of run required
             par[:ex_ctrl][i] = 1.0
        end
        if any(in.(par[:spp], Ref(i))) #invoke control parameter C if i equals species intended to be affected
        du[i] = (u[i] * (par[:r][i] + gen_interaction(u, i, par))) - (C * (u[i]^2 / (u[i]^2 + 1^2)))
            else 
            du[i] = (u[i] * (par[:r][i] + gen_interaction(u, i, par))) - ((u[i]^2 / (u[i]^2 + 1^2)))
            end
    end
   end
   
function stoc_gen_lv!(du, u, par, t)
    for i in eachindex(u)
        du[i] = u[i]*par[:ex_ctrl][i]*par[:W][i] #for multiplicative noise
        #du[i] = par[:ex_ctrl][i]*par[:W][i] #for additive noise
    end
end

function stoc_gen_lv_delayed!(du,u, par, t)
    if t < 150
        du .= 0.0
    else
        for i in 1:length(u)
        du[i] = par[:ex_ctrl][i]*par[:W][i] #for additive noise
        #du[i] = par[:ex_ctrl][i]*rand((0:0.001:0.05),1)[1] #for brownian noise in combination with W
        #w = ((0*rand(Distributions.Normal(0,0.25),1))+(sqrt(1-0^2)*rand(Distributions.Normal(0,0.25),1)))
        #du[i] = par[:ex_ctrl][i]*w[1]
        end
    end
end

#function extinction!(out,u,t,integrator) # Event when event_f(u,t) == 0
#    for i in eachindex(u)
#    out[i] = 0.001-u[i] #0.01 is the lowest tolerable density
#    end
#end
 
#function extinction_affect!(integrator,idx)
#    integrator.u[idx] = 0 #from the idx of extinction!, set those integrator u's to 0
#    integrator.p[:ex_ctrl][idx] = 0 #repeat for noise toggle (prevents de-extinction)
#  end

# extinction_event = VectorContinuousCallback(extinction!,extinction_affect!,4)

#A_4 = [-1 -1.09 -1.52 -0; 
 #    -0 -0.72 -0.3168 -0.9792; 
 #    -3.5649 -0 -1.53 -0.7191;
 #    -1.5367 -0.6477 -0.4445 -1.27]
#r_4 = [1 0.72 1.53 1.27] #1.52
#ex = [1.0 1.0 1.0 1.0]

#u0_4 = 1.2*rand(4)

#u0_4 = [0.1, 0.5, 0.3, 0.4]
#tspan = (0.0,100)
#u0_4 = [0.4 0.5 0.115 0.4] #0.27
#ctrl = [0.5 0.5 0.5 0.5]
#params_4 = Dict(:A => A_4, :r => r_4,:ex_ctrl => ex)
#W = WienerProcess(0.0,0.0,0.0)

#stoc_gen_prob = SDEProblem(gen_lotka_volterra!, stoc_gen_lv!, u0_4, tspan, params_4,callback = extinction_event,dt=0.1)
#stoc_gen_sol = solve(stoc_gen_prob,EM())
#plot(stoc_gen_sol)

#df_out = DataFrame(stoc_gen_sol)
#CSV.write("/Users/duncanobrien/Downloads/lv_run1.csv",df_out)
