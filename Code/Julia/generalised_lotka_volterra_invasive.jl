using DifferentialEquations
using Plots
using LinearAlgebra
using Interpolations #interpolate for harvesting/control parameters
using DataFrames #create dataframes
using CSV #export CSV files
using Distributions #uniform function


function gen_interaction_k(u, i, p)
    return sum(u.*p[:A][i, :])    #identical
    #return dot(u, p[:A][i, :])
end

function gen_lotka_volterra_invasive_harvest!(du, u, par, t)
    M = par[:signal](t)
    #for i in 1:length(u)
    for i in eachindex(u)
        # if t == 0 #due to issues with extinction_affect! modifying input ex, resetting all ex_ctrl to 1.0 at beginning of run required
        #     par[:ex_ctrl][i] = 1.0
        # end
        if i in par[:spp]
            if t<200.0
              Kt = par[:k][i]
            else
        Kt = par[:k][i]*(1.0 + (par[:eta][i]*M))
            end
        else
            Kt = par[:k][i]
        end
            du[i] = (u[i] * (par[:r][i]) * (1.0-(gen_interaction_k(u, i, par)/Kt))) + par[:imm][i]
            #du[i] = (u[i] * (par[:r][i]) * (1.0-(gen_interaction_k(u, i, par)))) 
            #du[i] = (u[i] * (par[:r][i]) * (Kt-(gen_interaction_k(u, i, par))))/Kt + par[:imm][i]

        end
end

function gen_lotka_volterra_invasive_remake!(du, u, par, t)
    M= par[:signal]
    for i in eachindex(u)
        # if t == 0 #due to issues with extinction_affect! modifying input ex, resetting all ex_ctrl to 1.0 at beginning of run required
        #     par[:ex_ctrl][i] = 1.0
        # end
        #if any(in.(par[:spp], Ref(i))) #invoke control parameter C if i equals species intended to be affected
        if i in par[:spp]
            if t<200
              Kt = par[:k][i]
            else
        Kt = par[:k][i]*(1 + (par[:eta][i]*M))
            end
        else
            Kt = par[:k][i]
        end
            du[i] = (u[i] * (par[:r][i]) * (1-(gen_interaction_k(u, i, par)/Kt))) + par[:imm][i]
            #du[i] = (u[i] * (par[:r][i]) * (Kt-(gen_interaction_k(u, i, par))))/Kt + par[:imm][i]

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
        du[i] = par[:ex_ctrl][i]*par[:W] #for additive noise
        #du[i] = par[:ex_ctrl][i]*rand((0:0.001:0.05),1)[1] #for brownian noise in combination with W
        #w = ((0*rand(Distributions.Normal(0,0.25),1))+(sqrt(1-0^2)*rand(Distributions.Normal(0,0.25),1)))
        #du[i] = par[:ex_ctrl][i]*w[1]
        end
    end
end

function build_A_dakos!(n::Int64)
    n = convert(UInt8,n)
    A = zeros(n, n)
    A[ UpperTriangular(trues(size(A))) ] = rand(Distributions.Uniform(0.0,1.5),convert(Int,(n * (n - 1)/ 2) + n)) #randomly populate upper triangle
    #A = Symmetric(A,:U)
    A[ LowerTriangular(trues(size(A))) ] = A[ UpperTriangular(trues(size(A))) ]
    A[ diagind(A) ] .= 1
    return A
end