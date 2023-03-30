#### Extinction Event Callback ####
using DifferentialEquations
using LinearAlgebra

function univariate_extinction!(out,u,t,integrator) # Event when event_f(u,t) == 0
    out = 0.01-u #0.01 is the lowest tolerable density
end
 
function univariate_extinction_affect!(integrator)
    integrator.u = 0 #from the idx of extinction!, set those integrator u's to 0
    integrator.p[:ex_ctrl] = 0 #repeat for noise toggle (prevents de-extinction)
  end

univariate_extinction_event = ContinuousCallback(univariate_extinction!,univariate_extinction_affect!)

function extinction!(out,u,t,integrator) # Event when event_f(u,t) == 0
  #for i in 1:length(u)
  for i in eachindex(u)

  out[i] = 0.001-u[i] #0.01 is the lowest tolerable density
  end
end

function extinction_affect!(integrator,idx)
  integrator.u[idx] = 0 #from the idx of extinction!, set those integrator u's to 0
  integrator.p[:ex_ctrl][idx] = 0 #repeat for noise toggle (prevents de-extinction)
end

function extinction_affect_buffer!(integrator,idx)
  integrator.u[idx] = integrator.uprev[idx] #from the idx of extinction!, set those integrator u's to 0
end

extinction_event = VectorContinuousCallback(extinction!,extinction_affect!,4)

extinction_event_buffer = VectorContinuousCallback(extinction!,extinction_affect_buffer!,4)
