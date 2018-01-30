%% Alpha function and synaptic inputs 
function out = alpha_synapse(Amplitude, time_constant, time)

out = Amplitude*time*exp(-time/time_constant);
end