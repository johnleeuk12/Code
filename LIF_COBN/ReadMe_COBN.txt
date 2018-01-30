Recurrent random network with excitatory and inhibitory Leaky Integrate-and-Fire (LIF) neurons with 
conductance-based synapses from the paper “Comparison of the dynamics of neural interactions between 
current-based and conductance-based integrate-and-fire recurrent networks” written by S.Cavallari, S.Panzeri 
and A.Mazzoni and published in Frontiers in Neural Circuits (2014), 8:12. doi:10.3389/fncir.2014.00012. The 
paper compares the activity of this conductance-based network (i.e. code_COBN.c) with the activity of a 
comparable network of LIF neurons with current-based synapses (whose source code is in the “LIF_CUBN” 
folder). 

The function code_COBN.c is a mex source code. You have to compile this routine in Matlab to generate the 
mex file (e.g.: code_COBN.mexw64). Note that you have to include the functions ran1.c and gasdev.c in the 
compiling instruction in the Matlab workspace, in the following way: 
mex code_COBN.c ran1.c gasdev.c 

After you compiled the function, you can call it as specified in the help. 
For more information use the help of the function (and see the example below): 
help code_COBN 

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
In the following an easy example for setting the arguments to generate the data used in figures 4I; 7 (i.e. LFP) 
and 6A (i.e. average firing rate) once you compiled the mex file.

In the Matlab workspace:

1.	parameters_COBN; % to generate the structure net_COBN with all the parameters of the network
2.	simulation_length = 4500; % units: (ms)
3.	M = simulation_length/net_COBN.Dt; % length of the simulation in time steps
4.	external_signal_intensity = 2; % units; (spikes/ms)/cell 
5.	external_signal = ones(M,1) * external_signal_intensity * net_COBN.Dt;
6.	SEED_OU = 1; % positive integer number
7.	external_noise = OU_process(M, net_COBN.Dt, 16, 0.16*net_COBN.Dt, SEED_OU);
8.	INPUT2E = external_signal + external_noise;
9.	INPUT2I = INPUT2E;
10.	SEED_connections = 2; % positive integer number
11.	SEED_poisson = 3; % positive integer number
12.	[E2EI,I2EI,eFR,iFR] = code_COBN(net_COBN, INPUT2E, INPUT2I, SEED_connections, SEED_poisson);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

Question on how to use the model should be addressed to:
ste.cavallari@gmail.com

Please cite the paper if you use the code.
