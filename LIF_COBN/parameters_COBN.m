% Conventions for variable names:
% - prefix "N" stands for "number of..."
%
% - prefix "e" stands for "excitatory"
% - prefix "i" stands for "inhibitory"
% - prefix "x" stands for "external"
%
% - prefix "e2e" stands for "excitatory to excitatory"
% - prefix "e2i" stands for "excitatory to inhibitory"
% - prefix "x2e" stands for "external to excitatory"
% - prefix "e2i" stands for "external to inhibitory"
%
% - "T" stands for greek letter "tau"



% THE NETWORK

% Time resolution [ms]
net_COBN.Dt = 0.05;


% NETWORK PROPERTIES -----------------------------------------------------

% Number of excitatory and inhibitory neurons, eNnrn and iNnrn, and total number
% of neurons, totNnrn:
net_COBN.eNnrn = 4000;
net_COBN.iNnrn = 1000;

% Connectivity, i.e., connection probability, p:
net_COBN.p = 0.2; 

% Membrane time constant, [ms]
net_COBN.eTm = 20; 
net_COBN.iTm = 10; 

% Leak membrane potential, [mV]
net_COBN.V_leaky = -70.; 

% Membrane potential firing threshold, [mV]
net_COBN.Vthr = -52; 
% the neuron fires following the following sequence:
% 1. the neuron potential is reset to a value Vres, [mV]: 
net_COBN.eVres = -59; 
net_COBN.iVres = -59; 
% 2. the neuron cannot fire again for a refractory period, Trp, equal to
%    2ms for E neurons and 1ms for I neurons, [ms]:
net_COBN.eTrp = 2; 
net_COBN.iTrp = 1; 
% 3. all post-synaptic neurons receive a spike with a delay, Tl, equal to
%    of 1ms after the time of threshold crossing, [ms]:
net_COBN.eTl = 1;
net_COBN.iTl = 1; 

% Rise and deacy times, Tr and Td [ms]:
% - of E => E synaptic currents:
net_COBN.e2eTr = 0.4;
net_COBN.e2eTd = 2.; 
% - of E => I synaptic currents:
net_COBN.e2iTr = 0.2; 
net_COBN.e2iTd = 1.; 
% - of I synaptic currents: rise and decay times are assumed idependent of
%   the type of neuron they act on (excitatory or inhibitory)
net_COBN.iTr = 0.25; 
net_COBN.iTd = 5.; 

% Synaptic reversal potentials, [mV]
net_COBN.VsynAMPA = 0; 
net_COBN.VsynGABA = -80;

% Synaptic conductances [nS]
% - on inhibitory neurons:
net_COBN.gi2i = 2.698602679456193;
net_COBN.ge2i = 0.233373613159943; 
net_COBN.gx2i = 0.316721332145637; 
% - on excitatory neurons:
net_COBN.gi2e = 2.008771996214003; 
net_COBN.ge2e = 0.178441556321102; 
net_COBN.gx2e = 0.233673466610967; 

% Membrane resistances, [GOhm]
net_COBN.eRm = 0.04; 
net_COBN.iRm = 0.05; 




