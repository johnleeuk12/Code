function Conductance_IFmodel(IE_delay,E_input,E_strength)




if IE_delay>=0
    Ge=E_input*E_strength;
    Gi=[zeros(1,delay) I_input(1:(length(I_input)-delay))]*I_strength;
elseif IE_delay<0
    Gi=I_input*I_strength;
    Ge=[zeros(1,delay) E_input(1:(length(E_input)-delay))]*E_strength;
end

%add pre  stim time of 500 ms
Ge=[zeros(size(0:step:PREstimulus_duration)) Ge];
Gi=[zeros(size(0:step:PREstimulus_duration)) Gi];
[spikes,V]=LIFmodel(Ge,Gi);
plot(V) 
end


function [spikes,V]=LIFmodel(Ge,Gi)
spikes=[]; V=[]; t=1; i=1;
step=.0001; %.1 ms duration  (temporal increment for running simulation)
C=0.25*1e-9; %0.25 nF, 10 ms time constant
Grest=25*1e-9; %25 nS
Erest=-0.065; %-65 mV
Ee=0; % 0 mV
Ei=-0.085;  %-85 mV
V(1)=Erest;
noise_magnitude=4e-8; %default noise level in conductance

%avoid negative conductances
Ge=Ge+noise_magnitude*randn(1,length(Ge));
Gi=Gi+noise_magnitude*randn(1,length(Gi));
Ge(find(Ge<0))=0;
Gi(find(Gi<0))=0;

while(t<length(Ge))
    V(t+1)=(-step*(Ge(t)*(V(t)-Ee)+Gi(t)*(V(t)-Ei)+Grest*(V(t)-Erest))/C)+V(t);
    if V(t+1)>(Erest+0.020) %20 mV above Erest %artificial threshold
        V(t+1)=0.050; %spike to 50 mV
        spikes(i)=step*(t+1);
        i=i+1;
        t=t+1;
        V(t+1)=Erest;
    end
    t=t+1;
end
