function Rec_Network()






%
% fr_ind = [2.6 2.7 2.8 2.9];
% parfor i = 1:4
%     fr = fr_ind(i);
COBN();
% Sync();




function [V,spikes] = Sync()

%% Parameters
Grest(1)    =   25*1e-9;            %25 nS
Grest(2)    =   20*1e-9;
C(1)        =   0.25*1e-9;          %0.25 nF, 10 ms time constant excitatory
C(2)        =   Grest(2)*tau_m(2);
tau_m(1)    =   0.01;
tau_m(2)    =   0.02;

Erest       =   -0.070;             %-65 mV
Ereset      =   -0.059;
Vth         =   -0.052;             %-50 mV
% refract(1)  =   .002;               %refractory period for excitatory and inhibitory
% refract(2)  =   .001;
E(1)        =   0;                  % 0 mV
E(2)        =   -0.085;             %-85 mV
sigma       =   0.01;
% noise_magnitude=4e-8; %default noise level in conductance
step        =   .0001;              %.1 ms duration  (temporal increment for running simulation)
N           =   10;                 % total number of neurons
Ne          =   floor(N*0.8);       %number of excitatory neurons. 
total_time  =   1.5;                % seconds
CX          =   0.1;                %connectivity 

tau_l       =   1*1e-3;                  %Synaptic latency
tau_r(1)    =   0.2*1e-3;           %for AMPA on inhibitory
tau_r(2)    =   0.25*1e-3;          %Synaptic rise time for GABA
tau_r(3)    =   0.4*1e-3;           %for AMPA on excitatory
tau_d(1)    =   1*1e-3; 
tau_d(2)    =   5*1e-3;             %Synaptic decay time
tau_d(3)    =   2*1e-3;

% tau_s       =   1.7*1e-3;

g_syn(1,1)  =   0.178*1e-9;         %AMPArec on excitatory 
g_syn(2,1)  =   0.233*1e-9;         %AMPArec on inhibitory 
g_syn(1,2)  =   2.01*1e-9;          %Synaptic conductance GABA on excitatory
g_syn(2,2)  =   2.70*1e-9;          %GABA on inhibitory
g_syn(1,3)  =   0.234*1e-9;         %AMPAext on excitatory
g_syn(2,3)  =   0.317*1e-9;         %AMPAext on inhibitory




stimulus_duration=0.5;  %half second
PREstimulus_duration=0.5;  %half second
POSTstimulus_duration=0.5;  %half second (0.1 second is included in stimulus)

nreps=10; % for each ICI
ICI_list= [125 83.3333 62.5 50 41.6667 35.7143 31.25 27.7778 25 22.7273 20.8333];

IE_delay = 5.; % in ms
E_strength = 4.5; % in nS
IE_ratio = 6;

I_strength=8.5;  %IE_ratio*E_strength;
kernel_time_constant=.005;  %time constant of 5 ms
jitter_magnitude=1; % Temporal jitter
step=.0001; %.1 ms duration  (temporal increment for running simulation)
stimulus_duration=0.5;  %half second
PREstimulus_duration=0.5;  %half second
POSTstimulus_duration=0.5;  %half second (0.1 second is included in stimulus)
latency_time=0.01;
latency=length(0:step:latency_time); %10 ms latency (auditory nerve to auditory cortex)

%% acoustic pulse train stimulus
%%%%%%%%%acoustic pulse train stimulus%%%%%%%%%%%%%%%%%

nreps=10; % for each ICI
ICI_list= [125 83.3333 62.5 50 41.6667 35.7143 31.25 27.7778 25 22.7273 20.8333]; %[3 5 7.5 10 12.5 15 20 25 30 35 40 45 50 55 60 65 70 75];

freq_list=1000./ICI_list; % Envelope Frequency
spike_distribution=NaN(length(ICI_list),nreps);
time_distribution=NaN(length(ICI_list),nreps);
spont_distribution=NaN(1,length(ICI_list)*nreps);
raster.stim=[];  raster.rep=[];  raster.spikes=[];

kernel = zeros(1,501);
j = 3;
for t = 1:501
    if t>round(tau_l/step)
        
        kernel(1,t) = (exp(-(t*step-tau_l)/tau_d(j))-exp(-(t*step-tau_l)/tau_r(j)));
    end
end




for f=1:length(freq_list)
    spikes_pooled=[];
    freq=freq_list(f);
%     t=0:step:(kernel_time_constant*10);
%     kernel=t.*exp(-t/kernel_time_constant);
%     kernel=1e-9*kernel/max(kernel); %amplitude of 1 nS
    input=zeros(size(0:step:(POSTstimulus_duration+stimulus_duration)));
    stimulus_input_length=length(0:step:(stimulus_duration));
    ipi=round(1/(freq*step)); %ipi=interpulse interval
    freq2=1/(step*ipi);
    rate_total = [];
    for r=1:nreps
        E_input=input;
        I_input=input;
        for j=1:10  %10 jitter excitatory and inhibitory inputs
            % for i=1:ipi:(stimulus_input_length-(length(kernel)/2))
            for i=1:ipi:(stimulus_input_length-250)
                jitter=round(randn(1)/(1000*step)); %1 ms jitter
                if (jitter+i)<1 || (jitter+i)>(length(input)-length(kernel))
                    jitter=0;
                end
                E_input((latency+i+jitter):(latency+i+jitter+length(kernel)-1))=E_input((latency+i+jitter):(latency+i+jitter+length(kernel)-1))+tau_m(1)/(tau_d(3)-tau_r(3))*kernel;
                
                
                jitter=round(randn(1)/(1000*step)); %1 ms jitter
                
                if (jitter+i)<1 || (jitter+i)>(length(input)-length(kernel))
                    jitter=0;
                end
                I_input((latency+i+jitter):(latency+i+jitter+length(kernel)-1))=I_input((latency+i+jitter):(latency+i+jitter+length(kernel)-1))+tau_m(2)/(tau_d(1)-tau_r(1))*kernel;
            end
        end
        
        delay=round(abs(IE_delay)/(1000*step));  %delay in steps
        if IE_delay>=0
            Ge=E_input*E_strength*1e-9;
            Gi=[zeros(1,delay) I_input(1:(length(I_input)-delay))]*I_strength*1e-9;
        elseif IE_delay<0
            Gi=I_input*I_strength*1e-9;
            Ge=[zeros(1,delay) E_input(1:(length(E_input)-delay))]*E_strength*1e-9;
        end
        
        %add pre  stim time of 500 ms
        Ge=[zeros(size(0:step:PREstimulus_duration)) Ge];
        Gi=[zeros(size(0:step:PREstimulus_duration)) Gi];
        [spikes,V]=run_LIFmodel(Ge,Gi);
        rate = zeros(1,1500);
        for st = spikes
            rate(1,round(st/(step*10))) = rate(1,round(st/(step*10)))+1;
        end
        rate_total = [rate_total ; rate*1000];
        spikes=spikes-PREstimulus_duration;
        spikes_pooled=[spikes_pooled spikes];
        
        spike_distribution(f,r)=length(spikes(find(spikes>0 & spikes<=(stimulus_duration+0.1))))/(stimulus_duration+0.1);
        spont_distribution(r+nreps*(f-1))=length(spikes(find(spikes>-PREstimulus_duration & spikes<0)))/PREstimulus_duration;
        
        raster.stim=[raster.stim f*ones(size(spikes))];
        raster.rep=[raster.rep r*ones(size(spikes))];
        raster.spikes=[raster.spikes spikes];
    end
    
    rate_av = mean(rate_total,1);
    spikes_pooled_for_vector_strength=spikes_pooled(find(spikes_pooled>0.05 & spikes_pooled<=(stimulus_duration+0.05)));
   
    %for calculating vector strength, subtract the first 50 ms and
    %include 50 ms post stimulus
    if ~isempty(spikes_pooled_for_vector_strength)
        total_spikes=length(spikes_pooled_for_vector_strength);
        x=0;
        y=0;
        if total_spikes>0
            x=sum(cos(2*pi*(spikes_pooled_for_vector_strength*freq2)));
            y=sum(sin(2*pi*(spikes_pooled_for_vector_strength*freq2)));
        end
        if total_spikes==0
            vector(f)=0;
        else
            vector(f)=sqrt(x^2+y^2)/total_spikes;
        end
        rayleigh(f)=2*total_spikes*vector(f)^2;
    else
        vector(f)=0;
        rayleigh(f)=0;
    end
    
    if rayleigh(f)<13.8
        vector(f)=0;
    end
    %spike rate calculated over stimulus duration plus 100 ms post stimulus
    spike_rate(f)=length(find(spikes_pooled>0 & spikes_pooled<=(stimulus_duration+0.1)))/(nreps*(stimulus_duration+0.1));
end




%% LIF model spikes and MP simulation


function [spikes,V]=run_LIFmodel(Ge,Gi)
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



function [V, spikes] = COBN()
%% Parameters
Grest(1)    =   25*1e-9;            %25 nS
Grest(2)    =   20*1e-9;
C(1)        =   0.25*1e-9;          %0.25 nF, 10 ms time constant excitatory
tau_m(1)    =   0.01;
tau_m(2)    =   0.02;
C(2)        =   Grest(2)*tau_m(2);
Erest       =   -0.070;             %-65 mV
Ereset      =   -0.059;
Vth         =   -0.052;             %-50 mV
% refract(1)  =   .002;               %refractory period for excitatory and inhibitory
% refract(2)  =   .001;
E(1)        =   0;                  % 0 mV
E(2)        =   -0.085;             %-85 mV
sigma       =   0.01;
% noise_magnitude=4e-8; %default noise level in conductance
step        =   .0001;              %.1 ms duration  (temporal increment for running simulation)
N           =   10;                 % total number of neurons
Ne          =   floor(N*0.8);       %number of excitatory neurons. 
total_time  =   1.5;                % seconds
CX          =   0.1;                %connectivity 

tau_l       =   1*1e-3;                  %Synaptic latency
tau_r(1)    =   0.2*1e-3;           %for AMPA on inhibitory
tau_r(2)    =   0.25*1e-3;          %Synaptic rise time for GABA
tau_r(3)    =   0.4*1e-3;           %for AMPA on excitatory
tau_d(1)    =   1*1e-3; 
tau_d(2)    =   5*1e-3;             %Synaptic decay time
tau_d(3)    =   2*1e-3;

% tau_s       =   1.7*1e-3;

g_syn(1,1)  =   0.178*1e-9;         %AMPArec on excitatory 
g_syn(2,1)  =   0.233*1e-9;         %AMPArec on inhibitory 
g_syn(1,2)  =   2.01*1e-9;          %Synaptic conductance GABA on excitatory
g_syn(2,2)  =   2.70*1e-9;          %GABA on inhibitory
g_syn(1,3)  =   0.234*1e-9;         %AMPAext on excitatory
g_syn(2,3)  =   0.317*1e-9;         %AMPAext on inhibitory

%% Initializing Variables

timebin  = 0:step:total_time;
V = zeros(N,length(timebin));
Itot = zeros(N,length(timebin)); %total synaptic current
I = zeros(N,3,length(timebin)); %AMPArec %GABA %AMPAext
s = zeros(N,3,length(timebin));
Poisson = zeros(N,length(timebin));

O = zeros(N,length(timebin));
V(:,1) = Erest;
W = zeros(N,N);
for i = 1:N
    for j = 1:N
        if rand() <= CX
            if i<=Ne && j<=Ne % AMPArecurrent on Excitatory
                W(i,j) = 1;
            elseif i >=Ne && j<=Ne %AMPArec on Inhibitory
                W(i,j) = 2;
            elseif i <=Ne && j>=Ne %GABA on excitatory
                W(i,j) = 3;
            else
                W(i,j) = 4;  
            end
        end
    end
end




%% Network

% fr = fr*1e3; %spikes/sec
for n = 1:N
    spikes{n} = []; %spike times 
    %generate Poisson input for all 100 neurons 
%     Poisson(n,:) = rand(1, length(timebin)) < fr*step;
end


%spont rate, noise

noise_magnitude=40; %default noise level in conductance
%avoid negative conductances
for n = 1:N
    for j = 1:3
        s(n,j,:) = s(n,j,:) + noise_magnitude*randn(1,1,length(timebin));
        for t = 1:length(timebin)
            if s(n,j,t) < 0
                s(n,j,t) = 0;
            end
        end
    end
end



%define synaptic kinetics kernels
% hold on
kernel = zeros(1,4,501);
for j = 1:3
    for t = 1:501
        if t>round(tau_l/step)
            
            kernel(1,j,t) = (exp(-(t*step-tau_l)/tau_d(j))-exp(-(t*step-tau_l)/tau_r(j)));
        end
    end
end
  


% figure
% plot(kernel(2,:))

tic




for t = 1:length(timebin)-501 %taking into account the synaptic input kernel. 
    for n = 1:N
        if n <= Ne % Excitatory or Inhibitory neurons 
            i = 1;
        else
            i = 2;
        end
        
        for j = 1:2  %1 = AMPArec 2 = GABA
            n_pre_list = find(W(n,:) == i + j - mod(i,2));
            
            if ~isempty(n_pre_list)
                for n_pre = n_pre_list
                    if ~isempty(spikes{n_pre})
                        for nn = 1: length(spikes{n_pre})
                            if t == max(spikes{n_pre}(nn))

                            s(n,j,t:t+500) = s(n,j,t:t+500) + tau_m(i)/(tau_d(j)-tau_r(j))*kernel(1,j,1:501);
                            end

                        end
                    end
                end
            end
            I(n,j,t) = g_syn(i,j)*s(n,j,t)*(V(n,t)-E(i));
        end
%         t_pre_list = find(Poisson(n,:) ==1);
        %         t_pre_list(find(t_pre_list>t)) =[];
        
        j = 3; %AMPA ext. (stim and spont rate)
%         if ismember(t, t_pre_list)==1
%             s(n,j,t:t+500) = s(n,j,t:t+500) +tau_m(i)/(tau_d(j)-tau_r(j))*kernel(1,j,1:501);
%         end
        
        
        %         if isnan(s(n,j,t)) ==1 
        %             test =1;
        %         end
        
        I(n,3,t) = g_syn(i,3)*s(n,3,t)*(V(n,t)-E(1));
        
        Itot(n,t) =  I(n,1,t) + I(n,2,t) + I(n,3,t);
        if ~isempty(spikes{n})
             if t-spikes{n}(end) >=2 %refractory period
                V(n,t+1) = -step/C(i)*(Grest(i)*(V(n,t)-Erest) + Itot(n,t)) +V(n,t)+sigma*randn*sqrt(step) ;
            else
                V(n,t+1) = Ereset + sigma*randn*sqrt(step);
            end
        else
            V(n,t+1) = -step/C(i)*(Grest(i)*(V(n,t)-Erest) + Itot(n,t)) +V(n,t)+sigma*randn*sqrt(step) ;
        end
        if V(n,t+1) > Vth
            if ~isempty(spikes{n})
                if t > spikes{n}(end)+(mod(i,2)+1)*1e-3/step %refractory period
                    
                    V(n,t+1) = 0.050; %spike to 50 mV
                    spikes{n} = [spikes{n} t];
                    O(n,t) = 1;
                end
            else
                V(n,t+1) = 0.050; %spike to 50 mV
                spikes{n} = [spikes{n} t];
                O(n,t) = 1;
            end
            
        end
    end
end
toc

%Statistics

rate = zeros(N,1500);
for n = 1:N    
    for st = spikes{n}
        rate(n,round(st/10)) = rate(n,round(st/10))+1000;
    end
end

rate_av = mean(rate(1:Ne,:),1);
xs = 1:1500;
h = 10; %kernal bandwidth. determines the shape of the function
ys = zeros(1,1500);
for i = 1:1500
    ys(i)=gaussian_kern_reg(xs(i),xs,rate_av,h);
end

% figure
% % title(num2str(fr))
% plot(ys);
% % namefig = ['ys_' num2str(fr)];
% savefig(namefig)


figure
for i = 1:10
subplot(10,1,i)
plot(V(i,:))
end

test = 1




















    