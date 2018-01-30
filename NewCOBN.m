function NewCOBN()


% spont = 9; %spont rate = 5spikes/sec


stim_time = 200*1e-3;
% stim_list = [3 4 5 6 7 8 9 10]*1e3;
% for i = 1:length(stim_list)
%     stim = stim_list(i);
%     disp(spont)
stim = 0.2*1e3;
% parfor i = [1 2 3 4]
% spont = 7.5;
% param2_list = 0.16:0.01:0.17;
% parfor i = 1:length(param2_list)
%     param2 = param2_list(i);
%     for spont = 10.6:0.1:11.5
%         for param2 = 0.15:0.02:0.3
% 
%             network(spont,stim,stim_time,param2);
%                     pause(0.1)
%         end
%     end
% end

% param1 = 9;
param2 = 0.2;
spont = 7.5;
% 
network(spont,stim,stim_time,param2);

% end
% end
%     pause(0.1);
% end

% network(spont,3*1e3,5*1e-3);


function [V,spikes] = network(spont,stim_fr,input_time,param2)
%% Parameters
Grest(1)    =   25*1e-9;            %25 nS
Grest(2)    =   20*1e-9;
tau_m(1)    =   0.02;%0.01;
tau_m(2)    =   0.01;%0.02;
C(1)        =   0.5*1e-9;%0.25*1e-9;          %0.25 nF, 10 ms time constant excitatory
C(2)        =   0.2*1e-9;%Grest(2)*tau_m(2);


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
N           =   1000;                 % total number of neurons
Ne          =   floor(N*0.9);       %number of excitatory neurons.
total_time  =   1.5;                % seconds
CX          =   0.05; %0.1;                %connectivity



tau_l       =   3*1e-3;             %Synaptic latency
tau_rp      =   3*1e-3;             %refractory period. 
%Alpha function model.

% tau_s(1)    =   1.5*1e-3;           %for AMPA on inhibitory. all fitted to exp kernels to imitate shape.
% tau_s(2)    =   1.7*1e-3;           %Synaptic time constant for GABA
% tau_s(3)    =   2.0*1e-3;           %for AMPA on excitatory

% % %Exponential model
tau_r(1)    =   0.2*1e-3;           %for AMPA on inhibitory
tau_r(2)    =   0.25*1e-3;          %Synaptic rise time for GABA
tau_r(3)    =   0.4*1e-3;           %for AMPA on excitatory
tau_d(1)    =   1*1e-3;
tau_d(2)    =   5*1e-3;             %Synaptic decay time
tau_d(3)    =   2*1e-3;


g_syn(1,1)  =   0.178*1e-9;         %AMPArec on excitatory
g_syn(2,1)  =   0.233*1e-9;         %AMPArec on inhibitory 
g_syn(1,2)  =   2.01*1e-9;          %Synaptic conductance GABA on excitatory  g_syn(1,1)*param1;  %1.8*1e-9; 
g_syn(2,2)  =   2.70*1e-9;          %GABA on inhibitory g_syn(2,1)*param1;  %3*1e-9; 
g_syn(1,3)  =   0.234*1e-9;         %AMPAext on excitatory
g_syn(2,3)  =   param2*1e-9;        %0.317*1e-9;         %AMPAext on inhibitory


% adaptation

tauP(1) = 5*1e-3;
tauP(2) = 30*1e-3;
f_D = 0.6;
P0= 1;


%% Initializing Variables

timebin  = 0:step:total_time;
V = zeros(N,length(timebin));
Itot = zeros(N,length(timebin)); %total synaptic current
I = zeros(N,3,length(timebin)); %AMPArec %GABA %AMPAext
s = zeros(N,3,length(timebin));

O = zeros(N,length(timebin));
V(:,1) = Erest;
W = zeros(N,N);
P_rel = ones(N,length(timebin));
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

% imagesc(W)

%% Network

kernel = zeros(1,3,501);
%exp kernel
for j = 1:3
    for t = 1:501
        if t>round(tau_l/step)
            kernel(1,j,t) = (exp(-(t*step-tau_l)/tau_d(j))-exp(-(t*step-tau_l)/tau_r(j)));
        end
    end
end
%
%
% test = 1;
% for t = 1:501
%     test1(t) = kernel(1,1,t);
%     test2(t) = kernel(1,2,t);
%     test3(t) = kernel(1,3,t);
% end
%
% figure
% plot(test1)
% hold on
% plot(test2)
% plot(test3)

% test = 1;
% Alpha kernel.

% for j = 1:3
%     for t = 1:501
%         if t>round(tau_l/step)
%             kernel(1,j,t) = (t*step).*exp(-t*step/tau_s(j));
%         end
%     end
%     kernel(1,j,:) = kernel(1,j,:)/max(kernel(1,j,:));
% end




%spont rate
% fr = fr*1e3; %spikes/sec
for n = 1:N
    spikes{n} = []; %spike times
    %     Poisson(n,:) = rand(1, length(timebin)) < fr*step; %generate Poisson input for all 100 neurons
end

noise_magnitude = spont*10;
for n = 1:N
    if n<=Ne
        s(n,3,:) = s(n,3,:) +  g_syn(1,3)*noise_magnitude*randn(1,1,length(timebin));
    else
        s(n,3,:) = s(n,3,:) +  g_syn(2,3)*noise_magnitude*randn(1,1,length(timebin));
    end
    for t = 1:length(timebin)
        if s(n,3,t) < 0
            s(n,3,t) = 0;
        end
    end
end

%end spont rate

%input
% input_time = 10*1e-3;
% fr = 3*1e3;
input_list = randsample(1:Ne,round((Ne)/10));
input_timebin = 0:step:input_time;
Poisson = zeros(N,length(input_timebin));


for n = 250:500
    if ismember(n,input_list)
        Poisson(n,:) = rand(1,length(input_timebin))< stim_fr*step;
        t_pre_list = find(Poisson(n,:) ==1);
        t_pre_list = t_pre_list + 5*1e3;
        for t = t_pre_list
            for nrep = 1:20
                s(n,3,t+1:t+501) = s(n,3,t+1:t+501) + g_syn(1,3)*(tau_m(1)/(tau_d(3)-tau_r(3)))*kernel(1,3,1:501);
            end
        end
    else
        Poisson(n,:) = zeros(1,length(input_timebin));
    end
end



tic
for t = 1:length(timebin)-501
    for n = 1:N
        if n<=Ne
            i = 1;
        else
            i = 2;
        end
        
        for j = 1:2
            %             n_pre_list = find(W(n,:) == i + j - mod(i,2));
            %
            %              if ~isempty(n_pre_list)
            %                 for n_pre = n_pre_list
            %                     if ~isempty(spikes{n_pre})
            %                         for nn = 1: length(spikes{n_pre})
            %                             if t == max(spikes{n_pre}(nn))
            %
            %                             s(n,j,t:t+500) = s(n,j,t:t+500) + g_syn(i,j)*kernel(1,j,1:501);
            %                             end
            %
            %                         end
            %                     end
            %                 end
            %              end
            I(n,j,t) = s(n,j,t)*(V(n,t)-E(j));
        end
        %         t_pre_list = find(Poisson(n,:) ==1);
        %                 t_pre_list(find(t_pre_list>t)) =[];
        
        %         j = 3; %AMPA ext. (stim and spont rate)
        %         if ismember(t, t_pre_list)==1
        %             s(n,j,t:t+500) = s(n,j,t:t+500) +g_syn(i,j)*kernel(1,j,1:501);
        %         end
        
        I(n,3,t) = s(n,3,t)*(V(n,t)-E(1));
        Itot(n,t) =  I(n,1,t) + I(n,2,t) + I(n,3,t);
        
        %membrane voltage change
        if ~isempty(spikes{n})
            if t-spikes{n}(end) >=2 %refractory period
                V(n,t+1) = -step/C(i)*(Grest(i)*(V(n,t)-Erest) + Itot(n,t)) +V(n,t)+sigma*randn*sqrt(step) ;
            else
                V(n,t+1) = Ereset + sigma*randn*sqrt(step);
            end
        else
            V(n,t+1) = -step/C(i)*(Grest(i)*(V(n,t)-Erest) + Itot(n,t)) +V(n,t)+sigma*randn*sqrt(step) ;
        end
        
        P_rel(n,t+1) = -step*(P_rel(n,t)-P0)/tauP(i) + P_rel(n,t); %P release for neuron n to excitatory neuron
        
        %spike?
        if V(n,t+1) > Vth
            if ~isempty(spikes{n})
                if t > spikes{n}(end)+ tau_rp/step%(mod(i,2)+1)*1e-3/step %refractory period
                    
                    V(n,t+1) = 0.050; %spike to 50 mV
                    spikes{n} = [spikes{n} t];
                    O(n,t) = 1;
                    P_rel(n,t+1) = f_D*P_rel(n,t);
                    for n_post = find(mod(W(:,n),2) == 1) %excitatory input to post synaptic neuron
                        if W(n_post,n) ~=3
                            s(n_post,i,t+1:t+501) = s(n_post,i,t+1:t+501) + g_syn(1,i)*(tau_m(i)/(tau_d(i)-tau_r(i)))*kernel(1,i,1:501).*P_rel(n_post,t);
                        else
                            s(n_post,i,t+1:t+501) = s(n_post,i,t+1:t+501) + g_syn(1,i)*(tau_m(i)/(tau_d(i)-tau_r(i)))*kernel(1,i,1:501);
                        end
                    end
                    for n_post = find(W(:,n) ~= 0 & mod(W(:,n),2) == 0) %inhibitory input to post synaptic neuron
                        s(n_post,i,t+1:t+501) = s(n_post,i,t+1:t+501) + g_syn(2,i)*(tau_m(i)/(tau_d(i)-tau_r(i)))*kernel(1,i,1:501).*P_rel(n_post,t);
                    end
                    
                    
                end
            else
                V(n,t+1) = 0.050; %spike to 50 mV
                spikes{n} = [spikes{n} t];
                O(n,t) = 1;
                P_rel(n,t+1) = f_D*P_rel(n,t);
                for n_post = find(mod(W(:,n),2) == 1) %excitatory input to post synaptic neuron
                    s(n_post,i,t+1:t+501) = s(n_post,i,t+1:t+501) + g_syn(1,i)*(tau_m(i)/(tau_d(i)-tau_r(i)))*kernel(1,i,1:501).*P_rel(n_post,t);
                end
                for n_post = find(W(:,n) ~= 0 & mod(W(:,n),2) == 0) %inhibitory input to post synaptic neuron
                    s(n_post,i,t+1:t+501) = s(n_post,i,t+1:t+501) + g_syn(2,i)*(tau_m(i)/(tau_d(i)-tau_r(i)))*kernel(1,i,1:501).*P_rel(n_post,t);
                end
            end
            
        end
    end
end


raster.spikesNe = [];
raster.neuronsNe = [];
raster.spikesNi = [];
raster.neuronsNi = [];
exci_list = randsample(1:Ne,100);
inhi_list = randsample(Ne+1:N,20);
% exci_list = 1:Ne;
% inhi_list = Ne+1:N;


i = 1;
for n = exci_list
    raster.spikesNe = [raster.spikesNe spikes{n}];
    raster.neuronsNe = [raster.neuronsNe i*ones(size(spikes{n}))];
    i = i+1;
end

for n = inhi_list
    raster.spikesNi = [raster.spikesNi spikes{n}];
    raster.neuronsNi = [raster.neuronsNi (i+1)*ones(size(spikes{n}))];
    i = i+1;
end
toc


rate = zeros(N,1500);
for n = 1:N
    for st = spikes{n}
        rate(n,round(st/10)) = rate(n,round(st/10))+1000;
    end
end

rate_av_Ne = mean(rate(1:Ne,:),1);
rate_av_Ni = mean(rate(Ne+1:end,:),1);
xs = 1:1500;
h = 10; %kernal bandwidth. determines the shape of the function
ys_Ne = zeros(1,1500);
ys_Ni = zeros(1,1500);
for i = 1:1500
    ys_Ne(i)=gaussian_kern_reg(xs(i),xs,rate_av_Ne,h);
    ys_Ni(i)=gaussian_kern_reg(xs(i),xs,rate_av_Ni,h);
end


figure
namefig = ['spont=' num2str(spont) ', param2=' num2str(param2) '.fig'];% ', stim =' num2str(stim_fr) ', stim_time=' num2str(input_time) ', param1=' num2str(param1)
% Rasterplot
% x = freq_list;
subplot(2,1,1)
xlabel('time (s)')
% ylabel('IPI (ms)')
% ylabel('Repetition rate (Hz)')
area([5000 5000 5000+input_time/step 5000+input_time/step],[0 120 120 0],'LineStyle','none','FaceColor',[.85 .85 1]);
hold on
plot(raster.spikesNe,raster.neuronsNe,'b.','MarkerSize',9);
plot(raster.spikesNi,raster.neuronsNi,'r.','MarkerSize',9);
% axis([0 stimulus_duration+POSTstimulus_duration 0 length(x)*nreps+1])

% Rasterplot end

subplot(2,1,2)
plot(ys_Ne);
hold on
plot(ys_Ni);
savefig(namefig)

% figure
% for i = 1:10
% subplot(10,1,i)
% plot(V(i,:))
% end

test = 1;



