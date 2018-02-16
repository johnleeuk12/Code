function Out = Conductance_LIF()

% clear all

%% Parameters
global ICI_list
ICI_list = [125 83.3333 62.5 50 41.6667 35.7143 31.25 27.7778 25 22.7273 20.8333]; % [20.8333 25 31.25 41.6667 62.5 125];%[20 30 40 50 60 80 100]; %[250 125 83.3333 62.5 50 41.6667 35.7143 31.25 27.7778 25 22.7273 20.8333];
IE_delay = 5; %ms
E_strength = 4.5; %:0.3:4.8; %1.5:0.3:6.; % in nS 3.5 for Sync+
I_strength = 8.5; %8.5; %:0.1:1.9;
% f_DE = 0.7;
% e = 2;%4.5;
% i = 2;%8.5;
n = 0;
%analysis of parameter space
% for e = 4.5 %3:5 %E_strength
%     for i = 4.5*1.7 %5:7%I_strength
%% global puretone
global pureT
pureT = 0;

global tau_pE tau_pI
tau_pE = 0.15;
tau_pI = 0.10;
pp = 0;
% figure
% f_DI = 0.9;
% f_DE = 0.6;
global kernel_time_constant
kernel_time_constant = 0.005;

% for tau_pE = 0.15


for f_DI = 1.2
    %         for kernel_time_constant = [0.005 0.010 0.020 0.040]
    %     for  f_DE = 0.7
    for f_DE = 0.9
        %             for tau_pI = 0.10:0.01:0.2
        tic
        pp = pp+1;
        disp(pp)
        output = {};
        output.raster.stim = [];
        output.raster.rep = [];
        output.raster.spikes = [];
        output.spiketime = {};
        for f = 1:length(ICI_list)
            
            out = run_model(IE_delay,E_strength,I_strength,f,f_DE,f_DI);
            output.raster.stim = [output.raster.stim out.raster.stim];
            output.raster.rep = [output.raster.rep out.raster.rep];
            output.raster.spikes = [output.raster.spikes out.raster.spikes];
            output.spiketime{f} = out.raster.spikes;
            output.VS(f) = out.VS;
%             output.VS_pop(f,:) = out.vector;
            output.mean_discharge_rate.mean(f) = out.discharge_rate.mean;
            output.mean_discharge_rate.error(f)  = out.discharge_rate.error;
            output.rate{f} = out.rate;
            output.rate_brut{f} = out.rate_brut;
            output.adaptation.E{f} = out.E_strength;
            output.adaptation.I{f} = out.I_strength;
            output.adaptation.E_I{f} = out.I_strength-out.E_strength;
            output.adaptation.E_I{f} = output.adaptation.E_I{f}(1:end-1);
            output.spikes_per_click{f} = out.spikes_per_click;
            output.time_period{f} = out.time_period;
            output.net_positivePclick{f} = out.net_positivePclick;
            output.net_negativePclick{f} = out.net_negativePclick;
            output.Fanofactor(f) = out.Fanofactor;
            output.var_ISI(f) = out.var_ISI;
            %                     for t = 1:length(output.adaptation.E_I{f})
            %                         output.product{f}(t) = output.adaptation.E_I{f}(t)*output.time_period{f}.mean(t);
            %                     end
        end
        %         figure
        %         plot(output.spikes_per_click{1}.xaxis,output.product{1})
        %         hold on
        %         plot(output.spikes_per_click{2}.xaxis,output.product{2})
        %         plot(output.spikes_per_click{3}.xaxis,output.product{3})
        
        
        %                 plot(output.rate{1},'linewidth', 2.0)
        %                 title(['f_DE = ' num2str(f_DE) ', f_DI = ' num2str(f_DI)])
        %                 pause(0.3);
        
        n = n+1;
        
        UnitInfo.List(n,1) = f_DE;
        UnitInfo.List(n,2) = f_DI;
        UnitInfo.List(n,3) = tau_pE;
        UnitInfo.List(n,4) = tau_pI;
        UnitInfo.List(n,5) = E_strength;
        UnitInfo.List(n,6) = I_strength;
        UnitInfo.Info(n).Output = output;
        %Analysis
        counter1 = 0;
        counter2 = 0;
        UnitInfo.Info(n).Sync = 0;
        UnitInfo.Info(n).Significant_rate = 0;
        for f = 1:length(ICI_list)
            if output.VS(f) >0.1
                counter1 = counter1 + 1;
            else
                counter1 = 0;
            end
            if counter1 == 3
                UnitInfo.Info(n).Sync = 1; %Sync
            end
            if output.mean_discharge_rate.mean(f) - 2*output.mean_discharge_rate.error(f)>0
                counter2 = counter2 + 1;
            else
                counter2 = 0;
            end
            if counter2 == 3
                UnitInfo.Info(n).Significant_rate = 1;% significant rate response
            end
        end
        %         if length(ICI_list) == 12
        %             new_ICI_list = ICI_list; %(1:end-1);
        %             new_all_mean_rate_stim = output.mean_discharge_rate.mean(1:end-1);
        %             %     else
        %             %         new_ICI_list = ICI_list(4:end);
        %             %         new_all_mean_rate_stim = all_mean_rate_stim(4:end);
        %         end
        new_ICI_list = ICI_list;
        new_all_mean_rate_stim = output.mean_discharge_rate.mean; %(1:end-1);
        
        [RHO,PVAL] = corr(new_ICI_list.',new_all_mean_rate_stim.','Type','Spearman')
        UnitInfo.Info(n).Rho = RHO;
        UnitInfo.Info(n).pval = PVAL;
        if PVAL <0.05
            if RHO > 0
                UnitInfo.Info(n).Positive = -1; % negative Sync for some reason
            elseif RHO < 0
                UnitInfo.Info(n).Positive = +1;
            end
        else
            UnitInfo.Info(n).Positive = 0;
        end
        %         toc
        %     end
        % end
        
        
        % figure
        % nreps = 10;
        % plot(output.raster.spikes,nreps*(output.raster.stim-1)+output.raster.rep,'k.','MarkerSize',9);
        
        test = 1 ;
        Hz_list = [];
        for i = 1:length(ICI_list)
            Hz_list = [Hz_list round(1000/ICI_list(i))];
        end
        
        save('model_datasyncP.mat')
        
        test;
        % figure
        % plot(UnitInfo.Info.Output.rate{1},'linewidth', 2.0)
        % pause
        
        %% plotting figures
        
        cmapp = [[0.1 0.7 0.1]; [0.9 0.6 0.1]; [0 0 0]   ];%    [0.26 0.5 0.9]   ]; %; [0.9 0.3 0.26]; ];
        
        cmap = colormap(jet(length(ICI_list)+1));
        figure
        
        for p = n
            norm_mean = [];
            subplot(2,2,[1 3])
            
            hold off
            for n =1:2:length(UnitInfo.Info(p).Output.rate)
                plot(UnitInfo.Info(p).Output.rate{n}, 'linewidth', 1.7,'color',cmap(n+1,:),'DisplayName', ...
                    [num2str(Hz_list(n)) 'Hz'])
                hold on
                %         norm_mean = [norm_mean UnitInfo.Info(p).Output.mean_discharge_rate.mean(n)/Hz_list(n)];
            end
                axis([300,1200,0,80]);
            legend('show')
            title(['// f_DE = ' num2str(UnitInfo.List(p,1))...
                '// f_DI = ' num2str(UnitInfo.List(p,2)) ...
                '// tau_pE = ' num2str(UnitInfo.List(p,3)) ...
                '// tau_pI = ' num2str(UnitInfo.List(p,4))])
            set(gca, 'FontSize', 16)
            hold on
            subplot(2,2,2)
            shadedErrorBar(Hz_list,UnitInfo.Info(p).Output.mean_discharge_rate.mean,UnitInfo.Info(p).Output.mean_discharge_rate.error,{'--','Color',cmapp(2,:)})
            %     errorbar(Hz_list,UnitInfo.Info(p).Output.mean_discharge_rate.mean,UnitInfo.Info(p).Output.mean_discharge_rate.error)
            
            
            hold on
            axis([0,50,-5,50]);
            subplot(2,2,4)
            
            plot(Hz_list,UnitInfo.Info(p).Output.VS,'linewidth',2.0);
            set(gca, 'FontSize', 16)
            axis([0,50,0.5,1]);
            %                 hold on
            %                 pause
            pause(0.1)
        end
        %             end
        %         end
    end
end
% test = 1;
% pause


% save('modeldata2.mat',UnitInfo)
% save('modeldata_new2.mat','UnitInfo')
% load('modeldata2.mat')
% for n = 1:length(UnitInfo.List)
%     if UnitInfo.Info(n).Sync == 1
%         scatter(UnitInfo.List(n,1),UnitInfo.List(n,2),[],'b')
%     else
%         scatter(UnitInfo.List(n,1),UnitInfo.List(n,2),[],'r')
%     end
%     hold on
% end

% negative = find([UnitInfo.Info.Positive]==-1);
% ICI_list= [20 30 40 50 60 80 100];


test = 1;


function out = run_model(IE_delay,E_strength,I_strength,f,f_DE,f_DI)

global ICI_list
global pureT
ICI = ICI_list(f);
global kernel_time_constant
% IE_ratio = I_strength/E_strength;
kernel_time_constant=.005;  %time constant of 5 ms
% jitter_magnitude=1; % Temporal jitter
step=.0001; %.1 ms duration  (temporal increment for running simulation)
stimulus_duration=0.5;  %half second
PREstimulus_duration=0.5;  %half second
POSTstimulus_duration=0.5;  %half second (0.1 second is included in stimulus)
latency_time=0.01;
latency=length(0:step:latency_time); %10 ms latency (auditory nerve to auditory cortex)
stepfn = zeros(1,15000);
stepfn(1,1:5000) = 0.01;
tau_d = 2*1e-3;
tau_r = 0.4*1e-3;
tau_l = 1e-3;
% tau_m(1) = 0.01;
% tau_m(2) = 0.01;
spikes_pooled=[];
freq = 1000./ICI;

%Alpha kernel
t=0:step:(kernel_time_constant*10);
kernel=t.*exp(-t/kernel_time_constant);
kernel=1e-9*kernel/max(kernel); %amplitude of 1 nS


%exponential kernel
% kernel = zeros(1,501);
% for t = 1:501
%     if t>round(tau_l/step)
%         
%         kernel(1,t) = (exp(-(t*step-tau_l)/tau_d)-exp(-(t*step-tau_l)/tau_r));
%     end
% end
% 
% kernel = kernel*1e-9;


input=zeros(size(0:step:(POSTstimulus_duration+stimulus_duration)));
stimulus_input_length=length(0:step:(stimulus_duration));
ipi=round(1/(freq*step)); %ipi=interpulse interval
freq2=1/(step*ipi);

puretone_input = conv(stepfn,kernel);


P_0 = 1;
E_str(1) = E_strength;
I_str(1) = I_strength;
E_strength_mean = [];
I_strength_mean = [];



% %adaptation parameters
%
global tau_pE tau_pI
% tau_pE = 0.12;
% tau_pI = 0.05;
% tau_pE = adapt.tau_pE;
% tau_pI = adapt.tau_pI;
% f_DI = adapt.f_DI;

% if f_DE ~= 0
%     f_DI = 1.6-f_DE;
% else    % No adaptation
%     f_DE = 1;
%     f_DI = 1;
% end

%% Model (Conductance, Spikes etc)

nb_rep = 30;
rate_total= [];
Ge_total = [];
Gi_total = [];
Net_excit_total = [];
raster.stim=[];  raster.rep=[];  raster.spikes=[];
vector_pop = [];


for r = 1:nb_rep
    E_input=input;
    I_input=input;
    
    for j=1:10  %10 jitter excitatory and inhibitory inputs
        % for i=1:ipi:(stimulus_input_length-(length(kernel)/2))
        p = 0;
        P_relE(1) = P_0;
        P_relI(1) = P_0;
        for i=1:ipi:(stimulus_input_length)
            p = p+1; %Click number (starts with 1)
            jitter=round(randn(1)/(1000*step)); %1 ms jitter
            
            if (i+jitter)<1 || (i+jitter)>(length(input)-length(kernel))
                jitter=1;
            end
            t0 = i+jitter;
            if i == 1
                P_relE(1:t0) = P_0;
                P_relI(1:t0) = P_0;
            end
            for t = t0:t0+2*ipi-1
                P_relE(t+1) = P_0 + (f_DE*P_relE(t0)-P_0)*exp(-(t+1-t0)*step/tau_pE);
                P_relI(t+1) = P_0 + (f_DI*P_relI(t0)-P_0)*exp(-(t+1-t0)*step/tau_pI);
                %                 if P_relE(t) <0.65;
                %                     ppp = 1;
                %                                         pause
                %                 end
            end
            
            E_str(p) = P_relE(t0)*E_strength;
            I_str(p) = P_relI(t0)*I_strength;
            E_input((latency+t0):(latency+t0+length(kernel)-1))=E_input((latency+t0):(latency+t0+length(kernel)-1))+ kernel*E_str(p); % exponential  +tau_m(1)/(tau_d-tau_r)*kernel*E_str(p);
            %             jitter=round(randn(1)/(1000*step)); %1 ms jitter
            %
            %             if (t0)<1 || (t0)>(length(input)-length(kernel))
            %                 jitter=0;
            %             end
            I_input((latency+t0):(latency+t0+length(kernel)-1))=I_input((latency+t0):(latency+t0+length(kernel)-1))+ kernel*I_str(p); %exponential +tau_m(2)/(tau_d-tau_r)*kernel*I_str(p);
            
            %             plot(P_relE)
            %             hold on
        end
        E_strength_mean = [E_strength_mean ; E_str];
        I_strength_mean = [I_strength_mean ; I_str];
    end
    
    %   %
    delay=round(abs(IE_delay)/(1000*step));  %delay in steps
    if pureT == 1
        E_input(1:length(P_relE)) = puretone_input(1:length(P_relE)).*P_relE;
        I_input(1:length(P_relI)) = puretone_input(1:length(P_relI)).*P_relI;
        Ge = E_input*E_strength;
        Gi=[zeros(1,delay) I_input(1:(length(I_input)-delay))*I_strength];
        
    else
        if IE_delay>=0
            Ge=E_input;
            Gi=[zeros(1,delay) I_input(1:(length(I_input)-delay))];
        elseif IE_delay<0
            Gi=I_input;
            Ge=[zeros(1,delay) E_input(1:(length(E_input)-delay))];
        end
    end
    
    
    %     if IE_delay>=0
    %         Ge=E_input;
    %         Gi=[zeros(1,delay) I_input(1:(length(I_input)-delay))];
    %     elseif IE_delay<0
    %         Gi=I_input;
    %         Ge=[zeros(1,delay) E_input(1:(length(E_input)-delay))];
    %     end
    %
    %     if IE_delay>=0
    %         Ge=E_input*E_strength;
    %         Gi=[zeros(1,delay) I_input(1:(length(I_input)-delay))]*I_strength;
    %     elseif IE_delay<0
    %         Gi=I_input*I_strength;
    %         Ge=[zeros(1,delay) E_input(1:(length(E_input)-delay))]*E_strength;
    %     end
    
    %add pre  stim time of 500 ms
    Ge=[zeros(size(0:step:PREstimulus_duration)) Ge];
    Gi=[zeros(size(0:step:PREstimulus_duration)) Gi];
    Ge_total = [Ge_total; Ge];
    Gi_total = [Gi_total; Gi];
    Net_excit_total = [Net_excit_total; Ge-Gi];
    [spikes,V]=run_LIFmodel(Ge,Gi,f);
    
    %rate
    rate = zeros(1,1500);
    for st = spikes
        rate(1,round(st/(step*10))) = rate(1,round(st/(step*10)))+1;
    end
    rate_total = [rate_total ; rate*1000];
    
    %ISI
    spikes3 = spikes(find(spikes<PREstimulus_duration+0.5));
    spikes4 = [0 spikes3(1:end-1)];
    isi = spikes3-spikes4;
    isi = isi(2:end);
    %pool spikes
    spikes=spikes-PREstimulus_duration;
    spikes_pooled=[spikes_pooled spikes];
%     spikes_pooled_for_VS = spikes(find(spikes>0.05 & spikes<=(stimulus_duration+0.05)));
%     if ~isempty(spikes_pooled_for_VS)
%         total_spikes=length(spikes_pooled_for_VS);
%         x=0;
%         y=0;
%         if total_spikes>0
%             x=sum(cos(2*pi*(spikes_pooled_for_VS*freq2)));
%             y=sum(sin(2*pi*(spikes_pooled_for_VS*freq2)));
%         end
%         if total_spikes==0
%             vector=0;
%         else
%             vector=sqrt(x^2+y^2)/total_spikes;
%         end
%         rayleigh=2*total_spikes*vector^2;
%     else
%         vector=0;
%         rayleigh=0;
%     end
%     
%     if rayleigh<13.8
%         vector=0;
%     end
%     
%     vector_pop = [vector_pop;vector];
    %     spike_distribution(f,r)=length(spikes(find(spikes>0 & spikes<=(stimulus_duration+0.1))))/(stimulus_duration+0.1);
    %     spont_distribution(r+nreps*(f-1))=length(spikes(find(spikes>-PREstimulus_duration & spikes<0)))/PREstimulus_duration;
    %
    raster.stim=[raster.stim f*ones(size(spikes))];
    raster.rep=[raster.rep r*ones(size(spikes))];
    raster.spikes=[raster.spikes spikes];
end


out.raster = raster;
% out.vector = vector_pop;
%% Analysis

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
        vector=0;
    else
        vector=sqrt(x^2+y^2)/total_spikes;
    end
    rayleigh=2*total_spikes*vector^2;
else
    vector=0;
    rayleigh(f)=0;
end

if rayleigh<13.8
    vector=0;
end

%VS pop analysis. 



%  average rate

PRE = PREstimulus_duration*1000;
POST = POSTstimulus_duration*1000;
STIM = stimulus_duration*1000;
total_time = PRE+POST+STIM;
rate_av = mean(rate_total,1);
spont_rate = mean2(rate_total(:,1:PRE));
discharge_rate.mean = mean2(rate_total(:,PRE+1:PRE+STIM+100))-spont_rate;
SEM = std2(rate_total(:,PRE+1:PRE+STIM+100))/sqrt(nb_rep*(STIM+100));
ts = tinv([0.025  0.975],nb_rep*(STIM+100)-1);      % T-Score 95%
discharge_rate.error = ts(2)*SEM;

%

% CI = mean(x) + ts*SEM;                      % Confidence Intervals



xs = 1:total_time;
h = 10; %kernal bandwidth. determines the shape of the function
for i = 1:total_time
    ys(i)=gaussian_kern_reg(xs(i),xs,rate_av,h);
end

% spikes per click
spikes_per_click = {};
p = floor(STIM/ICI);
for q = 1:p
    clicktime = round((q-1)*ICI)+15+PRE; %plus input latency + kernel peak
    spikes_per_click.mean(q) = mean2(rate_total(:,clicktime - 10 : clicktime + 10));
    spikes_per_click.std(q) = std2(rate_total(:,clicktime - 10 : clicktime + 10))/sqrt(20*nb_rep);
    spikes_per_click.xaxis(q) = clicktime;
end


% Fanofactor
SpikeCount = sum(rate_total,2)/1000.;
out.Fanofactor = std(SpikeCount)^2/mean(SpikeCount);


%Var ISI
out.var_ISI = std(isi);

out.rate_brut = rate_av;
out.rate = ys;
out.VS = vector;
out.spikes_per_click = spikes_per_click ;
out.discharge_rate = discharge_rate;

% evolution of E_strength and I_strength

out.E_strength = mean(E_strength_mean,1);
out.I_strength = mean(I_strength_mean,1);
for ind = 1:length(out.I_strength)
    out.IE_ratio(ind) = out.I_strength(ind)/out.E_strength(ind);
end

% time period when E > I
% time resolution is 10 times higher than rate.
% net excitation per click (as area under curve)

times = [];
total_time_period = [];
total_net_positivePclick = zeros(nb_rep,p);
total_net_negativePclick = zeros(nb_rep,p);
for k = 1:nb_rep
    for q = 1:p
        clicktime0 = round((q-2)*ICI)+PRE;
        clicktime1 = round((q-1)*ICI)+PRE;
        clicktime2=  round((q)*ICI)+PRE;
        times = [];
        for t = clicktime1*10 : clicktime2*10
            if Net_excit_total(k,t) > 0
                times(t) = t;
                total_net_positivePclick(k,q) = total_net_positivePclick(k,q) + Net_excit_total(k,t);
            end
        end
        for t = clicktime0*10 : clicktime1*10
            if Net_excit_total(k,t) < 0
                total_net_negativePclick(k,q) = total_net_negativePclick(k,q) +  Net_excit_total(k,t);
            end
        end
        
        if ~isempty(times)
            total_time_period(k,q) = round(max(times)-min(times(times>0)));
        else
            total_time_period(k,q) = 0;
        end
        
        
    end
end
time_period.mean = mean(total_time_period,1);
time_period.std = std(total_time_period,1);
net_positivePclick.mean = mean(total_net_positivePclick,1);
net_positivePclick.std = std(total_net_positivePclick,1);
net_negativePclick.mean = mean(total_net_negativePclick,1);
net_negativePclick.std = std(total_net_negativePclick,1);

out.time_period = time_period;
out.net_positivePclick = net_positivePclick;
out.net_negativePclick = net_negativePclick;

% net excitation per click (as area under curve)


function [spikes,V]=run_LIFmodel(Ge,Gi,f)
spikes=[]; V=[]; t=1; i=1;
step=.0001; %.1 ms duration  (temporal increment for running simulation)
C=0.25*1e-9; %0.25 nF, 10 ms time constant
Grest=25*1e-9; %25 nS
Erest=-0.065; %-65 mV
Ee=0; % 0 mV
Ei=-0.085;  %-85 mV
Ek = -0.075;
% Ereset = -0.065; %     =   -0.059;
% Vth = -0.045;   %      =   -0.052;
% V = zeros(1,length(Ge));
V(1)=Erest;

noise_magnitude=4*1e-8; %default noise level in conductance

%avoid negative conductances
Ge=Ge+noise_magnitude*randn(1,length(Ge));
Gi=Gi+noise_magnitude*randn(1,length(Gi));
Ge(find(Ge<0))=0;
Gi(find(Gi<0))=0;
sigma = 0.01    ;

%spike rate adaptation
Gsra = zeros(1,length(Ge));
tau_sra = 0.1; %100ms
delta_sra = 1*1e-9;
f_sra = 0.998;

%synaptic depression
% global ICI_list
% global tau_pE tau_pI

% ICI = ICI_list(f);
% freq = 1000./ICI;
% ipi=round(1/(freq*step)); %ipi=interpulse interval

while(t<length(Ge))
    V(t+1)=(-step*( Ge(t)*(V(t)-Ee) + Gi(t)*(V(t)-Ei) + Grest*(V(t)-Erest))/C)+V(t) + sigma*randn*sqrt(step); %+Gsra(t)*(V(t)-Ek))/C)
    Gsra(t+1) = Gsra(t)*f_sra; % + sigma*randn*sqrt(step);
    if V(t+1)>(Erest+0.020) %20 mV above Erest %artificial threshold
        V(t+1)=0.050; %spike to 50 mV
        spikes(i)=step*(t+1);
        Gsra(t+1) = Gsra(t+1) +delta_sra;
        i=i+1;
        t=t+1;
        V(t+1)=Erest;
        Gsra(t+1) = Gsra(t)*f_sra;
        
    end
    t =t+1;
end













% % % % expontential model, testing
% for t = 1:length(Ge)
%     
%     
%     if ~isempty(spikes)
%         if t-spikes(end)/step <2
%             V(t+1) = Ereset + sigma*randn*sqrt(step);
%         else
%             V(t+1)=(-step*( Ge(t)*(V(t)-Ee) + Gi(t)*(V(t)-Ei) + Grest*(V(t)-Erest)+Gsra(t)*(V(t)-Ek))/C)+V(t) + sigma*randn*sqrt(step);
%         end
%         Gsra(t+1) = Gsra(t)*f_sra; % + sigma*randn*sqrt(step);
%     else
%         V(t+1)=(-step*( Ge(t)*(V(t)-Ee) + Gi(t)*(V(t)-Ei) + Grest*(V(t)-Erest)+Gsra(t)*(V(t)-Ek))/C)+V(t) + sigma*randn*sqrt(step);
%         
%     end
%     Gsra(t+1) = Gsra(t)*f_sra;% + sigma*randn*sqrt(step);
%     
%     if V(t+1)>Vth  %20 mV above Erest %artificial threshold
%         if ~isempty(spikes)
%             if t > spikes(end)/step+(mod(1,2)+1)*1e-3/step %refractory period
%                 V(t+1)=0.050; %spike to 50 mV
%                 spikes(i)=step*(t+1);
%                 i=i+1;
%                 %         t=t+1;
%                 %         V(t+1)=Ereset;
%                 Gsra(t+1) = Gsra(t+1)+ delta_sra ;
%             end
%         else
%             V(t+1)=0.050; %spike to 50 mV
%             spikes(i)=step*(t+1);
%             i=i+1;
%             %         t=t+1;
%             %         V(t+1)=Ereset;
%             Gsra(t+1) = Gsra(t+1)+ delta_sra ;
%         end
%         
%     end
% end
% 
% test =1;

%% if need spike rate adaptation
% delta_sra = 20*1e-9; %40ns
% G_sra(1) = 0;
% f_sra = 0.998;
% Ek = -0.075;
% while(t<length(Ge))
%
%     V(t+1)=(-step*( Ge(t)*(V(t)-Ee) + Gi(t)*(V(t)-Ei) + Grest*(V(t)-Erest))/C)+V(t) + sigma*randn*sqrt(step);
%     %     V(t+1)=(-step*(Ge(t)(V(t)-Ee)+Gi(t)(V(t)-Ei)+Grest*(V(t)-Erest))/C)+V(t)...
%     %     + sigma*randn*sqrt(step);- G_sra(t)*(V(t)-Ek) + G_sra(t)*(V(t)-Ek)
% %     G_sra(t+1) = G_sra(t)*f_sra;
%     if V(t+1)>(Erest+0.020) %20 mV above Erest %artificial threshold
%         V(t+1)=0.050; %spike to 50 mV
%         spikes(i)=step*(t+1);
% %         G_sra(t+1) = G_sra(t+1) +delta_sra;
%         i=i+1;
%         t=t+1;
%         V(t+1)=Erest;
% %         G_sra(t+1) = G_sra(t)*f_sra;
%
%     end
%
%     t=t+1;
%
% end

