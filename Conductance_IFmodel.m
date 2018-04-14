function Data = Conductance_IFmodel() %IE_delay,E_strength,IE_ratio)

%% Parameters
%E_strength: strength of excitation
%I_strength:strength of inhibition
%IE_ratio:ratio of Excitation strength to Inhibition strength
%IE_delay:delay in milliseconds between excitation and inhibition
clear all
tic
IE_delay  = 5.; % in ms
E_strength_range = 4.5; %:0.3:4.8; %1.5:0.3:6.; % in nS
IE_ratio_range = 1.5; %:0.1:1.9;
% ICI = 20.; %20 is 50hz, where we could observe  synaptic depression in some cases


global ICI_list
ICI_list= [20 25 30 40 60 100];% [3 5 7.5 10 12.5 15 20 25 30 35 40 45 50 55 60 65 70 75];

Data = {};
i = 1;
% global tau_pE tau_pI f_DE f_DI
% 
% tau_pE = 0.12;
% tau_pI = 0.05; 
% f_DE = 0.85; %need to fit to data, decay function
% f_DI = 0.7;


figure
cmap = colormap(hsv(length(ICI_list)));
variables = 0.05:0.01:0.2;


% for var = variables
% for var1 = 0.75 % 0.6:0.05:0.9
%     %     tau_pE = var;
%     %     tau_pI = var;
%     f_DE = var1 ;
%     for var2 = 0.85% 0.6:0.05:0.9
%         f_DI = var2 ;
%         %         subplot(7,7,i)
%         for E_strength = E_strength_range
%             tic
%             for j = 1:length(IE_ratio_range)
%                 for k = 1:length(ICI_list)
%                     
%                     Data.IE_ratio{j} = IE_ratio_range(j);
%                     Data.E_strength{i} = E_strength;
%                     out = run_model(IE_delay,E_strength,IE_ratio_range(j),k);
%                     Data.rate{i,j,k} = out.rate;
%                     Data.VS{i,j,k} = out.VS;
%                     subplot(2,1,2)
%                     plot(out.E_strength,out.IE_ratio,'color',cmap(k,:)) %,'LineWidth',10,'Marker','o')
%                     axis([1.5,4.6,1.5,3]);
%                     hold on
% %                     set(gca,'FontSize',20)
%                     
%                     subplot(2,1,1)
%                     plot(out.rate,'color',cmap(k,:),'LineWidth',1.7,'DisplayName', ...
%                         [num2str(ceil(1000/ICI_list(k))) 'Hz'] )
%                     axis([0,1500,0,60]);
%                     hold on
% %                     set(gca,'FontSize',20)
%                     
%                     %                 disp(var)
%                     %                                 pause(0.3)
%                 end
%             end
%             
%         end
% %         set(gca,'FontSize',20)
% %         legend('show')
%         %         title([num2str(var1) ',' num2str(var2)])
%         i = i+1;
%         %     pause(0.3)
%         toc
%     end
% end

%
% save('Data')
toc

% load('Data.mat')
% % 
% ICI = 5;
% data_plot = [4.5,1.8,ICI]; % [E strength,IE_ratio]
% indice.x = find(E_strength_range == data_plot(1));
% indice.y = find(IE_ratio_range == data_plot(2));
% indice.z = find(ICI_list == data_plot(3));
% 
% % % 
% % plot_VS = zeros(1,18);
% % for i = 1:length(ICI_list)
% % plot_VS(i) = Data.VS{indice.x,indice.y,i};
% % end
% % 
% % plot(ICI_list,plot_VS)
% % 
% Data.time = -499:1000;
% Data.ICI = ICI;
% Data.ICI_vector = 0:ICI:500;
% 
% figure
% xlabel('ms')
% ylabel('rate, spikes per s')
% plot(Data.time,Data.rate{indice.x,indice.y,indice.z})
% hold on
% scatter(Data.ICI_vector,ones(size(Data.ICI_vector))*1)
% grid on 


% two neurons

adapt1.tau_pE = 0.15;
adapt1.tau_pI = 0.1;
adapt1.f_DI = 0.6;
adapt1.f_DE = 0.9;

adapt2.tau_pE = 0.15;
adapt2.tau_pI = 0.1;
adapt2.f_DI = 0.9;
adapt2.f_DE = 0.6;

adapt3.tau_pE = 0.15;
adapt3.tau_pI = 0.1;
adapt3.f_DI = 1.0;
adapt3.f_DE = 0.5;

cmap = colormap(hsv(length(ICI_list)));
for k = [1 2 3 6] % 1:length(ICI_list)
    out1 = run_model(IE_delay,4.5,1.7,k,adapt1);
    out2 = run_model(IE_delay,4.5,1.7,k,adapt2);
    out3 = run_model(IE_delay,4.5,1.7,k,adapt3);
    rate_decal = zeros(1,1500);
    rate_decal(1,31:1480) = out3.rate_brut(1,1:1450);
    rate = out2.rate_brut ;%+ rate_decal;
    rate2 = out3.rate_brut;% + rate_decal;
    rate3 = out1.rate_brut;
%     rate3 = 1.3*rate - rate2;
    xs = 1:1500;
    h = 15; %kernal bandwidth. determines the shape of the function
%     for i = 1:1500
%         ys1(i)=gaussian_kern_reg(xs(i),xs,rate,h);
%         ys2(i)=gaussian_kern_reg(xs(i),xs,rate3,h);
%         ys3(i)=gaussian_kern_reg(xs(i),xs,rate2,h);
%     end
%     
%     plot(out1.rate_brut,'color',cmap(k,:),'LineWidth',1.7,'DisplayName', ...
%         [num2str(ceil(1000/ICI_list(k))) 'Hz'] );
%     hold on
%     
%     plot(out2.rate_brut,'color',cmap(k+2,:),'LineWidth',1.7,'DisplayName', ...
%         [num2str(ceil(1000/ICI_list(k))) 'Hz'] );
%     
%     hold on
%     plot(out3.rate_brut,'color',cmap(k+4,:),'LineWidth',1.7,'DisplayName', ...
%         [num2str(ceil(1000/ICI_list(k))) 'Hz'] );
    subplot(2,2,1)
    plot(out1.rate,'color',cmap(k,:),'LineWidth',1.7,'DisplayName', ...
        [num2str(ceil(1000/ICI_list(k))) 'Hz'] );
    hold on
    subplot(2,2,2)
    plot(out3.rate,'color',cmap(k,:),'LineWidth',1.7,'DisplayName', ...
        [num2str(ceil(1000/ICI_list(k))) 'Hz'] );
    hold on
    subplot(2,2,4)
    plot(out2.rate,'color',cmap(k,:),'LineWidth',1.7,'DisplayName', ...
        [num2str(ceil(1000/ICI_list(k))) 'Hz'] );
    hold on
        subplot(2,2,3)
%     plot(ys2,'color',cmap(k,:),'LineWidth',1.7,'DisplayName', ...
%         [num2str(ceil(1000/ICI_list(k))) 'Hz'] );
%     hold on
end

test = 1

function out = run_model(IE_delay,E_strength,IE_ratio,f,adapt)

global ICI_list
ICI = ICI_list(f);

I_strength=IE_ratio*E_strength;
kernel_time_constant=.005;  %time constant of 5 ms
% jitter_magnitude=1; % Temporal jitter
step=.0001; %.1 ms duration  (temporal increment for running simulation)
stimulus_duration=0.5;  %half second
PREstimulus_duration=0.5;  %half second
POSTstimulus_duration=0.5;  %half second (0.1 second is included in stimulus)
latency_time=0.01;
latency=length(0:step:latency_time); %10 ms latency (auditory nerve to auditory cortex)

spikes_pooled=[];
freq = 1000./ICI;
t=0:step:(kernel_time_constant*10);
kernel=t.*exp(-t/kernel_time_constant);
kernel=1e-9*kernel/max(kernel); %amplitude of 1 nS
input=zeros(size(0:step:(POSTstimulus_duration+stimulus_duration)));
stimulus_input_length=length(0:step:(stimulus_duration));
ipi=round(1/(freq*step)); %ipi=interpulse interval
freq2=1/(step*ipi);

% tau_pE = 0.02;
% tau_pI = 0.1;

P_0 = 1;

E_str(1) = E_strength;
I_str(1) = I_strength;
E_strength_mean = [];
I_strength_mean = [];

% global tau_pE tau_pI f_DI f_DE
%adaptation parameters

tau_pE = adapt.tau_pE;
tau_pI = adapt.tau_pI;
f_DI = adapt.f_DI;
f_DE = adapt.f_DE;


%% Model (Conductance, Spikes etc)

nb_rep = 30;
rate_total = [];




for r = 1:nb_rep
    E_input=input;
    I_input=input;
    
    for j=1:10  %10 jitter excitatory and inhibitory inputs
        % for i=1:ipi:(stimulus_input_length-(length(kernel)/2))
        p = 1;
        P_relE(1) = P_0;
        P_relI(1) = P_0;
        for i=1:ipi:(stimulus_input_length-250)
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
            E_input((latency+t0):(latency+t0+length(kernel)-1))=E_input((latency+t0):(latency+t0+length(kernel)-1))+kernel*E_str(p);
            %             jitter=round(randn(1)/(1000*step)); %1 ms jitter
            %
%             if (t0)<1 || (t0)>(length(input)-length(kernel))
%                 jitter=0;
%             end
            I_input((latency+t0):(latency+t0+length(kernel)-1))=I_input((latency+t0):(latency+t0+length(kernel)-1))+kernel*I_str(p);
            p = p+1;
%             plot(P_relE)
%             hold on
        end
        E_strength_mean = [E_strength_mean ; E_str];
        I_strength_mean = [I_strength_mean ; I_str];
    end
    
%     
    delay=round(abs(IE_delay)/(1000*step));  %delay in steps
    if IE_delay>=0
        Ge=E_input;
        Gi=[zeros(1,delay) I_input(1:(length(I_input)-delay))];
    elseif IE_delay<0
        Gi=I_input;
        Ge=[zeros(1,delay) E_input(1:(length(E_input)-delay))];
    end
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
    test = cumsum(Ge-Gi);
    
    [spikes,V]=run_LIFmodel(Ge,Gi,f);
    
    %rate
    rate = zeros(1,1500);
    for st = spikes
        rate(1,round(st/(step*10))) = rate(1,round(st/(step*10)))+1;
    end
    rate_total = [rate_total ; rate*1000];
    %pool spikes
    spikes=spikes-PREstimulus_duration;
    spikes_pooled=[spikes_pooled spikes];
    
    %     spike_distribution(f,r)=length(spikes(find(spikes>0 & spikes<=(stimulus_duration+0.1))))/(stimulus_duration+0.1);
    %     spont_distribution(r+nreps*(f-1))=length(spikes(find(spikes>-PREstimulus_duration & spikes<0)))/PREstimulus_duration;
    %
    %     raster.stim=[raster.stim f*ones(size(spikes))];
    %     raster.rep=[raster.rep r*ones(size(spikes))];
    %     raster.spikes=[raster.spikes spikes];
end

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


%  average rate
rate_av = mean(rate_total,1);
xs = 1:1500;
h = 15; %kernal bandwidth. determines the shape of the function 
for i = 1:1500
    ys(i)=gaussian_kern_reg(xs(i),xs,rate_av,h);
end
% plot(ys)
out.rate_brut = rate_av;
out.rate = ys;
out.VS = vector;

% evolution of E_strength and I_strength

out.E_strength = mean(E_strength_mean,1);
out.I_strength = mean(I_strength_mean,1);
for ind = 1:length(out.I_strength)
    out.IE_ratio(ind) = out.I_strength(ind)/out.E_strength(ind);
end


function [spikes,V]=run_LIFmodel(Ge,Gi,f)
spikes=[]; V=[]; t=1; i=1; P_relE = [];
step=.0001; %.1 ms duration  (temporal increment for running simulation)
C=0.25*1e-9; %0.25 nF, 10 ms time constant
Grest=25*1e-9; %25 nS
Erest=-0.065; %-65 mV
Ee=0; % 0 mV
Ei=-0.085;  %-85 mV
Ek = -0.075;

V(1)=Erest;
G_sra(1) = 0;
noise_magnitude=4e-8; %default noise level in conductance
delta_sra = 20*1e-9; %40ns

%avoid negative conductances
Ge=Ge+noise_magnitude*randn(1,length(Ge));
Gi=Gi+noise_magnitude*randn(1,length(Gi));
Ge(find(Ge<0))=0;
Gi(find(Gi<0))=0;
sigma = 0.01    ;

%synaptic depression
global ICI_list
global tau_pE tau_pI

ICI = ICI_list(f);
freq = 1000./ICI;
ipi=round(1/(freq*step)); %ipi=interpulse interval

% f_sra = 0.998;
while(t<length(Ge))
    
    V(t+1)=(-step*( Ge(t)*(V(t)-Ee) + Gi(t)*(V(t)-Ei) + Grest*(V(t)-Erest))/C)+V(t) + sigma*randn*sqrt(step);
    %     V(t+1)=(-step*(Ge(t)*P_relE_1*(V(t)-Ee)+Gi(t)*P_relI_1*(V(t)-Ei)+Grest*(V(t)-Erest))/C)+V(t)
    %     + sigma*randn*sqrt(step);- G_sra(t)*(V(t)-Ek) + G_sra(t)*(V(t)-Ek)
    %     *P_rel(t)
%     G_sra(t+1) = G_sra(t)*f_sra;
    if V(t+1)>(Erest+0.020) %20 mV above Erest %artificial threshold
        V(t+1)=0.050; %spike to 50 mV
        spikes(i)=step*(t+1);
%         G_sra(t+1) = G_sra(t+1) +delta_sra;
        i=i+1;
        t=t+1;
        V(t+1)=Erest;
%         G_sra(t+1) = G_sra(t)*f_sra;
        
    end
    
    t=t+1;
    
end
%
%
% plot(G_sra)
t  = 1;%






