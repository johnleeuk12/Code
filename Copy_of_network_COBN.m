function copy_of_network_COBN()

%Parameters, variables
% CX2,EIbalance,VIP, AtoL2, BtoL2,VIPtoL2, L2toVIP, RecEx, RecIn
CX2 = 0.02; %layer2 network connectivity
EIbalance = 0.80;
VIP =0;% (initial inhibition strength)
AtoL2 = 0.0;
% BtoL2 = -1.5;
L1toIn = 0.;
VIPtoL2Ex = 0;
L1toVIP = 0;
% AtoVIP = 0 ;
BtoVIP = 0 ;
L2toVIP = 0;
RecEx = 0.6;
RecIn = -2.;
Posi = 0; % 1 for positive, 0 for negative
sim = 0;% 1 for indv neuron analysis, 0 for population analysis
iii = 2; % (number of simulations)

% for AtoL2 = 0.35
%     for Posi = 0
%         iii = 8;
%         run_network(CX2,EIbalance,VIP, AtoL2, L1toIn,VIPtoL2Ex,L1toVIP,BtoVIP, L2toVIP, RecEx, RecIn,Posi,sim,iii)
%         pause(0.1)
%     end
% end

% for AtoL2 = 0.
%     for Posi = 0
%         for sim = 1
%             for RecEx = 0.8
%                 for RecIn = -2.3
%                     iii = 7;
%                     run_network(CX2,EIbalance,VIP, AtoL2, L1toIn,VIPtoL2Ex,L1toVIP,BtoVIP, L2toVIP, RecEx, RecIn,Posi,sim,iii)
%                     pause(0.1)
%                 end
%             end
%         end
%     end
% end

for CX2 = 0.2
    for RecEx = 0.4
        RecIn = -2.3;
        run_network(CX2,EIbalance,VIP, AtoL2, L1toIn,VIPtoL2Ex,L1toVIP,BtoVIP, L2toVIP, RecEx, RecIn,Posi,sim,iii)
        pause(0.1)
    end
end


% for CX2 = 0:0.05:0.2
%     run_network(CX2,EIbalance,VIP, AtoL2, BtoL2,VIPtoL2Ex,AtoVIP,BtoVIP, L2toVIP, RecEx, RecIn,Posi,sim,iii)
%     pause(0.2)
% end


function run_network(CX2,EIbalance,VIP, AtoL2, L1toIn,VIPtoL2Ex,L1toVIP,BtoVIP, L2toVIP, RecEx, RecIn,Posi,sim,iii)
stimulus_duration=0.5;  %half second
PREstimulus_duration=0.5;  %half second
POSTstimulus_duration=0.5;  %half second (0.1 second is included in stimulus)

nb_rep = 4;
Nb_neurons.Layer1 = 50; %Nb of first layer neurons with A number of positive and B number of negative sync neurons
Nb_neurons.A = round(Nb_neurons.Layer1*0.5);
Nb_neurons.B = Nb_neurons.Layer1-Nb_neurons.A;
% Nb_neurons.NbSelf = 1;
NN = 1400;
Nb_neurons.Ex = NN*EIbalance;
Nb_neurons.In = NN*(1-EIbalance);
Nb_neurons.Layer2 = Nb_neurons.In + Nb_neurons.Ex; %Nb of second layer neurons

Nb_neurons.Layer3 = 50; %VIP inhibitory neuron
Nb_neurons.total = Nb_neurons.Layer1 + Nb_neurons.Layer2 + Nb_neurons.Layer3;



global tau_pE tau_pI IE_delay ICI_list stepp
stepp=.0001; %.1 ms duration  (temporal increment for running simulation)
% ICI_list = [125 62.5 41.66667 31.25 25 20.8333]; % 62.5 41.6667 31.25 25 20.8333];% [125 20.83333];%  ;%%[20.8333 25 31.25 41.6667 62.5 125]; %[20 30 40 70 100]; %[20 40 100];
ICI_list = 125 % [125 83.3333 62.5 50 41.6667 35.7143 31.25 27.7778 25 22.7273 20.8333];
tau_pE = 0.15;
% tau_pI = 0.1; %Sync+
tau_pI = 0.10; %sync-
IE_delay = 5; %ms

% W = ones(Nb_neurons.total,Nb_neurons.total)*0.6;
% W(1:20,:) = 0;
W = zeros(Nb_neurons.total,Nb_neurons.total);
W(Nb_neurons.Layer1+1:Nb_neurons.Layer1+Nb_neurons.Ex, 1:Nb_neurons.A) = AtoL2; %input from A to layer2
% W(Nb_neurons.Layer1+Nb_neurons.Ex+1:end-Nb_neurons.Layer3, 1:Nb_neurons.Layer1) = L1toIn; %input from B to layer 2
% W(Nb_neurons.Layer1+Nb_neurons.Layer2+1:end, 1:Nb_neurons.Layer1) = L1toVIP; 
% W(Nb_neurons.Layer1+Nb_neurons.Layer2+1:end, 1:Nb_neurons.A) = AtoVIP; %input from A to layer2
% W(Nb_neurons.Layer1+Nb_neurons.Layer2+1:end, Nb_neurons.A+1:Nb_neurons.Layer1) = BtoVIP; %input from B to layer 2

% W(Nb_neurons.Layer1+Nb_neurons.Layer2+1:end, 1:Nb_neurons.Layer1) = L1toVIP;% -2; %input from layer 1 to VIP
% W(Nb_neurons.Layer1+1:Nb_neurons.Layer1+Nb_neurons.Ex, Nb_neurons.Layer1+Nb_neurons.Layer2+1:end) =VIPtoL2Ex;% -3; % VIP to layer 2
% W(Nb_neurons.Layer1+Nb_neurons.Layer2+1:end, Nb_neurons.Layer1+1:Nb_neurons.total-Nb_neurons.Layer3) = L2toVIP; %-0.5; %layer 2 to 2VIP
% W(Nb_neurons.Layer1+1:Nb_neurons.total-Nb_neurons.Layer3, Nb_neurons.Layer1+1:Nb_neurons.total-Nb_neurons.Layer3) = rand(Nb_neurons.Layer2,Nb_neurons.Layer2)<CX2; %recurrent connections layer 2
% W(Nb_neurons.Layer1+1:end, Nb_neurons.Layer1+1:end) = W(Nb_neurons.Layer1+1:end, Nb_neurons.Layer1+1:end);
W(Nb_neurons.Layer1+1:Nb_neurons.total-Nb_neurons.Layer3, Nb_neurons.Layer1+1:Nb_neurons.Layer1+Nb_neurons.Ex) = RecEx;
W(Nb_neurons.Layer1+1:Nb_neurons.total-Nb_neurons.Layer3, Nb_neurons.Layer1+Nb_neurons.Ex+1:Nb_neurons.total-Nb_neurons.Layer3) = RecIn;

columnformat = {'char', 'numeric'};
d = {'AtoL2',AtoL2; 'BtoL2', L1toIn; 'VIP', VIP; 'VIPtoL2Ex', VIPtoL2Ex;...
    'L1toVIP', L1toVIP; 'BtoVIP', BtoVIP; 'L2toVIP', L2toVIP; 'RecEx', RecEx; 'RecIn', RecIn; ...
    'CX2',CX2; 'EIbalance', EIbalance; 'f_D',0.6; 'f_I', 0.9;'L1',Nb_neurons.Layer1;'tau_pI',tau_pI };


for ind = Nb_neurons.Layer1+1:Nb_neurons.total-Nb_neurons.Layer3
    for ind2 = Nb_neurons.Layer1+1:Nb_neurons.total-Nb_neurons.Layer3
        W(ind,ind2) = W(ind,ind2)*(rand<CX2);
    end
end



% ind = find(W==1);
% for w = 1:length(ind)
%     p = rand;
%     if p >EIbalance
%         W(ind(w)) = RecIn;
%     else
%         W(ind(w)) = RecEx;
%     end
%
% end

W = W.';
%
% imagesc(W)
% colorbar

% Estim = 4.5;
% Istim = 4.5*1.7;
% Esyn= ones(Nb_neurons,1)*4.5;
% Isyn= ones(Nb_neurons,1)*4.5*1.7;
% listOFneurons = 20; %[1 20 21]

out.vector =[];
out.dis_rate.mean = [];
out.dis_rate.std = [];
if sim ==0
    L12neurons = Nb_neurons.Layer1:Nb_neurons.Layer1+Nb_neurons.Ex;% 1:Nb_neurons.total-Nb_neurons.Layer3; % Nb_neurons.Layer1+1:Nb_neurons.total-Nb_neurons.Layer3;
else
    L12neurons = 1:Nb_neurons.total-Nb_neurons.Layer3;
end


for neuronNB = L12neurons
    output.VS{neuronNB} = [];
    output.DRstd{neuronNB} = [];
    output.DRmean{neuronNB} = [];
    raster.stim{neuronNB}=[];
    raster.rep{neuronNB}=[];
    raster.spikes{neuronNB}=[];
end



for f = 1 : length(ICI_list)
    disp(f)
    tic
    Ge = [];
    Gi = [];
    if Posi ==1
        [Ge1, Gi1] = model_stim(4.5,8.5,0.9,0.6,10,f); %Sync+
    else
        [Ge1, Gi1] = model_stim(4.5,8.5,0.6,0.9,10,f); %Sync-
    end
    %     [Ge2, Gi2] = model_stim(4.5,8.5,0.5,1.0,10,f);
%     [Ge2, Gi2] = model_stim(4.5,8.5,0.9,0.6,10,f);
    Ge = zeros(Nb_neurons.total,length(Ge1));
    Gi = zeros(Nb_neurons.total,length(Gi1));
%     for n = 1:Nb_neurons.total
%         if n <= Nb_neurons.A
%             Ge(n,:) = Ge1;% = [Ge;Ge1];
%             Gi(n,:) = Gi1;%[Gi;Gi1];
% %         elseif n > Nb_neurons.A && n <= Nb_neurons.Layer1
% %             Ge = [Ge;Ge2];
% %             Gi = [Gi;Gi2];
% %         elseif n > Nb_neurons.Layer1 && n <= Nb_neurons.Layer2
% %             Ge = [Ge;zeros(1,length(Ge1))];
% %             Gi = [Gi;zeros(1,length(Gi1))];
% %             %             Gi = [Gi;Gi2];
% %             %             GiStim
% %         else
% %             Ge = [Ge;zeros(1,length(Ge1))]; %[ones(1,5000)*VIP*1e-9 zeros(1,5050) ones(1,4952)*VIP*1e-9]];
% %             Gi = [Gi;zeros(1,length(Gi1))];
%         end
%     end
    toc
    
    
    rate_total = [];
    spikes_pooled=[];
    indvec = [];
    for NB = L12neurons
        spikes_pooled_neuron{NB}= [];
    end
    for n = 1:nb_rep
        tic
        [V,spikes] = COBNmodel(Ge,Gi,Nb_neurons,W);
        
        %rate
        indvec = [indvec; L12neurons.'];
        for neuronNB = L12neurons% L12neurons%1:Nb_neurons.Layer1 %Nb_neurons.total-Nb_neurons.Layer3 % L2neurons %1:length(listOFneurons) 1:Nb_neurons.A%
            rate = zeros(1,1500);
            for st = spikes{neuronNB}
                rate(1,round(st/(stepp*10))) = rate(1,round(st/(stepp*10)))+1;
            end
            rate_total = [rate_total ; rate*1000];
            %pool spikes
            spikes{neuronNB}=spikes{neuronNB}-PREstimulus_duration;
            spikes_pooled=[spikes_pooled spikes{neuronNB}];
            spikes_pooled_neuron{neuronNB} = [spikes_pooled_neuron{neuronNB} spikes{neuronNB}] ;
            
            raster.stim{neuronNB}=[raster.stim{neuronNB} f*ones(size(spikes{neuronNB}))];
            raster.rep{neuronNB}=[raster.rep{neuronNB} n*ones(size(spikes{neuronNB}))];
            raster.spikes{neuronNB}=[raster.spikes{neuronNB} spikes{neuronNB}];
            
            
        end
        toc
    end

    
    %%% ind neurons
    
    if sim ==1
        for neuronNB = L12neurons %1:Nb_neurons.total-Nb_neurons.Layer3% L2neurons
            
            %% Analysis
            spikes_pooled_for_vector_strength=spikes_pooled_neuron{neuronNB}(find(spikes_pooled_neuron{neuronNB}>0.05 & spikes_pooled_neuron{neuronNB}<=(stimulus_duration+0.05)));
            freq2 = 1000./ICI_list(f);
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
            output.VS{neuronNB} = [output.VS{neuronNB} vector];
            
            %  average rate
            
            PRE = PREstimulus_duration*1000;
            POST = POSTstimulus_duration*1000;
            STIM = stimulus_duration*1000;
            total_time = PRE+POST+STIM;
            rate_av = mean(rate_total(find(indvec == neuronNB),:),1);
            spont_rate = mean2(rate_total(find(indvec == neuronNB),PRE-100:PRE));
            discharge_rate{neuronNB}.mean = mean2(rate_total(find(indvec == neuronNB),PRE+1:PRE+STIM+100))-spont_rate;
            discharge_rate{neuronNB}.std = std2(rate_total(find(indvec == neuronNB),PRE+1:PRE+STIM+100))/sqrt(nb_rep*(STIM+100));
            output.DRmean{neuronNB} = [output.DRmean{neuronNB} discharge_rate{neuronNB}.mean];
            output.DRstd{neuronNB} = [output.DRstd{neuronNB} discharge_rate{neuronNB}.std];
            output.ratebrut{neuronNB,f} = rate_av;
            xs = 1:total_time;
            h = 10; %kernal bandwidth. determines the shape of the function
            for i = 1:total_time
                ys(i)=gaussian_kern_reg(xs(i),xs,rate_av,h);
            end
            output.rate{neuronNB,f} = ys;
        end
        
        
        %%% population Analysis
    else
        spikes_pooled_for_vector_strength=spikes_pooled(find(spikes_pooled>0.05 & spikes_pooled<=(stimulus_duration+0.05)));
        freq2 = 1000./ICI_list(f);
        % %     for calculating vector strength, subtract the first 50 ms and
        % %     include 50 ms post stimulus
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
        
        out.vector = [out.vector vector]
        
        % %      average rate
        
        PRE = PREstimulus_duration*1000;
        POST = POSTstimulus_duration*1000;
        STIM = stimulus_duration*1000;
        total_time = PRE+POST+STIM;
        rate_av = mean(rate_total,1);
        spont_rate = mean2(rate_total(:,PRE-100:PRE));
        discharge_rate.mean = mean2(rate_total(:,PRE+1:PRE+STIM+100))-spont_rate;
        discharge_rate.std = std2(rate_total(:,PRE+1:PRE+STIM+100))/sqrt(nb_rep*(STIM+100));
        
        
        
        xs = 1:total_time;
        h = 10; %kernal bandwidth. determines the shape of the function
        for i = 1:total_time
            ys(i)=gaussian_kern_reg(xs(i),xs,rate_av,h);
        end
        out.dis_rate.mean = [out.dis_rate.mean discharge_rate.mean];
        out.dis_rate.std = [out.dis_rate.std discharge_rate.std];
        out.rate_brut{f} = rate_av;
        out.rate{f} = ys;
        
    end
end

output.raster = raster;
%
% figure
% plot(output.raster.spikes{446},20*(output.raster.stim{446}-1)+output.raster.rep{446},'k.','MarkerSize',9);
if sim ==1 %Ind neuron analysis
    output.Analysis = struct('Rho',{},'Pval',{},'Slope',{});
    for neuronNB = L12neurons
        counter1 = 0;
        [RHO,PVAL] = corr(ICI_list.',output.DRmean{neuronNB}.','Type','Spearman');
        for f = 1:length(ICI_list)
            if output.VS{neuronNB}(f) >0.1
                counter1 = counter1 + 1;
            else
                counter1 = 0;
            end
        end
        if counter1 < 3 && PVAL <0.05
            
            %         if RHO > 0
            %             disp('Positive')
            %         elseif RHO < 0
            %             disp('Negative')
            %         end
            mdl = fitlm(ICI_list,output.DRmean{neuronNB});
            slope = mdl.Coefficients.Estimate(2);
        else
            slope = 0;
        end
        output.Analysis(neuronNB).Rho = RHO;
        output.Analysis(neuronNB).Pval = PVAL;
        output.Analysis(neuronNB).Slope = slope;
        %     output{neuronNB}.Analysis = [RHO,PVAL, slope];
        %     output{neuronNB}.parameters = [RecEx RecIn CX2 slope];
    end
    NsyncIdx = find([output.Analysis.Slope] ~= 0);
    title = ['Network' int2str(iii) '.mat'];
    save(title,'output', 'ICI_list', 'L12neurons','W','NsyncIdx','d')
    % test = 1 ;
    
    total_DR = [];
    
    total_DRerror = [];
    for n = NsyncIdx
        total_DR = [total_DR ; output.DRmean{n}];
        total_DRerror = [total_DRerror;output.DRstd{n}];
    end
    
    
    Hz_list = [];
    for i = 1:length(ICI_list)
        Hz_list = [Hz_list round(1000/ICI_list(i))];
    end
    
    new_mean = mean(total_DR,1);
    new_error = mean(total_DRerror,1);
    
    errorbar(Hz_list, new_mean,new_error)
    test = 1;
    
else %%% population Analysis
    
    [RHO,PVAL] = corr(ICI_list.',out.dis_rate.mean.','Type','Spearman')
    if PVAL <0.05
        if RHO > 0
            disp('Positive')
        elseif RHO < 0
            disp('Negative')
        end
    else
        disp(PVAL)
    end
    %
    disp(d)
    % output.VS = out.vector
    % output.Spearman = [RHO,PVAL]
    % output.parameters = [RecEx RecIn CX2 slope];
    
    
    Hz_list = [];
    for i = 1:length(ICI_list)
        Hz_list = [Hz_list round(1000/ICI_list(i))];
    end
%     figure('units','normalized','outerposition',[0 0 1 1])
    cmap = colormap(jet(length(ICI_list)));
    subplot(2,2,1)
    for f = 1:length(ICI_list)
        plot(out.rate{f},'linewidth', 2.0, 'DisplayName', ...
            [num2str(ceil(1000/ICI_list(f))) 'Hz'] );
        lgd{f} = [num2str(ceil(1000/ICI_list(f))) 'Hz'];
        hold on
        % % % x = 510:ICI_list(f):1110;
        % % % x = [x; ones(1,length(x))*(f-1)];
        % % % scatter(x(1,:),x(2,:),'k','filled')
    end
 axis([0,1300,0 inf])
    set(gca, 'FontSize', 16)

    legend('show')
    
    subplot(2,2,2)
    uitable('Data',d,'ColumnFormat', columnformat)

% 
%     subplot(2,2,3)
%     plot(Hz_list,out.vector)
%     axis([0,50,0 1])
%     set(gca, 'FontSize', 16)
%     
%     
%     subplot(2,2,4)
%     %     errorbar(Hz_list,out.dis_rate.mean,out.dis_rate.std)
%     cmapp = [[0.9 0.1 0.1]; [0.1 0.2 0.9]; [0 0 0]   ];
%     shadedErrorBar(Hz_list,out.dis_rate.mean,out.dis_rate.std,{'--','Color',cmapp(2,:)})
%     axis([0,50,0 inf])
%     test = 1;
%     set(gca, 'FontSize', 16)
%     
end


%
% %         testing for monotonicity take between 8 and 50Hz


% [spikes1,V1,O1] = run_LIFmodel(Ge1,Gi1);
% [spikes2,V2,O2] = run_LIFmodel(Ge2,Gi2)
% test = 1


%%% /population Analysis

%% Create Conductance from stimuli
function [GeStim,GiStim] = model_stim(Estim, Istim,f_DE,f_DI,NbInputs,f)
global kernel
global tau_pE tau_pI IE_delay ICI_list stepp

kernel_time_constant=.005;  %time constant of 5 ms
% jitter_magnitude=1; % Temporal jitter

stimulus_duration=0.5;  %half second
PREstimulus_duration=0.5;  %half second
POSTstimulus_duration=0.5;  %half second (0.1 second is included in stimulus)
latency_time=0.01;
latency=length(0:stepp:latency_time); %10 ms latency (auditory nerve to auditory cortex)



ICI = ICI_list(f);
freq = 1000./ICI;
t=0:stepp:(kernel_time_constant*10);
kernel=t.*exp(-t/kernel_time_constant);
kernel=1e-9*kernel/max(kernel); %amplitude of 1 nS
input=zeros(size(0:stepp:(POSTstimulus_duration+stimulus_duration)));
stimulus_input_length=length(0:stepp:(stimulus_duration));
ipi=round(1/(freq*stepp)); %ipi=interpulse interval
% freq2=1/(stepp*ipi);

P_0 = 1;
E_strStim(1) = 1;
I_strStim(1) = 0.7;
E_strength_mean = [];
I_strength_mean = [];

%adaptation parameters


% f_DE = 0.9;
% if f_DE ~= 0
%     f_DI = 1.6-f_DE;
% else    % No adaptation
%     f_DE = 1;
%     f_DI = 1;
% end


Ge_total = [];
Gi_total = [];
E_input=input;
I_input=input;


for j=1:NbInputs %10 jitter excitatory and inhibitory inputs
    % for i=1:ipi:(stimulus_input_length-(length(kernel)/2))
    p = 0;
    P_relE(1) = P_0;
    P_relI(1) = P_0;
    for i=1:ipi:(stimulus_input_length)
        p = p+1; %Click number (starts with 1)
        jitter=round(randn(1)/(1000*stepp)); %1 ms jitter
        
        if (i+jitter)<1 || (i+jitter)>(length(input)-length(kernel))
            jitter=1;
        end
        t0 = i+jitter;
        if i == 1
            P_relE(1:t0) = P_0;
            P_relI(1:t0) = P_0;
        end
        for t = t0:t0+2*ipi-1
            P_relE(t+1) = P_0 + (f_DE*P_relE(t0)-P_0)*exp(-(t+1-t0)*stepp/tau_pE);
            P_relI(t+1) = P_0 + (f_DI*P_relI(t0)-P_0)*exp(-(t+1-t0)*stepp/tau_pI);
        end
        
        E_strStim(p) = P_relE(t0)*Estim;
        I_strStim(p) = P_relI(t0)*Istim;
        E_input((latency+t0):(latency+t0+length(kernel)-1))=E_input((latency+t0):(latency+t0+length(kernel)-1))+kernel*E_strStim(p);
        %             jitter=round(randn(1)/(1000*stepp)); %1 ms jitter
        %
        %             if (t0)<1 || (t0)>(length(input)-length(kernel))
        %                 jitter=0;
        %             end
        I_input((latency+t0):(latency+t0+length(kernel)-1))=I_input((latency+t0):(latency+t0+length(kernel)-1))+kernel*I_strStim(p);
        
    end
    E_strength_mean = [E_strength_mean ; E_strStim];
    I_strength_mean = [I_strength_mean ; I_strStim];
end

delay=round(abs(IE_delay)/(1000*stepp));  %delay in stepps
if IE_delay>=0
    GeStim=E_input;
    GiStim=[zeros(1,delay) I_input(1:(length(I_input)-delay))];
elseif IE_delay<0
    GiStim=I_input;
    GeStim=[zeros(1,delay) E_input(1:(length(E_input)-delay))];
end
GeStim=[zeros(size(0:stepp:PREstimulus_duration)) GeStim];
GiStim=[zeros(size(0:stepp:PREstimulus_duration)) GiStim];
% Ge = GeStim;
% Gi = GiStim;


%% neuron-network model
function [V,spikes] = COBNmodel(Ge,Gi,Nb_neurons,W) % NbSelf = 0 if no connections


C=0.25*1e-9; %0.25 nF, 10 ms time constant
Grest=25*1e-9; %25 nS
Erest=-0.065; %-65 mV
Ee=0; % 0 mV
Ei=-0.085;  %-85 mV
noise_magnitude=4e-8; %default noise level in conductance
P_0 = 1;

V = zeros(Nb_neurons.total,length(Ge));
O = zeros(Nb_neurons.total,length(Ge));
V(1:Nb_neurons.total,1)=Erest;
for n = 1:Nb_neurons.total
    spikes{n} = [];
    t0(n) = 1;
end

global tau_pE kernel IE_delay stepp

sigma = 0.01;

P_syn(1:Nb_neurons.total,1) = P_0;
tau_p(1:Nb_neurons.total,1) = tau_pE;
tau_p(Nb_neurons.Layer1+1:end,1) =0.5; %L2neurons
tau_p(Nb_neurons.Layer1+Nb_neurons.Layer2+1:end,1) = tau_pE; % VIP neurons
f_D(1:Nb_neurons.total,1) = 1;
% f_D(1:end-Nb_neurons.Layer3,1) = 0.7; %L2neurons changed from 0.7 10/05/2017
% f_D(Nb_neurons.Layer1+1:Nb_neurons.Layer1+Nb_neurons.Ex,1) =0.7;

f_D(1:Nb_neurons.Layer1) = 1; %inhibitory VIP neuronsS
% f_D(1:Nb_neurons.A,1) = 1;
% W = ones(Nb_neurons,Nb_neurons)*0.1; %Connectivity Weight matrix
% O = zeros(Nb_neurons,length(Ge)); %spike count vector
% t0(n) = 1;
% p = 0; %nb of spikes]
latency2_mean = length(0:stepp:0.003); %3ms, reccurent connection latency
latency2 = rand(Nb_neurons.total,Nb_neurons.total)*20+latency2_mean-20;
latency2 = round(latency2.');


for t = 1:length(Ge)
    for n = 1:Nb_neurons.total
        P_syn(n,t+1) = P_0 + (P_syn(n,t0(n))-P_0)*exp(-(t+1-t0(n))*stepp/tau_p(n));
        %slow recovery of VIP neurons
        if n > Nb_neurons.Layer1 + Nb_neurons.Layer2
            if t<13000
                Gi(n,t+1+max(latency2(:,n))) = Gi(n,1) + (Gi(n,t)-Gi(n,1))*exp(-(1+max(latency2(:,n))*stepp/tau_p(n)));
            end
        end
        %Add noise to conductance
        if n ~= Nb_neurons.total
            Ge(n,t) = Ge(n,t)+noise_magnitude*randn;
            Gi(n,t) = Gi(n,t)+noise_magnitude*randn;
        end
        if Ge(n,t) <0
            Ge(n,t) = 0;
        end
        if  Gi(n,t) <0
            Gi(n,t) = 0;
        end
        
        %voltage update
        V(n,t+1)=(-stepp*( Ge(n,t)*(V(n,t)-Ee) + Gi(n,t)*(V(n,t)-Ei) + Grest*(V(n,t)-Erest))/C)+V(n,t) + sigma*randn*sqrt(stepp);
        
        if V(n,t+1)>(Erest+0.020) %20 mV above Erest %artificial threshold
            if O(n,t) > 0 % refractory period : if there was a spike before
                V(n,t+1)=Erest;
                P_syn(n,t+1) = P_syn(n,t);
            else
                V(n,t+1)=0.050; %spike to 50 mV
                spikes{n}=[spikes{n} stepp*(t+1)];
                O(n,t+1) = 1;
                %                 i=i+1;
                t0(n)=t+1;
                %             if NbSelf ~= 0
                P_syn(n,t+1) = P_syn(n,t)*f_D(n);
                %                 p = p+1;
                %                 J(p) = P_syn(n,t+1);
                %add kernel to all Ge
                for nb = 1:Nb_neurons.total
                    if 1<t0(n) && t0(n)<13000
                        if W(n,nb)>0
                            Ge(nb,(latency2(n,nb)+t0(n)):(latency2(n,nb)+t0(n)+length(kernel)-1)) ...
                                = Ge(nb,(latency2(n,nb)+t0(n)):(latency2(n,nb)+t0(n)+length(kernel)-1))+ kernel*W(n,nb)*P_syn(n,t+1);% +noise_magnitude*randn(1,length(kernel));
                            %                             Gi(nb,(latency2+t0(n)):(latency2+t0(n)+length(kernel)-1))...
                            %                                 = Gi(nb,(latency2+t0(n)):(latency2+t0(n)+length(kernel)-1))+ kernel*W(n,nb)*P_syn(n,t+1);
%                         elseif W(n,nb) == 1.1
%                             Ge(nb,(latency2(n,nb)+500+t0(n)):(latency2(n,nb)+t0(n)+500+length(kernel)-1)) ...
%                                 = Ge(nb,(latency2(n,nb)+500+t0(n)):(latency2(n,nb)+500+t0(n)+length(kernel)-1))+ kernel*W(n,nb)*P_syn(n,t+1);% +noise_magnitude*randn(1,length(kernel));
                        elseif W(n,nb)<0
                            Gi(nb,(latency2(n,nb)+t0(n)):(latency2(n,nb)+t0(n)+length(kernel)-1)) = Gi(nb,(latency2(n,nb)+t0(n)):(latency2(n,nb)+t0(n)+length(kernel)-1))+ kernel*-W(n,nb)*P_syn(n,t+1);
                        end
                    end
                end
                %             Ge(find(Ge<0))=0;
                %             Gi(find(Ge<0))=0;
                %             end
            end
        end
        
        
    end
end
% figure
% subplot(2,2,1)
% plot(V(21,:))
% subplot(2,2,2)
% plot(V(20,:))
% subplot(2,2,3)
% plot(V(1,:))
% subplot(2,2,4)
% plot(Ge(21,:)-Gi(21,:))
% figure
% plot(Gi(121,:))
% test = 1;


% function out = SpkAnalysis(V,spikes)





%
% %% Input neurons (Sync neurons?)
% function [spikes,V,O]=run_LIFmodel(Ge,Gi) %,f)
% spikes=[]; V=[]; t=1; i=1;
% stepp=.0001; %.1 ms duration  (temporal increment for running simulation)
% C=0.25*1e-9; %0.25 nF, 10 ms time constant
% Grest=25*1e-9; %25 nS
% Erest=-0.065; %-65 mV
% Ee=0; % 0 mV
% Ei=-0.085;  %-85 mV
%
%
% V(1)=Erest;
% O = zeros(1,length(Ge));
% noise_magnitude=4e-8; %default noise level in conductance
%
% %avoid negative conductances
% Ge=Ge+noise_magnitude*randn(1,length(Ge));
% Gi=Gi+noise_magnitude*randn(1,length(Gi));
% Ge(find(Ge<0))=0;
% Gi(find(Gi<0))=0;
% sigma = 0.01    ;
%
% %synaptic depression
% % global ICI_list
% global tau_pE tau_pI
%
% % ICI = ICI_list(f);
% ICI = 20;
% freq = 1000./ICI;
% % ipi=round(1/(freq*stepp)); %ipi=interpulse interval
%
% while(t<length(Ge))
%
%     V(t+1)=(-stepp*( Ge(t)*(V(t)-Ee) + Gi(t)*(V(t)-Ei) + Grest*(V(t)-Erest))/C)+V(t) + sigma*randn*sqrt(stepp);
%     if V(t+1)>(Erest+0.020) %20 mV above Erest %artificial threshold
%         V(t+1)=0.050; %spike to 50 mV
%         spikes(i)=stepp*(t+1);
%         O(t+1) = 1;
%         i=i+1;
%         t=t+1;
%         V(t+1)=Erest;
%     end
%     t=t+1;
% end
%
% %% kernel function
% function kernel = knl(t,t0)
% max_kernel = 0.0018;
% kernel=t.*exp(-t/kernel_time_constant);
% kernel=1e-9*kernel/max_kernel;
%
%







