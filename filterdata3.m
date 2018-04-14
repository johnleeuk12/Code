function filterdata3()

clear all

% will be needed when taking all data
%load('infoset.mat');
%animals = unique(cellfun(@char,{output.animal},'unif',0));
% % % 

%Extracting data for all neurons of a certain category, regardless of stim
%set type. 


%% file info

Sync = 1; %or 0 for Nsync
Positive = 1; % or 1 for Negative
StimType = 12; % or 20

ICI_list1 = [2 2.5 3 5 7.5 10 12.5 15 20 25 30 35 40 45 50 55 60 65 70 75]; %ms, ICI
ICI_list2 = [250 125 83.3333 62.5 50 41.6667 35.7143 31.25 27.7778 25 22.7273 20.8333];

% if StimType == 12
%     ICI_list = ICI_list2;
% else % UnitInfo.Info(2,n) == 20
%     ICI_list = ICI_list1;
% end


% createFileInfo('m2p') %animal namecode
animal_list = {'m36n','m2p','m41o','m32q'};

load('Analysis.mat')
% Analysis.SyncPosi = {};
% Analysis.SyncNega = {};
% Analysis.NsyncPosi = {};
% Analysis.NsyncNega = {};
ICI_list = [25 50];

for f = 1:2
    output.rates_stim{f} = [];
    output.rates_pre{f} = [];
    output.rates_post{f} = [];
    output.isi_total{f} = [];
    output.mean_neuron_rates{f} = [];
end
output.names = [];
nn = 1;

% freq_list=round(1000./ICI_list);


for ind = 1:4
    animal = animal_list{ind};
    load([animal '_List3']);
    
    %% spikes, VS and raster
%     directory = ['U:\Neural and Behavioural Data' '\marmoset\' animal];
    directory = ['C:\Users\John\Documents\Marmoset\' animal];
    indxList = find( [UnitInfo.Info.Sync] == Sync & [UnitInfo.Info.Positive] == Positive ...
         & [UnitInfo.Info.Pre_stim_Duration] == 500 & [UnitInfo.Info.Post_stim_Duration] == 500); %& [UnitInfo.Info.Significant_rate] ==1 );
    
    for n = indxList
        disp(UnitInfo.List{1,n})
        output.names(nn,1) = ind;
        spiketable = load([directory filesep UnitInfo.List{1,n}]);
        output.names(nn,2) = str2num(UnitInfo.List{1,n}(length(animal)+1:end-4));
        output.names(nn,3) = UnitInfo.Info(n).Channel_Nb;
        
        PREstim_duration = UnitInfo.Info(n).Pre_stim_Duration/1000; %seconds
        stim_duration = 0.5;
        POSTstim_duration = UnitInfo.Info(n).Post_stim_Duration/1000;
        TrialLength = UnitInfo.Info(n).Pre_stim_Duration + UnitInfo.Info(n).Post_stim_Duration +500;
        
        raster.stim=[];  raster.rep=[];  raster.spikes=[];
        emptycount = 0;
        
        
        if UnitInfo.Info(n).Stimuli_Nb == 12
            ICI_ind = [find(ICI_list2 == 25) find(ICI_list2 == 50)];
        else
            ICI_ind = [find(ICI_list1 == 25) find(ICI_list1 == 50)];
        end
        
        for f = ICI_ind
            ff = find(ICI_ind == f);
            indstim = find(spiketable(:,1)==f);
            nbrep = spiketable(find(spiketable(:,1)==f),2);
            nreps = max(nbrep.');
            spikes_pooled = [];
            rate_total = [];
            isi_total = [];
            if ~isempty(nbrep)
                for r = unique(nbrep.')
                    spikes1 = []; %channel 1
                    spikes2 = []; %channel 2
                    for i = indstim.'
                        if spiketable(i,2) == r && spiketable(i,3)== 1 && spiketable(i,4) > 0 % channel 1
                            spikes1 = [spikes1 (spiketable(i,4)*1e-6 - PREstim_duration)]; % spikes rounded up to 0.1ms cad 100microseconds. converted to seconds
                        end
                        if spiketable(i,2) == r && spiketable(i,3)== 2 && spiketable(i,4) > 0 % channel 1
                            spikes2 = [spikes2 (spiketable(i,4)*1e-6 - PREstim_duration)];
                        end
                    end
                    if length(spikes1) > length(spikes2)
                        spikes =  spikes1;
                    else
                        spikes = spikes2;
                    end
                    if isempty(spikes(find(spikes>PREstim_duration & spikes< PREstim_duration + 0.5)))
                        emptycount = emptycount+1;
                    end
                    
                    spikes_pooled = [spikes_pooled spikes];
                    
                    raster.stim=[raster.stim f*ones(size(spikes))];
                    raster.rep=[raster.rep r*ones(size(spikes))];
                    raster.spikes=[raster.spikes spikes];
                    
                    % rate
                    rate = zeros(1,TrialLength);
                    spikes4rate = spikes(find(spikes<PREstim_duration+0.5)) + PREstim_duration;
                    
                    
                    for st = spikes4rate
                        if ceil(st*1e3) <= length(rate)
                            rate(1,ceil(st*1e3)) = rate(1,ceil(st*1e3))+1;
                        end
                        
                    end
                    rate_total = [rate_total ; rate*1000];
%                     spikes3 = spikes(find(spikes<PREstim_duration+0.5));
%                     spikes4 = [0 spikes3(1:end-1)];
%                     isi = spikes3-spikes4;
%                     isi = isi(2:end);
%                     output.isi_total{f} = [output.isi_total{f} isi];
                end
                spikes_pooled_for_vector_strength=spikes_pooled(find(spikes_pooled>0.05 & spikes_pooled<=(stim_duration+0.05)));
                freq2 = 1000./ICI_list(ff);
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
                    rayleigh(ff)=0;
                end
                
                if rayleigh<13.8
                    vector=0;
                end
                output.VS(nn,ff) = vector;
                
                %  average rate
                
                PRE = PREstim_duration*1000;
                POST = POSTstim_duration*1000;
                STIM = stim_duration*1000;
                total_time = PRE+POST+STIM;
                %         rate_av = mean(rate_total,1);
                %         spont_rate = mean2(rate_total(find(indvec == neuronNB),PRE-100:PRE));
                %         discharge_rate{neuronNB}.mean = mean2(rate_total(find(indvec == neuronNB),PRE+1:PRE+STIM+100))-spont_rate;
                %         discharge_rate{neuronNB}.std = std2(rate_total(find(indvec == neuronNB),PRE+1:PRE+STIM+100))/sqrt(nb_rep*(STIM+100));
                %         output.DRmean{neuronNB} = [output.DRmean{neuronNB} discharge_rate{neuronNB}.mean];
                %         output.DRstd{neuronNB} = [output.DRstd{neuronNB} discharge_rate{neuronNB}.std];
                
            end
            %             spikes_pooled_for_VS=spikes_pooled(find(spikes_pooled>0.05 & spikes_pooled<=(stim_duration+0.05)));
            
            %Interspike interval
            
            %             var_isi = [var_isi std(isi)];
            output.rates_pre{ff} = [output.rates_pre{ff}; ...
                rate_total(:,1:PRE)];
            output.rates_stim{ff} = [output.rates_stim{ff}; ...
                rate_total(:,PRE+1:PRE+STIM+50)];
            output.mean_neuron_rates{ff} = [output.mean_neuron_rates{ff}; ...
                mean(rate_total(:,PRE+1:PRE+STIM+50),1)];
            output.rates_post{ff} = [output.rates_post{ff}; ...
                rate_total(:,PRE+STIM+51:end)];
        end
        nn = nn+1;
    end
end

disp(nn)
% Rasterplot 
% x = freq_list;
% xlabel('time (s)')
% % ylabel('IPI (ms)')
% ylabel('Repetition rate (Hz)')
% % area([0 500 stimulus_duration 0],[length(x)*nreps+1 length(x)*nreps+1 0 0],'LineStyle','none','FaceColor',[.85 .85 1]);
% hold on
% plot(raster.spikes,nreps*(raster.stim-1)+raster.rep,'k.','MarkerSize',9);
% axis([0 stimulus_duration+POSTstimulus_duration 0 length(x)*nreps+1])

% Rasterplot end

%     test =1;
%     all_rate_total = [];
output.Fanofactor = [];
% vector = zeros(1,length(ICI_list));

output.meanVS = mean(output.VS,1);
tsVS = (tinv([0.025  0.975],nn));
output.errorVS = std(output.VS,1)/sqrt(nn)*tsVS(2);



output.meanDR = [];
output.errorDR = [];
output.totalDR = [];
output.totalDR_yaxis = [];
output.var_isi = [];
% output.meanspt = [];
for f = 1:length(ICI_list)
    output.mean_rate_stim{f} = mean(output.rates_stim{f},1);
    output.std_rate_stim{f} = std(output.rates_stim{f},0,1);
    output.mean_rate_pre{f} = mean(output.rates_pre{f},1);
    output.std_rate_pre{f} = std(output.rates_pre{f},0,1);
    output.mean_rate_post{f} = mean(output.rates_post{f},1);
    SpikeCount = sum(output.rates_stim{f},2)/1000.;
    output.Fanofactor = [output.Fanofactor std(SpikeCount)^2/mean(SpikeCount)];
    output.var_isi = [output.var_isi std(output.isi_total{f})];
    
    if Sync == 1
        %spikes per click
        output.spikes_per_click{f} = {};
        p = floor(500/ICI_list(f));
        for q = 1:p
            clicktime = round((q-1)*ICI_list(f))+50; %plus input latency + kernel peak
            output.spikes_per_click{f}.mean(q) = mean2(output.rates_stim{f}(:,clicktime - 10 : clicktime + 10));
            output.spikes_per_click{f}.sem(q) = std2(output.rates_stim{f}(:,clicktime - 10 : clicktime + 10))/sqrt(128*20);  %460 for SYnc + 355 for SYnc-
            output.spikes_per_click{f}.xaxis(q) = clicktime+500;
            output.spikes_per_click{f}.brut(:,q) = mean(output.mean_neuron_rates{f}(:,clicktime - 20 : clicktime + 20),2); %-mean(output.rates_pre{f},1);
            output.spikes_per_click{f}.brut(:,q) = (output.spikes_per_click{f}.brut(:,q)-mean2(output.rates_pre{f}))/(mean2(output.rates_stim{f})-mean2(output.rates_pre{f}));
        end
        
    end
    
    
    
    %Gaussian smoothing
    xs = 1:1500;
    h = 10;
    for i = 1:1500
        ys(i) = gaussian_kern_reg(xs(i),xs,[output.mean_rate_pre{f} output.mean_rate_stim{f} output.mean_rate_post{f}],h);
    end
    smooth_rate{f} = ys;
    output.meanDR = [ output.meanDR mean2(output.rates_stim{f})-mean2(output.rates_pre{f})];
%     output.meanspt = [output.meanspt mean2(output.rates_pre{f})];
    SEM = std2(output.rates_stim{f})/sqrt(size(output.rates_stim{1},1)*size(output.rates_stim{1},2));
    ts = tinv([0.025  0.975],size(output.rates_stim{1},1)*size(output.rates_stim{1},2)-1);      % T-Score
    %     CI = mean2(output.rates_stim{f}) + ts(2)*SEM;                      % Confidence Intervals
    output.errorDR = [output.errorDR  ts(2)*SEM];
%     for g = 1:46
%         output.totalDR = [output.totalDR; mean2(output.rates_stim{f}((g-1)*10+1:g*10,:))];
%     end
%     output.totalDR_yaxis = [output.totalDR_yaxis; ones(46,1)*f];
end


output.onset = mean(output.mean_rate_stim{1,1}(1,20:30),2);
Hz_list = [];
for i = 1:length(ICI_list)
    Hz_list = [Hz_list round(1000/ICI_list(i))];
end

[RHO,PVAL] = corr(ICI_list(2:end).',output.meanDR(2:end).','Type','Spearman')

% save('testoutput.mat', 'output')

cmapp = [[0.26 0.5 0.9]; [0.9 0.3 0.26]; [0 0 0];  [0.1 0.7 0.1]; [0.9 0.6 0.1] ]; %Green Yellow Black Blue Red : Nsync+ Nsync- Model Sync+ Sync-

color = cmapp(3*Sync + 1.5 + Positive/2,:);

cmap = colormap(jet(length(ICI_list)));


% ff = 12;
% % norm48Hz = output.spikes_per_click{ff}.mean/output.spikes_per_click{ff}.mean(1);
% % normSEM = output.spikes_per_click{ff}.sem/sqrt(output.spikes_per_click{ff}.mean(1));
% % plot(output.spikes_per_click{ff}.xaxis,norm48Hz,'color',cmapp(4,:),'LineWidth',1.7)
% errorbar(output.spikes_per_click{ff}.xaxis,norm48Hz,normSEM,'color',cmapp(5,:),'LineWidth',1.7)
% axis([500,1100,0.2,1.2])

set(gca, 'FontSize', 16)
hold on

figure
for f = 1: length(ICI_list)    
errorbar(output.spikes_per_click{f}.xaxis,output.spikes_per_click{f}.mean,output.spikes_per_click{f}.sem, 'color',cmap(f,:),'LineWidth',1.7,'DisplayName', ...
        [num2str(ceil(1000/ICI_list(f))) 'Hz'])
hold on
end
legend('show')
set(gca, 'FontSize', 16)
test = 1;
xlabel('Time(ms)')
ylabel('Average number of spikes per click (Spikes/sec)')


titre = ['RHO= ' num2str(RHO) ', p value = ' num2str(PVAL)];
figure 
subplot(2,3,2)
shadedErrorBar(Hz_list,output.meanVS,output.errorVS,{'Color',color})
xlabel('Repetition rate (Hz)')
ylabel('Firing rate (Spikes/sec)')
title('Vector Strength')
axis([0,50,-0,1])
set(gca, 'FontSize', 16)
subplot(2,3,3)

shadedErrorBar(Hz_list,output.meanDR,output.errorDR,{'Color',color})
ylabel('Firing rate (Spikes/sec)')
xlabel('Repetition rate (Hz)')
title('Discharge Rate')
axis([0,50,0,40])
set(gca, 'FontSize', 16)

cmap = colormap(jet(length(ICI_list)));
subplot(2,3,[1 4])

for f = 1: length(ICI_list)
    plot(smooth_rate{f},'color',cmap(f,:),'LineWidth',1.7,'DisplayName', ...
        [num2str(ceil(1000/ICI_list(f))) 'Hz'] );
    %     lgd{f} = [num2str(ceil(1000/ICI_list(f))) 'Hz'];
    hold on
end
legend('show')
axis([300,1200,0,80])
title(titre)
set(gca, 'FontSize', 16)
% 
% 
% 
subplot(2,3,5)
plot(Hz_list,output.Fanofactor,'Color',color)
title('FanoFactor')
xlabel('Repetition rate (Hz)')
subplot(2,3,6)
plot(Hz_list,output.var_isi,'Color',color)
title('Variance of ISI')
xlabel('Repetition rate (Hz)')


