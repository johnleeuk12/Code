%% test 



% posi = [];
% nega = [];
% null = [];
% animal_list = {'m36n','m2p','m41o','m32q'};
% for ind =2:4% 2:4
%     animal = animal_list{ind};
%     load([animal '_List3']);
%     
%     directory = ['U:\Neural and Behavioural Data' '\marmoset\' animal];
%     %     directory = ['C:\Users\John\Documents\Marmoset\' animal];
%     indxList = find([UnitInfo.Info.Stimuli_Nb] ==12 &... %[UnitInfo.Info.Positive] == Positive &...
%         [UnitInfo.Info.Sync] == 1 &  ...
%         [UnitInfo.Info.Pre_stim_Duration] == 500 &...
%         [UnitInfo.Info.Post_stim_Duration] == 500 & [UnitInfo.Info.Significant_rate] ==1 );
%     % & [UnitInfo.Info.Stimuli_Nb] ==StimType &[UnitInfo.Info.Positive] == Positive &
%     
%     for n =  indxList
%         if UnitInfo.Info(n).pval<=0.05
%             posi = [posi UnitInfo.Info(n).Rho];
% %         elseif UnitInfo.Info(n).Positive == 1
% %             nega = [nega UnitInfo.Info(n).Rho];
%         else
%             null = [null UnitInfo.Info(n).Rho];
%         end
%     end
% end
% 
%     




%% Mutual information population
% load('SyncN_new.mat')
% SyncNH = round([sum(output.rates_stim{3},2) sum(output.rates_stim{4},2) sum(output.rates_stim{5},2)]);
% SyncNH = reshape(SyncNH,[],1)/1000;
% SyncNL = round([sum(output.rates_stim{8},2) sum(output.rates_stim{9},2) sum(output.rates_stim{10},2)]);
% SyncNL = reshape(SyncNL,[],1)/1000;
% SyncN = [SyncNH SyncNL];
% SyncN = reshape(SyncN,[],1);
% 
% load('SyncP_new.mat')
% SyncPH = round([sum(output.rates_stim{3},2) sum(output.rates_stim{4},2) sum(output.rates_stim{5},2)]);
% SyncPH = reshape(SyncPH,[],1)/1000;
% SyncPL = round([sum(output.rates_stim{8},2) sum(output.rates_stim{9},2) sum(output.rates_stim{10},2)]);
% SyncPL = reshape(SyncPL,[],1)/1000;
% SyncP = [SyncPH SyncPL];
% SyncP = reshape(SyncP,[],1);
% rangeP = max(SyncP);
% rangeN = max(SyncN);
% 
% pdN = zeros(1,rangeN+1); %0 to rangeN)
% pdNH = zeros(1,rangeN+1); %0 to rangeN)
% pdNL = zeros(1,rangeN+1); %0 to rangeN)
% for i = 1:rangeN+1
%     pdN(i) = sum(SyncN == i);
%     pdNH(i) = sum(SyncNH ==i);
%     pdNL(i) = sum(SyncNL ==i);
% end
% 
% pdP = zeros(1,rangeP+1); %0 to rangeN)
% pdPH = zeros(1,rangeP+1); 
% pdPL = zeros(1,rangeP+1); 
% for i = 1:rangeP+1
%     pdP(i) = sum(SyncP == i);
%     pdPH(i) = sum(SyncPH == i);
%     pdPL(i) = sum(SyncPL == i);
% end
% pdN = pdN/sum(pdN);
% pdNH = pdNH/sum(pdNH);
% pdNL = pdNL/sum(pdNL);
% 
% pdP = pdP/sum(pdP);
% pdPH = pdPH/sum(pdPH);
% pdPL = pdPL/sum(pdPL);
% 
% pd = pdN.'*pdP;
% pdH = pdNH.'*pdPH;
% pdL = pdNL.'*pdPL;
% sumH = log2((pdH*2)./pd).*pdH;
% sumH(isnan(sumH)) = 0;
% sumL = log2((pdL*2)./pd).*pdL;
% sumL(isnan(sumL)) = 0;
% 
% sumPH = log2(pdPH*2./pdP).*pdPH;
% sumPH(isnan(sumPH)) = 0;
% 
% sumPL = log2(pdPL*2./pdP).*pdPL;
% sumPL(isnan(sumPL)) = 0;
% 
% InfoP = sum(sum(sumPH)) + sum(sum(sumPL));
% Info = sum(sum(sumH))+sum(sum(sumL));
% 
% % % clear all

%% Mutual information neuron by neuron
% 

% % h =  hist(output.isi_total{6});
% test = output.isi_total{2}*10^3;
% for i = length(test):-1:1
%     if test(i) < 10
%         test(i) = [];
%     end
% end
% edges = 0:5:500;
% h = histc(test,edges);
% % histogram(test,edges);
% % figure
% plot(h)
% set(gca,'XScale','log')
% hold on
% 
% [f, xi] = ksdensity(output.VS(1,:));
% plot(xi,f)
% edges = 0:0.1:1;
% histogram(output.VS(:,2),edges)


%% Calculating MI 07/11/2017
clear all
load('SyncP_new.mat')
NP = size(output.rates_stim{1, 1},1);
DR_trialP = [];
for f = 1:12
    for n = 1:NP
        DR_trialP(n,f) = mean(output.rates_stim{f}(n,:));
    end
end

% 08/11/2017 attempt to linearly extrapolate a pdf.
edgesTotalP = 0:1: ceil(max(max(DR_trialP)));
[NtotalP,edgesTotalP] = histcounts(DR_trialP,edgesTotalP);
for i = 1:length(NtotalP)
    if NtotalP(i) == 0
        NtotalP(i) = NaN;
    end
end
xaxisP = edgesTotalP(1:end-1);
% plot(xaxis,Ntotal)
% linear interpolation into infinite dataset
[xDataP, yDataP] = prepareCurveData( xaxisP, NtotalP );
vq1 = interp1(xDataP, yDataP, xaxisP);
pdfPall = vq1/sum(vq1);
% vq2 = interp1(xData, yData, xaxis,'spline');
% plot(xaxis,vq1)
% hold on 
% plot(xaxis,vq2)


InfoP = 0;
for f= 1:12
    [NcountP,edgesTotalP] = histcounts(DR_trialP(:,f),edgesTotalP);
    [xDataP, yDataP] = prepareCurveData( xaxisP, NcountP );
    vq2N = interp1(xDataP, yDataP, xaxisP);
    pdfP(f,:) = vq2N/sum(vq2N);
    for r = 1:length(pdfPall)
        if pdfP(f,r) ~=0 %  pdfPall(r) > 0.001 &&
        InfoP = InfoP + pdfP(f,r)*log2(pdfP(f,r)/pdfPall(r));
        end
    end
end
InfoP = InfoP/12;


load('SyncN_new.mat')

%for NM neurons, randomly select 22neurons



N = size(output.rates_stim{1, 1},1);
DR_trialN = [];
Nall = 300;
InfosN = [];
for iii = 1:100
    pdfN = [];
    Indice = randsample([1:N],Nall);
    
    for f = 1:12
        for n = 1:Nall
            DR_trialN(n,f) = mean(output.rates_stim{f}(Indice(n),:));
        end
    end
    
    % 08/11/2017 attempt to linearly extrapolate a pdf.
    edgesTotalN = 0:1: ceil(max(max(DR_trialN)));
    [NtotalN,edgesTotalN] = histcounts(DR_trialN,edgesTotalN);
    for i = 1:length(NtotalN)
        if NtotalN(i) == 0
            NtotalN(i) = NaN;
        end
    end
    xaxisN = edgesTotalN(1:end-1);
    % plot(xaxis,Ntotal)
    % linear interpolation into infinite dataset
    [xDataN, yDataN] = prepareCurveData( xaxisN, NtotalN );
    vq1 = interp1(xDataN, yDataN, xaxisN);
    pdfNall = vq1/sum(vq1);
    % vq2 = interp1(xData, yData, xaxis,'spline');
    % plot(xaxis,vq1)
    % hold on
    % plot(xaxis,vq2)
    
    
    InfoN = 0;
    for f= 1:12
        [NcountN,edgesTotalN] = histcounts(DR_trialN(:,f),edgesTotalN);
        [xDataN, yDataN] = prepareCurveData( xaxisN, NcountN );
        vq2N = interp1(xDataN, yDataN, xaxisN);
        pdfN(f,:) = vq2N/sum(vq2N);
        for r = 1:length(pdfNall)
            if  pdfN(f,r) ~=0  % pdfNall(r) > 0.001 &&
                InfoN = InfoN + pdfN(f,r)*log2(pdfN(f,r)/pdfNall(r));
            end
        end
    end
    InfoN = InfoN/12;
    
%     InfoPandN = 0;
%     for f = 1:12
%         %     [NcountN,edgesTotalN] = histcounts(DR_trialN(:,f),edgesTotalN);
%         %     [xDataN, yDataN] = prepareCurveData( xaxisN, NcountN );
%         %     vq2N = interp1(xDataN, yDataN, xaxisN);
%         %     pdfN(f,:) = vq2N/sum(vq2N);
%         %     [NcountP,edgesTotalP] = histcounts(DR_trialP(:,f),edgesTotalP);
%         %     [xDataP, yDataP] = prepareCurveData( xaxisP, NcountP );
%         %     vq2P = interp1(xDataP, yDataP, xaxisP);
%         %     pdfP(f,:) = vq2P/sum(vq2P);
%         for r1 = 1:length(pdfPall)
%             for r2 = 1:length(pdfNall)
%                 if   pdfP(f,r1)*pdfN(f,r2) ~=0 %pdfPall(r1)*pdfNall(r2) > 0.001 &&
%                     InfoPandN = InfoPandN + ...
%                         pdfP(f,r1)*pdfN(f,r2)*log2( pdfP(f,r1)*pdfN(f,r2)...
%                         /(pdfPall(r1)*pdfNall(r2)));
%                 end
%             end
%         end
%     end
%     
%     InfoPandN = InfoPandN/12;
    InfosN = [InfosN InfoN];
end
% 08/11/2017 attempt to calculate pdf by histograms
% % test = mean(DR_trial,1);
% % edges = 0:1: floor(max(max(DR_trial)));
% Nbins = 219;
% % Ntotal = histogram(DR_trial,edges);
% [Ntotal,edgesTotal] = histcounts(DR_trial,Nbins);
% 
% pdfPall =  Ntotal/sum(Ntotal);
% InfoP = 0;
% for f= 1:12
%     [Ncount,edges] = histcounts(DR_trial(:,f),Nbins);
%     pdfP(f,:) = Ncount/sum(Ncount);
%     for r = 1:length(pdfPall)
%         if pdfPall(r) > 0.001 && pdfP(f,r) ~=0 
%         InfoP = InfoP + pdfP(f,r)*log2(pdfP(f,r)/pdfPall(r));
%         end
%     end
% end
% % Ntotal = histogram(DR_trial,edges);
% 




%% Caculate MI with ISI and with VS
clear all

load('SyncN_new.mat')

N = size(output.VS,1);
Nall = 25;
InfosVS = [];
for iii = 1:100
Indice = randsample([1:N],Nall);
SampleVS = [];

    pdfVS = [];
    
    for f = 1:12
        for n = 1:Nall
            SampleVS(n,f) = output.VS(Indice(n),f);
        end
    end
    
    
    edgesTotalN = [0 0.2:0.1:1];
    xaxis =  edgesTotalN(1:end-1); %+0.05;
    [NtotalN,edgesTotalN] = histcounts(SampleVS,edgesTotalN);
    pdfVSall = NtotalN/sum(NtotalN);
%     plot(xaxis,pdfVSall);
%     pause
    InfoVS = 0;
    for f= 1:12
        [NcountN,edgesTotalN] = histcounts(SampleVS(:,f),edgesTotalN);
        pdfVS(f,:) = NcountN/sum(NcountN);
        for r = 1:length(pdfVSall)
            if  pdfVS(f,r) ~=0  % pdfNall(r) > 0.001 &&
                InfoVS = InfoVS + pdfVS(f,r)*log2(pdfVS(f,r)/pdfVSall(r));
            end
        end
    end
    
    InfoVS = InfoVS/12;
    InfosVS = [InfosVS InfoVS];
end


% ISI
% 
% clear all
% load('SyncP_new.mat')
% InfosISI = [];
% % for iii = 1:100
%     ISItotal = [];
%     for f = 1:12
%         ISItotal = [ISItotal output.isi_total{f}];
%     end
%     edgesISI = 0:1e-3:max(ISItotal);
%     [ISItotalhist edgesISI] = histcounts(ISItotal,edgesISI);
%     pdfISItotal = ISItotalhist/sum(ISItotalhist);
%     
%     InfoISI = 0;
%     pdfISI = [];
%     for f = 1:12
%         [ISIcount, edgesISI] = histcounts(output.isi_total{f},edgesISI);
%         pdfISI(f,:) = ISIcount/sum(ISIcount);
%         for r = 1:length(pdfISItotal)
%             if  pdfISI(f,r) ~=0  % pdfNall(r) > 0.001 &&
%                 InfoISI = InfoISI + pdfISI(f,r)*log2(pdfISI(f,r)/pdfISItotal(r));
%             end
%         end
%     end
%     
%     
%     InfoISI = InfoISI/12;
%     InfosISI = [InfosISI InfoISI]
% % end
% xaxis = edgesISI(1:end-1);
% % plot(xaxis, ISItotalhist);


%ISI ind neurons : 15/11/2017

clear all
load('SyncN_new.mat')
N = size(output.isi_rep,2);
Nall = 25;

InfosISI = [];
for iii =1:100
    Indice = randsample([1:N],Nall);
    SampleISI = [];
    ISIperStim = [];
    ISItotal = [];
    for f = 1:12
        ISIperStim{f} = [];
        for n = 1:Nall
            SampleISI{n,f} = [];
            for reps = 1:size(output.isi_rep{Indice(n)},2)
                SampleISI{n,f} = [SampleISI{n,f} output.isi_rep{Indice(n)}{f,reps}];
            end
            ISIperStim{f} = [ISIperStim{f} SampleISI{n,f}];
        end
        ISItotal = [ISItotal ISIperStim{f}];
    end
    edgesISI = 0:1e-3:max(ISItotal);
    [ISItotalhist edgesISI] = histcounts(ISItotal,edgesISI);
    
    pdfISItotal = ISItotalhist/sum(ISItotalhist);
    
    InfoISI = 0;
    pdfISI = [];
    for f = 1:12
        [ISIcount, edgesISI] = histcounts(ISIperStim{f},edgesISI);
        pdfISI(f,:) = ISIcount/sum(ISIcount);
        for r = 1:length(pdfISItotal)
            if  pdfISI(f,r) ~=0  % pdfNall(r) > 0.001 &&
                InfoISI = InfoISI + pdfISI(f,r)*log2(pdfISI(f,r)/pdfISItotal(r));
            end
        end
    end
    
    InfoISI = InfoISI/12;
    InfosISI = [InfosISI InfoISI];
end

%% plot MI

clear all
load('DataInfoStat.mat')
Infos = [ones(1,100)*InfoP;InfosN;InfosNM;InfosVSP;InfosVSN;InfosVSNM;InfosISIP;InfosISIN;InfosISINM];
Infomean = mean(Infos,2);
stdmean = [];
for i = 1:9
    stdmean(i) = std(Infos(i,:));
end
errorbar(Infomean,stdmean,'+','CapSize',30)
std(ones(1,10))
axis([0,10,0,1.2])

%%
% figure
% for f= 1:12
%     histogram(DR_trial(:,f),edges)
%     hold on
%     pause
% end
nn = max(output.spikecount_neuron{1});

VSlow = [output.VS{1}]
for i = 1:nn
    for f = 1:12
        SyncN{i}(:,f) =  output.spikecount{f}(find(output.spikecount_neuron{f}==i));
    end
    rangeN = max(max(SyncN{i}));
    pdNH{i} = zeros(1,(rangeN-mod(rangeN,10))/10 + 1);
    pdNM{i} = zeros(1,(rangeN-mod(rangeN,10))/10 + 1);
    pdNL{i} = zeros(1,(rangeN-mod(rangeN,10))/10 + 1);
    
    for j = 1:(rangeN-mod(rangeN,10))/10 + 1 
        pdNL{i}(j) = sum(sum(SyncN{i}(:,1:4)< j*10 & SyncN{i}(:,1:4)>= (j-1)*10));
        pdNM{i}(j) = sum(sum(SyncN{i}(:,5:8)< j*10 & SyncN{i}(:,5:8)>= (j-1)*10));
        pdNH{i}(j) = sum(sum(SyncN{i}(:,9:12)< j*10 & SyncN{i}(:,9:12)>= (j-1)*10));
        pdN{i}(j) = sum(sum(SyncN{i}< j*10 & SyncN{i}>= (j-1)*10));
    end
    pdNL{i} = pdNL{i}/sum(pdNL{i});
    pdNM{i} = pdNM{i}/sum(pdNM{i});
    pdNH{i} = pdNH{i}/sum(pdNH{i});
    pdN{i} = pdN{i}/sum(pdN{i});
    sumL = log2((pdNL{i}*2)./pdN{i}).*pdNL{i};
    sumL(isnan(sumL)) = 0;
    sumL(isinf(sumL)) = 0;
    sumM = log2((pdNL{i}*2)./pdN{i}).*pdNM{i};
    sumM(isnan(sumM)) = 0;
    sumM(isinf(sumM)) = 0;
    sumH = log2((pdNL{i}*2)./pdN{i}).*pdNH{i};
    sumH(isnan(sumH)) = 0;
    sumH(isinf(sumH)) = 0;
    InfoN(i) = sum(sumH)+sum(sumM)+sum(sumL);
end



load('SyncP_new.mat')
nn = max(output.spikecount_neuron{1});
for i = 1:nn
    for f = 1:12
        SyncP{i}(:,f) =  output.spikecount{f}(find(output.spikecount_neuron{f}==i));
    end
    rangeP = max(max(SyncP{i}));
    pdPH{i} = zeros(1,(rangeP-mod(rangeP,10))/10 + 1);
    pdPM{i} = zeros(1,(rangeP-mod(rangeP,10))/10 + 1);
    pdPL{i} = zeros(1,(rangeP-mod(rangeP,10))/10 + 1);
    
    for j = 1:(rangeP-mod(rangeP,10))/10 + 1
        pdPL{i}(j) = sum(sum(SyncP{i}(:,1:4)< j*10 & SyncP{i}(:,1:4)>= (j-1)*10));
        pdPM{i}(j) = sum(sum(SyncP{i}(:,5:8)< j*10 & SyncP{i}(:,5:8)>= (j-1)*10));
        pdPH{i}(j) = sum(sum(SyncP{i}(:,9:12)< j*10 & SyncP{i}(:,9:12)>= (j-1)*10));
        pdP{i}(j) = sum(sum(SyncP{i}< j*10 & SyncP{i}>= (j-1)*10));
    end
    pdPL{i} = pdPL{i}/sum(pdPL{i});
    pdPM{i} = pdPM{i}/sum(pdPM{i});
    pdPH{i} = pdPH{i}/sum(pdPH{i});
    pdP{i} = pdP{i}/sum(pdP{i});
    
    sumL = log2((pdPL{i}*2)./pdP{i}).*pdPL{i};
    sumL(isnan(sumL)) = 0;
    sumL(isinf(sumL)) = 0;
    sumM = log2((pdPM{i}*2)./pdP{i}).*pdPM{i};
    sumM(isnan(sumM)) = 0;
    sumM(isinf(sumM)) = 0;
    sumH = log2((pdPH{i}*2)./pdP{i}).*pdPH{i};
    sumH(isnan(sumH)) = 0;
    sumH(isinf(sumH)) = 0;
    InfoP(i) = sum(sumH)+sum(sumM)+sum(sumL);
end




load('SyncNM_new.mat')
nn = max(output.spikecount_neuron{1});
for i = 1:nn
    for f = 1:12
        SyncNM{i}(:,f) =  output.spikecount{f}(find(output.spikecount_neuron{f}==i));
    end
    rangeNM = max(max(SyncNM{i}));
    pdNMH{i} = zeros(1,(rangeNM-mod(rangeNM,10))/10 + 1);
    pdNMM{i} = zeros(1,(rangeNM-mod(rangeNM,10))/10 + 1);
    pdNML{i} = zeros(1,(rangeNM-mod(rangeNM,10))/10 + 1);
    
    for j = 1:(rangeNM-mod(rangeNM,10))/10 + 1
        pdNML{i}(j) = sum(sum(SyncNM{i}(:,1:4)< j*10 & SyncNM{i}(:,1:4)>= (j-1)*10));
        pdNMM{i}(j) = sum(sum(SyncNM{i}(:,5:8)< j*10 & SyncNM{i}(:,5:8)>= (j-1)*10));
        pdNMH{i}(j) = sum(sum(SyncNM{i}(:,9:12)< j*10 & SyncNM{i}(:,9:12)>= (j-1)*10));
        pdN_M{i}(j) = sum(sum(SyncNM{i}< j*10 & SyncNM{i}>= (j-1)*10));
    end
    pdNML{i} = pdNML{i}/sum(pdNML{i});
    pdNMM{i} = pdNMM{i}/sum(pdNMM{i});
    pdNMH{i} = pdNMH{i}/sum(pdNMH{i});
    pdN_M{i} = pdN_M{i}/sum(pdN_M{i});
    sumL = log2((pdNML{i}*2)./pdN_M{i}).*pdNML{i};
    sumL(isnan(sumL)) = 0;
    sumL(isinf(sumL)) = 0;
    sumM = log2((pdNMM{i}*2)./pdN_M{i}).*pdNMM{i};
    sumM(isnan(sumM)) = 0;
    sumM(isinf(sumM)) = 0;
    sumH = log2((pdNMH{i}*2)./pdN_M{i}).*pdNMH{i};
    sumH(isnan(sumH)) = 0;
    sumH(isinf(sumH)) = 0;
    InfoNM(i) = sum(sumH)+sum(sumM)+sum(sumL);
end


% NbNeurons = size(output.spikes_per_click{1, 11}.brut,1);  

%% Vector Strength Histogram 

load('SyncN_new.mat')

% load('SyncP_new.mat')


edges = 0:0.2:1;
h1 = histc(output.VS(:,2),edges); 
h2 = histc(output.VS(:,6),edges);
h3 = histc(output.VS(:,10),edges);

figure
plot(edges,h1)
hold on
plot(edges,h2)
plot(edges,h3)

%% plot spikes per click for all neurons


% load('SyncN_new.mat')
load('SyncP_new.mat')

% %48Hz adaptation
% xaxis = [];
% brut = [];
% for i = 1:16
%     xaxis = [xaxis output.spikes_per_click{1, 11}.xaxis];
%     brut = [brut output.spikes_per_click{1, 11}.brut2(i,:)];%/output.spikes_per_click{1, 11}.brut(i,1)];
% end
% 
%   
% 
% figure
% j = 1;
% for j = 1:25
%     disp(j)
%         hold on
%         plot(output.spikes_per_click{1, 11}.xaxis(j,:),output.spikes_per_click{1, 11}.brut1(j,:),'linewidth',2.0)
%         pause
% %     if brut(j) ==0;
% %     brut(j) =1;
% %     end
% end


% logbrut = log(brut);


bin = [];
f = 10;
tets = 0;
for j = size(output.spikes_per_click{1, f}.brut1,1) :-1:1
    if output.spikes_per_click{1,f}.brut1(j,1)< 8 || output.spikes_per_click{1,f}.brut1(j,2)<1
        output.spikes_per_click{1,f}.brut1(j,:) = [];
        tets = tets +1;
    end
end

    




bin(1,:) = output.spikes_per_click{1, f}.brut1(:,1);
bin(2,:) = -output.spikes_per_click{1, f}.brut1(:,1)+output.spikes_per_click{1, f}.brut1(:,2);
bin(3,:) = -output.spikes_per_click{1, f}.brut1(:,1)+output.spikes_per_click{1, f}.brut1(:,3);

bin(4,:) = output.spikes_per_click{1, f}.brut1(:,2);
bin(5,:) = -output.spikes_per_click{1, f}.brut1(:,2)+output.spikes_per_click{1, f}.brut1(:,3);

bin(6,:) = -output.spikes_per_click{1, f}.brut1(:,2)+output.spikes_per_click{1, f}.brut1(:,5);
bin(7,:) = output.spikes_per_click{1, f}.brut1(:,5);
bin(8,:) = -output.spikes_per_click{1, f}.brut1(:,5)+output.spikes_per_click{1, f}.brut1(:,8);

for j = 1:size(output.spikes_per_click{1, f}.brut1,1)

    bin(2,j) = bin(2,j)/bin(1,j);
    bin(3,j) = bin(3,j)/bin(1,j);
    bin(5,j) = bin(5,j)/bin(4,j);
    bin(6,j) = bin(6,j)/bin(4,j);
    bin(8,j) = bin(8,j)/bin(7,j);
end



C1 = cov(bin(1,:),bin(2,:));
[U,S,V] = svd(C1);
V = V*sqrt(S)*3;
v1 = V(:,1);
v2 = V(:,2);
coM(1) = mean(bin(1,:));
coM(2) = mean(bin(2,:));

C2 = cov(bin(1,:),bin(3,:));
[U2, S2, V2] = svd(C2);
V2 = V2*sqrt(S2)*3;
v3 = V2(:,1);
v4 = V2(:,1);
coM(3) = mean(bin(3,:));

% C3 = cov(bin(4,:),bin(5,:));
% C4 = corr(bin(4,:).',bin(6,:).');
% C5 = corr(bin(7,:).',bin(8,:).');



figure
subplot(2,2,[1 3])
scatter(bin(1,:),bin(2,:),'filled')
hold on
scatter(bin(4,:),bin(6,:),'filled')
scatter(bin(7,:),bin(8,:),'filled')

xlimit = get(gca,'xlim');


line([0 250], [0 0]);
line([0 250], [1 1]);
line([0 250], [-1 -1]);
axis([0,250,-5,5])

legend1 = sprintf('1st to 2nd pulse. corr = %.2f',corr(bin(1,:).',bin(2,:).'));
legend2 = sprintf('2nd to 5th pulse. corr = %.2f',corr(bin(4,:).',bin(6,:).'));
legend3 = sprintf('5th to 8th pulse. corr = %.2f',corr(bin(7,:).',bin(8,:).'));
legend({legend1, legend2, legend3});
ylabel('factor of increase in activity')
xlabel('firing rate above spont rate')
title('Evolution of firing rate in response to intial activity')
% scatter(bin(1,:),bin(7,:),'filled')
% quiver(coM(1),coM(2),v1(1),v1(2), 'b-')
% quiver(coM(1),coM(2),v2(1),v2(2), 'b-')

% scatter(bin(1,:),bin(3,:),'filled')
% % quiver(coM(1),coM(3),v3(1),v3(2), 'r-')
% % quiver(coM(1),coM(3),v4(1),v4(2), 'r-')
% 
% % hold on 



subplot(2,2,2)
edges = -5:0.20:5;
% [h1,edges] = histcounts(bin(2,:), edges);
histogram(bin(2,:), edges);
hold on 
% [h2,edges] = histcounts(bin(3,:), edges);
histogram(bin(3,:), edges);
legend1 = sprintf('1st to 2nd pulse');%. corr = %.2f',corr(bin1,bin2));
legend2 = sprintf('1st to 3rd pulse');%. corr = %.2f',corr(bin1,bin3));
% legend3 = sprintf('mu = %.3f', mean(h3));
legend({legend1, legend2}) %, legend3});

title('plasticity near onset')

subplot(2,2,4)
edges = -5:0.20:5;
% [h3,edges]= histcounts(bin(6,:), edges);
histogram(bin(6,:), edges);
hold on
% [h4,edges] = histcounts(bin(8,:), edges);
histogram(bin(8,:), edges);
legend1 = sprintf('2nd to 5th pulse');%. corr = %.2f',corr(bin1,bin2));
legend2 = sprintf('5th to 8thpulse');%. corr = %.2f',corr(bin1,bin3));
% legend3 = sprintf('mu = %.3f', mean(h3));
legend({legend1, legend2}) %, legend3});
title('plasticity near middle of stimuli')
xlabel('factor of increase in activity')


%% fitting spike per click data to a curve function 01/23/2018

load('SyncN_new.mat');
Brut = output.spikes_per_click{1, 10}.brut1  ;
xaxis = output.spikes_per_click{1, 10}.xaxis  ;
yaxis = zeros(size(xaxis));
yaxis = yaxis(:,1:end-1);
for n = 1:length(Brut)
    Brut(n,:)= Brut(n,:)/Brut(n,1);
    t=1;
    while t < size(xaxis,2)
        yaxis(n,t) = Brut(n,t+1)-1;
        t = t+1;
    end
end

xaxis = xaxis(:,2:end);
xaxis = xaxis(:);
yaxis = yaxis(:);




%% analysis of individual neurons 
% 
% 
% Analysis.pearson = zeros(NbNeurons,2);
% Analysis.spearman = zeros(NbNeurons,2);


%% Stats
% 
[p,h,stats] = ranksum(bin(6,:),bin(8,:))
% figure
% scatter(bin(2,:),bin(6,:));




%% 01/23/2018 VS comparison between first half and second half. 
load('SyncN_new.mat');
% load('SyncP_new.mat');
[p,h,stats] = ranksum(output.VStime{1, 11}(:,1),output.VStime{1, 11}(:,2))


corr(output.VStime{1, 10}(:,1),output.VStime{1, 10}(:,2))



%% 01/25/2018 Spike latency calculation, see Bendor 2008
clear all
load('SyncN_new.mat');
% load('SyncP_new.mat');


nn = size(output.mean_neuron_spont{1},1);
time = size(output.mean_neuron_spont{1},2);
for f = 1:12
    spont{f} = zeros(nn,time/2);
    rate{f} = zeros(nn,time/2);
end

for f = 1:12
    for n = 1:nn
        t = 1;
        while t < time/2 + 1
            spont{f}(n,t) = (output.mean_neuron_spont{f}(n,2*t)...
                +(output.mean_neuron_spont{f}(n,2*t-1)))/2 ;
            rate{f}(n,t) = (output.mean_neuron_rates{f}(n,2*t)...
                +(output.mean_neuron_rates{f}(n,2*t-1)))/2 ;
            t = t+1;
        end
        spontstd{f}(n) = std(spont{f}(n,:),1);
    end
end


responsetime = zeros(nn,12);
for n = 1:nn
    for f = 1:12
        t = 1;
        while t < time/2-2
            if rate{f}(n,t) > spontstd{f}(n)*2 ...
                    && rate{f}(n,t+1) > spontstd{f}(n)*2 ...
                    && rate{f}(n,t+2) > spontstd{f}(n)*2
                responsetime(n,f) = t;
                t = time;
            else
                t = t+1;
            end
        end

    end
end


responsetime = responsetime*1e-3;
% responsetime(12,:) = []; %for SyncN neurons, neuron 12 doesnt react to
% stim

% edges = 0:5:max(responsetime);
% figure
% for n = 1:nn
%     [h{n},edges] = histcounts(responsetime(n,:),edges);
%     subplot(15,2,n)
%     histogram(responsetime(n,:),edges)
%     stdNeuron(n) = std(responsetime(n,:),1);
% end


ICI_list= [250 125 83.3333 62.5 50 41.6667 35.7143 31.25 27.7778 25 22.7273 20.8333]*1e-3;

for n = 1:nn
    for f = 1:8
        firsthalf{n,f} = output.ind_spike_time{n,f}(find(responsetime(n,f) <= output.ind_spike_time{n,f} & output.ind_spike_time{n,f} <= 0.275));
        if ~isempty(firsthalf{n,f})
        firsthalf{n,f} = mod(firsthalf{n,f},ICI_list(f));
        end
        secondhalf{n,f} = output.ind_spike_time{n,f}(find(0.275 < output.ind_spike_time{n,f} & output.ind_spike_time{n,f} < 0.550));
        if ~isempty(secondhalf{n,f})
        secondhalf{n,f} = mod(secondhalf{n,f},ICI_list(f));
        end
    end
end

% perihistogram
edges = 0:0.001:0.25;
for n = 1:nn
    figure('position',[800 100 800 900])
    for f = 1:12
        subplot(12,1,f)
        histogram(firsthalf{n,f},edges)
        hold on 
        histogram(secondhalf{n,f},edges)
        
    end
    pause
end
        
    





% figure 
% boxplot(responsetime.')























