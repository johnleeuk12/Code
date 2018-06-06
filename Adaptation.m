%% Model robustness 04/06/2018

clear all
load('modeldataSN_noise.mat')

noise = 0:2:16;
% noise = 0:10;
mean_Rho = zeros(1,length(noise));
error_Rho = zeros(1,length(noise));
mean_spont = zeros(1,length(noise));
error_spont = zeros(1,length(noise));
collect1 = {};
collect2 = {};
for i = 1:length(noise)
    collect1{i} = [];
    collect2{i} = [];
    for trial = 1:10
        collect1{i} = [collect1{i} UnitInfo.Info((i-1)*10+trial).Rho];
        collect3 = [];
        for f = 1:11
            collect3 = [collect3 UnitInfo.Info((i-1)*10+trial).Output.rate{f}(1:400)];
        end
           
        collect2{i} = [collect2{i} mean(collect3,2)]; %Spont rate across stim set
    end
    mean_Rho(i) = mean(collect1{i});
    error_Rho(i) = std(collect1{i});
    mean_spont(i) = mean(collect2{i});
    error_spont(i) = std(collect2{i});
end

% for i = 2:length(noise) 
% [h,p,ci,stats] = ttest(collect{i-1},collect{i})
% pause
% end


mean_Rho = -mean_Rho;
figure
bar(noise,mean_Rho)
hold on
e = errorbar(noise,mean_Rho,error_Rho,'.','CapSize',18,'LineWidth',2);
e.Color = 'black';
e.CapSize = 18;
axis([-1 11 -1.2 1.2])

figure
errorbar(noise,mean_spont,error_spont);

figure
plot(UnitInfo.Info(30).Output.mean_discharge_rate.mean)


%% Puretone  06/06/2018

SN = load('med_spk_PT_SN.mat');
SP = load('med_spk_PT_SP.mat');
SN.median_spike = SN.median_spike-0.2;
SP.median_spike = SP.median_spike-0.2;

[p,h,stats] = ranksum(SN.median_spike,SP.median_spike)
% [p,h,stats] = ranksum(SN.median_spike,SP.median_spike)
data_mean = mean(SN.median_spike);
data_error = std(SN.median_spike);
data_mean(2) = mean(SP.median_spike);
data_error(2) = std(SP.median_spike);
figure
bar(data_mean)
hold on 
e = errorbar(data_mean,data_error,'.','CapSize',18,'LineWidth',2);
e.Color = 'black';
e.CapSize = 18;

edges = [0:0.025:0.3];
figure
histogram(SN.median_spike,edges)
hold on
histogram(SP.median_spike,edges)


% model data
load('puretoneModel.mat');
mdlSN.spiketime = UnitInfo.Info(3).Output.spiketime{1}(find(UnitInfo.Info(3).Output.spiketime{1}>0));
mdlSN.spiketime = mdlSN.spiketime(find(mdlSN.spiketime<0.2));

median(mdlSN.spiketime)

mdlSP.spiketime = UnitInfo.Info(2).Output.spiketime{1}(find(UnitInfo.Info(1).Output.spiketime{1}>0));
mdlSP.spiketime = mdlSP.spiketime(find(mdlSP.spiketime<0.2));

median(mdlSP.spiketime)

[p,h,stats] = ranksum(mdlSN.spiketime,mdlSP.spiketime)

%% Vector strength robustness 04/06/2018


clear all
load('modeldataSN_Stability.mat')

% noise = 0:2:16;
noise = 0:10;
mean_VS = zeros(length(noise),11);
error_VS = zeros(length(noise),11);
collect1 = {};
collect2 = {};
for i = 1:length(noise)
    collect1{i} = [];
    collect2{i} = [];
    for trial = 1:10
        collect1{i} = [collect1{i}; UnitInfo.Info((i-1)*10+trial).Output.VS];
        collect3 = [];
    end
    mean_VS(i,:) = mean(collect1{i},1);
    error_VS(i,:) = std(collect1{i},1);
end

% for i = 2:length(noise) 
% [h,p,ci,stats] = ttest(collect{i-1},collect{i})
% pause
% end



figure
bar3(mean_VS)
% ylabel('noise (1e-8nS)')
ylabel('jitter (ms)')
xlabel('stim (Hz)')
zlabel('VS')
% 
% set(gca,'YTickLabel',{'0' '2' '4' '6' '8' '10' '12' '14' '16'});
set(gca,'YTickLabel',{'0' '1' '2' '3' '4' '5' '6' '7' '8' '9' '10'});
set(gca,'XTickLabel',{'8' '12' '16' '20' '24' '28' '32' '36' '40' '44' '48'});
% hold on
% e = errorbar(noise,mean_Rho,error_Rho,'.','CapSize',18,'LineWidth',2);
% e.Color = 'black';
% e.CapSize = 18;
% axis([-1 11 -1.2 1.2])
% 
% figure
% errorbar(noise,mean_spont,error_spont);
% 
% figure
% plot(UnitInfo.Info(30).Output.mean_discharge_rate.mean)


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


%%  Studying FR joint distribution for SyncP SyncN and SyncM neurons. 03/19/2018

clear all

% loading and reorganizing data
load('SyncP_new.mat')
NP = size(output.rates_stim{1, 1},1);
DR_trialP = [];
for f = 1:10
    for n = 1:NP
        DR_trialP(n,f) = mean(output.rates_stim{f+1}(n,:));
    end
end

load('SyncN_new.mat')
NN = size(output.rates_stim{1, 1},1);
DR_trialN = [];
for f = 1:10
    for n = 1:NN
        DR_trialN(n,f) = mean(output.rates_stim{f+1}(n,:));
    end
end

DR_trialN(randsample(1:305,5),:) = []; %randomly cutting data size to match SyncP neurons. 

DR_trialN(find(DR_trialN>150)) = 0;
DR_trialP(find(DR_trialP>150)) = 0;
% DR_trialN(find(DR_trialN == 0)) = NaN;
% DR_trialP(find(DR_trialP == 0)) = NaN;

SyncPall = reshape(DR_trialP,[],1);
SyncNall = reshape(DR_trialN,[],1);

X = [SyncNall,SyncPall];
edges = {0:2:ceil(max(max(X))) 0:2:ceil(max(max(X)))};
hist3(X,edges)
[CountData, Xedges,Yedges] = histcounts2(X(:,1),X(:,2),edges{1},edges{2});
% linear interpolarization of Distribution
% Vqall = interp2(interp2(CountData,'cubic'),'cubic');
Vqall = interp2(CountData,'cubic');
Vqall(find(Vqall<0)) = 0;
pdfall = Vqall/sum(sum(Vqall));
figure
mesh(Vqall)



% Do the same thing for all frequencies

Countdatas = {};
Vq = {};
pdfs = {};
for f = 1 :10
    Countdatas{f} = histcounts2(DR_trialN(:,f),DR_trialP(:,f),edges{1},edges{2});
    %     Vq{f} = interp2(interp2(Countdatas{f},'cubic'),'cubic');
    Vq{f} = interp2(Countdatas{f},'cubic');
    Vq{f}(find(Vq{f}<0)) = 0;
    pdfs{f} = Vq{f}/sum(sum(Vq{f}));
end

% Calculating Information content

InfoPN = 0;

for f = 1:10
    for r1 = 1:length(pdfall)
        for r2 = 1:length(pdfall)
            if pdfs{f}(r1,r2) ~= 0 && pdfall(r1,r2) ~=0
                InfoPN = InfoPN + pdfs{f}(r1,r2)*log2(pdfs{f}(r1,r2)/pdfall(r1,r2));
            end
        end
    end
end

InfoPN = InfoPN/10;
log2(10)


load('SyncNM_new.mat')
NNM = size(output.rates_stim{1, 1},1);
DR_trialNM1 = zeros(300,12);
DR_trialNM2 = zeros(300,12);
NM1 = randsample(NNM,600);
NM2 = NM1(301:600);
NM1 = NM1(1:300);

for f = 1:10
    for n = 1:300
        DR_trialNM1(n,f) = mean(output.rates_stim{f+1}(NM1(n),:));
        DR_trialNM2(n,f) = mean(output.rates_stim{f+1}(NM2(n),:));
    end
end


% DR_trialNM1(find(DR_trialNM1 == 0)) = NaN;
% DR_trialNM2(find(DR_trialNM2 == 0)) = NaN;
% DR_trialNM1(find(DR_trialNM1>150)) = 0;
% DR_trialNM2(find(DR_trialNM2>150)) = 0;



SyncNM1all = reshape(DR_trialNM1,[],1);
SyncNM2all = reshape(DR_trialNM2,[],1);

XNM = [SyncNM1all,SyncNM2all];
edges = {0:2:ceil(max(max(XNM))) 0:2:ceil(max(max(XNM)))};
hist3(XNM,edges)
[CountDataNM, Xedges,Yedges] = histcounts2(XNM(:,1),XNM(:,2),edges{1},edges{2});

% linear interpolarization of Distribution
% Vqall = interp2(interp2(CountData,'cubic'),'cubic');
VqallNM = interp2(CountDataNM,'cubic');
VqallNM(find(VqallNM<0)) = 0;
pdfallNM = VqallNM/sum(sum(VqallNM));
figure
mesh(VqallNM)


% Do the same thing for all frequencies

CountdatasNM = {};
VqNM = {};
pdfsNM = {};
for f = 1 :10
    CountdatasNM{f} = histcounts2(DR_trialNM1(:,f),DR_trialNM2(:,f),edges{1},edges{2});
    %     Vq{f} = interp2(interp2(Countdatas{f},'cubic'),'cubic');
    VqNM{f} = interp2(Countdatas{f},'cubic');
    VqNM{f}(find(VqNM{f}<0)) = 0;
    pdfsNM{f} = VqNM{f}/sum(sum(VqNM{f}));
end

% Calculating Information content

InfoNM = 0;

for f = 1:10
    for r1 = 1:length(pdfallNM)
        for r2 = 1:length(pdfallNM)
            if pdfsNM{f}(r1,r2) ~= 0 && pdfallNM(r1,r2) ~=0
                InfoNM = InfoNM + pdfsNM{f}(r1,r2)*log2(pdfsNM{f}(r1,r2)/pdfallNM(r1,r2));
            end
        end
    end
end

InfoNM = InfoNM/10;


%% Calculating MI 07/11/2017

% Recalculating MI using a more simple method, directly comparable to joint
% MI


clear all
load('SyncP_new.mat')
% NP = size(output.rates_stim{1, 1},1);
NP = 300;
DR_trialP = [];
for f = 1:10
    for n = 1:NP
        DR_trialP(n,f) = mean(output.rates_stim{f+1}(n,:));
    end
end

DR_trialP(find(DR_trialP>100)) = 0;

% 08/11/2017 attempt to linearly extrapolate a pdf.
edgesTotalP = 0:2: ceil(max(max(DR_trialP)));
[NtotalP,edgesTotalP] = histcounts(DR_trialP,edgesTotalP);
% for i = 1:length(NtotalP)
%     if NtotalP(i) == 0
%         NtotalP(i) = NaN;
%     end
% end

histogram(DR_trialP,edgesTotalP)

plot(NtotalP)
xq = 0:1: ceil(max(max(DR_trialP)));
% linear interpolation into infinite dataset
vq1 = interp1(edgesTotalP(2:end),NtotalP,xq(2:end),'pchip');
vq1(find(vq1 <0)) = 0;
vq1 = [vq1(2:end) 0];
pdfPall = vq1/sum(vq1);
% vq2 = interp1(xData, yData, xaxis,'spline');
% plot(xaxis,vq1)
% hold on 
% plot(xaxis,vq2)
% vq2 = interp1(edgesTotalP(2:end),NtotalP,xq(2:end),'pchip');
% plot(vq2)


%Gaussian kernel based estimation. 
edgesTotalPgauss = 0:2: ceil(max(max(DR_trialP)));
[NtotalPgauss,edgesTotalPgauss] = histcounts(DR_trialP,edgesTotalPgauss);
xs = edgesTotalPgauss(1:end-1)+0.5;
h = 1;

% [vqKS,xi,bw] =
% ksdensity(NtotalPgauss,xs,'Bandwidth',1,'BoundaryCorrection','reflection');
% Offcial kernel density estimation, but I don't know how to use it yet.

for i = 1:length(xs)
    ys(i) = gaussian_kern_reg2(xs(i),xs,NtotalPgauss,h);
end
vqgauss = ys/sum(ys);

%Compare
% figure
% plot(xs,vqgauss)
% hold on
% plot(pdfPall)
% % plot(xs,vqKS)


pdfPall = vqgauss;

InfoP = 0;
figure
for f= 1:10
    %linear interpolation
    %     [NcountP,edgesTotalP] = histcounts(DR_trialP(:,f),edgesTotalP);
    % %     NcountP(find(NcountP==0)) = NaN;
    %     vq2N = interp1(edgesTotalP(2:end),NcountP,xq(2:end),'pchip');
    
    % gaussian kernel
    [NcountPgauss,edgesTotalPgauss] = histcounts(DR_trialP(:,f),edgesTotalPgauss);
    for i = 1:length(xs)
        vq2N(i) = gaussian_kern_reg(xs(i),xs,NcountPgauss,h);
    end
    
    %below is common 
    vq2N(find(vq2N <0)) = 0;
    vq2N= [vq2N(2:end) 0];
    pdfP(f,:) = vq2N/sum(vq2N);
    plot(pdfP(f,:))
    hold on
    for r = 1:length(pdfPall)
        if pdfP(f,r) ~=0 && pdfPall(r) ~=0
            InfoP = InfoP + pdfP(f,r)*log2(pdfP(f,r)/pdfPall(r));
        end
    end
end
InfoP = InfoP/12;








% 
% clear all
% load('SyncP_new.mat')
% NP = size(output.rates_stim{1, 1},1);
% DR_trialP = [];
% for f = 1:12
%     for n = 1:NP
%         DR_trialP(n,f) = mean(output.rates_stim{f}(n,:));
%     end
% end
% 
% DR_trialP(find(DR_trialP>100)) = 0;
% 
% % 08/11/2017 attempt to linearly extrapolate a pdf.
% edgesTotalP = 0:1: ceil(max(max(DR_trialP)));
% [NtotalP,edgesTotalP] = histcounts(DR_trialP,edgesTotalP);
% for i = 1:length(NtotalP)
%     if NtotalP(i) == 0
%         NtotalP(i) = NaN;
%     end
% end
% xaxisP = edgesTotalP(1:end-1);
% % plot(xaxis,Ntotal)
% % linear interpolation into infinite dataset
% [xDataP, yDataP] = prepareCurveData( xaxisP, NtotalP );
% vq1 = interp1(xDataP, yDataP, xaxisP);
% pdfPall = vq1/sum(vq1);
% % vq2 = interp1(xData, yData, xaxis,'spline');
% % plot(xaxis,vq1)
% % hold on 
% % plot(xaxis,vq2)
% 
% 
% InfoP = 0;
% for f= 1:12
%     [NcountP,edgesTotalP] = histcounts(DR_trialP(:,f),edgesTotalP);
%     [xDataP, yDataP] = prepareCurveData( xaxisP, NcountP );
%     vq2N = interp1(xDataP, yDataP, xaxisP);
%     pdfP(f,:) = vq2N/sum(vq2N);
%     for r = 1:length(pdfPall)
%         if pdfP(f,r) ~=0 %  pdfPall(r) > 0.001 &&
%         InfoP = InfoP + pdfP(f,r)*log2(pdfP(f,r)/pdfPall(r));
%         end
%     end
% end
% InfoP = InfoP/12;
% 
% 
% load('SyncN_new.mat')
% 
% %for NM neurons, randomly select 22neurons
% 
% 
% 
% N = size(output.rates_stim{1, 1},1);
% DR_trialN = [];
% Nall = 300;
% InfosN = [];
% for iii = 1:100
%     pdfN = [];
%     Indice = randsample([1:N],Nall);
%     
%     for f = 1:12
%         for n = 1:Nall
%             DR_trialN(n,f) = mean(output.rates_stim{f}(Indice(n),:));
%         end
%     end
%     
%     % 08/11/2017 attempt to linearly extrapolate a pdf.
%     edgesTotalN = 0:1: ceil(max(max(DR_trialN)));
%     [NtotalN,edgesTotalN] = histcounts(DR_trialN,edgesTotalN);
%     for i = 1:length(NtotalN)
%         if NtotalN(i) == 0
%             NtotalN(i) = NaN;
%         end
%     end
%     xaxisN = edgesTotalN(1:end-1);
%     % plot(xaxis,Ntotal)
%     % linear interpolation into infinite dataset
%     [xDataN, yDataN] = prepareCurveData( xaxisN, NtotalN );
%     vq1 = interp1(xDataN, yDataN, xaxisN);
%     pdfNall = vq1/sum(vq1);
%     % vq2 = interp1(xData, yData, xaxis,'spline');
%     % plot(xaxis,vq1)
%     % hold on
%     % plot(xaxis,vq2)
%     
%     
%     InfoN = 0;
%     for f= 1:12
%         [NcountN,edgesTotalN] = histcounts(DR_trialN(:,f),edgesTotalN);
%         [xDataN, yDataN] = prepareCurveData( xaxisN, NcountN );
%         vq2N = interp1(xDataN, yDataN, xaxisN);
%         pdfN(f,:) = vq2N/sum(vq2N);
%         for r = 1:length(pdfNall)
%             if  pdfN(f,r) ~=0  % pdfNall(r) > 0.001 &&
%                 InfoN = InfoN + pdfN(f,r)*log2(pdfN(f,r)/pdfNall(r));
%             end
%         end
%     end
%     InfoN = InfoN/12;
%     
% %     InfoPandN = 0;
% %     for f = 1:12
% %         %     [NcountN,edgesTotalN] = histcounts(DR_trialN(:,f),edgesTotalN);
% %         %     [xDataN, yDataN] = prepareCurveData( xaxisN, NcountN );
% %         %     vq2N = interp1(xDataN, yDataN, xaxisN);
% %         %     pdfN(f,:) = vq2N/sum(vq2N);
% %         %     [NcountP,edgesTotalP] = histcounts(DR_trialP(:,f),edgesTotalP);
% %         %     [xDataP, yDataP] = prepareCurveData( xaxisP, NcountP );
% %         %     vq2P = interp1(xDataP, yDataP, xaxisP);
% %         %     pdfP(f,:) = vq2P/sum(vq2P);
% %         for r1 = 1:length(pdfPall)
% %             for r2 = 1:length(pdfNall)
% %                 if   pdfP(f,r1)*pdfN(f,r2) ~=0 %pdfPall(r1)*pdfNall(r2) > 0.001 &&
% %                     InfoPandN = InfoPandN + ...
% %                         pdfP(f,r1)*pdfN(f,r2)*log2( pdfP(f,r1)*pdfN(f,r2)...
% %                         /(pdfPall(r1)*pdfNall(r2)));
% %                 end
% %             end
% %         end
% %     end
% %     
% %     InfoPandN = InfoPandN/12;
%     InfosN = [InfosN InfoN];
% end
% % 08/11/2017 attempt to calculate pdf by histograms
% % % test = mean(DR_trial,1);
% % % edges = 0:1: floor(max(max(DR_trial)));
% % Nbins = 219;
% % % Ntotal = histogram(DR_trial,edges);
% % [Ntotal,edgesTotal] = histcounts(DR_trial,Nbins);
% % 
% % pdfPall =  Ntotal/sum(Ntotal);
% % InfoP = 0;
% % for f= 1:12
% %     [Ncount,edges] = histcounts(DR_trial(:,f),Nbins);
% %     pdfP(f,:) = Ncount/sum(Ncount);
% %     for r = 1:length(pdfPall)
% %         if pdfPall(r) > 0.001 && pdfP(f,r) ~=0 
% %         InfoP = InfoP + pdfP(f,r)*log2(pdfP(f,r)/pdfPall(r));
% %         end
% %     end
% % end
% % % Ntotal = histogram(DR_trial,edges);
% % 

%% Analyzing ind ISI 19/03/2018

clear all

load('SyncP_new.mat')
ISI_ind = {};
for n = 1:25
    for f = 1:12
        ISI_ind{n,f} = horzcat(output.isi_rep{n}{f,:})*1e3;
    end
    
end
C = {};

% edges = [[0:1:9] [10:5:90] [100:50:1000]];
edges = [0:5:500];
for n = 1:25
    figure
    for f = 2:2:12
        [C{n,f}, edges] = histcounts(ISI_ind{n,f},edges);
        plot(edges(2:end),C{n,f})
        set(gca,'XScale','log')
        hold on
    end
    legend('125ms', '62.5ms', '41.6ms','31.25ms', '25ms', '20ms');
    hold off
end

%conclusion, ISI analysis really doesnt give much



    
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
    stdmean(i) = deviation(Infos(i,:));
end
errorbar(Infomean,stdmean,'+','CapSize',30)
deviation(ones(1,10))
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



%% 07/03/2018 
clear all
load('SyncP_new.mat')
ICI_list= [125 83.3333 62.5 50 41.6667 35.7143 31.25 27.7778 25 22.7273 20.8333];
Hz_list = round(1./ICI_list*1e3);

SPC_mean = []; %average spikes per click
SPC_error = [];

for f = 1:length(ICI_list)
    SPC_mean = [SPC_mean mean(output.spikes_per_click{f}.mean)];
    SPC_error = [SPC_error mean(output.spikes_per_click{f}.error)];
end

[RHO,PVAL] = corr(ICI_list.',SPC_mean.','Type','Spearman')


errorbar(Hz_list, SPC_mean, SPC_error,'Linewidth',1.5)
hold on     

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
% clear all
% load('SyncN_new.mat');
load('SyncP_new.mat');


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

edges = 0:5*1e-3:max(responsetime);
figure
for n = 1:nn
    [h{n},edges] = histcounts(responsetime(n,:),edges);
    subplot(15,2,n)
    histogram(responsetime(n,:),edges)
    stdNeuron(n) = std(responsetime(n,:),1);
end
figure

edges = 0:2*1e-3:max(responsetime);
histogram(responsetime,edges)

ICI_list= [250 125 83.3333 62.5 50 41.6667 35.7143 31.25 27.7778 25 22.7273 20.8333]*1e-3;

for n = 1:nn
    for f = 1:12
        firsthalf{n,f} = output.ind_spike_time{n,f}(find(responsetime(n,f) <= output.ind_spike_time{n,f} & output.ind_spike_time{n,f} <= 0.275));
        if ~isempty(firsthalf{n,f})
        firsthalf{n,f} = 2*pi*mod(firsthalf{n,f},ICI_list(f))/ICI_list(f);
        end
        secondhalf{n,f} = output.ind_spike_time{n,f}(find(0.275 < output.ind_spike_time{n,f} & output.ind_spike_time{n,f} < 0.550));
        if ~isempty(secondhalf{n,f})
        secondhalf{n,f} = 2*pi*mod(secondhalf{n,f},ICI_list(f))/ICI_list(f);
        end
    end
end

% perihistogram
% edges = 0:0.001:0.25;
% edges = 0: pi/32: 2*pi;
% for n = 1:nn
%     figure('position',[800 100 800 900])
%     for f = 1:12
%         subplot(12,1,f)
%         
%         histogram(firsthalf{n,f},edges) %Blue
% %         pause
%         hold on
%         histogram(secondhalf{n,f},edges) %Orange
%         axis([0,2*pi,0,20])
%         if f <12
%             set(gca,'XTick',[]);
%         else
%             xticks([0 0.5*pi pi 1.5*pi 2*pi])
%             xticklabels({'0' '0.5\pi' '\pi' '1.5\pi' '2\pi'})
%         end
%         
%     end
% %     pause
% end
% 


figure 
boxplot(responsetime.')


% 29/03/2018    subtracting spike latency to spike time to pool all neurons

spiketime = {};
for f = 1:12
    spiketime{f} = [];
    for n = 1:nn
        spiketime{f} = [spiketime{f} output.ind_spike_time{n,f}-responsetime(n,f)];
    end
end


for f = 1:12
    firsthalf{f} = spiketime{f}(find(responsetime(n,f) <= spiketime{f} & spiketime{f} <= 0.275));
    if ~isempty(firsthalf{f})
        firsthalf{f} = 2*pi*mod(firsthalf{f},ICI_list(f))/ICI_list(f);
    end
    secondhalf{f} = spiketime{f}(find(0.275 < spiketime{f} & spiketime{f} < 0.550));
    if ~isempty(secondhalf{f})
        secondhalf{f} = 2*pi*mod(secondhalf{f},ICI_list(f))/ICI_list(f);
    end
end

edges = 0: pi/32: 2*pi;
% edges = 0:0.002:0.25;

 figure('position',[800 100 800 900])
 for f = 1:12
     subplot(12,1,f)
     
     histogram(firsthalf{f},edges) %Blue
     %         pause
     hold on
     histogram(secondhalf{f},edges) %Orange
     axis([0,2*pi,0,inf])
     if f <12
         set(gca,'XTick',[]);
%      else
         xticks([0 0.5*pi pi 1.5*pi 2*pi])
         xticklabels({'0' '0.5\pi' '\pi' '1.5\pi' '2\pi'})
     end
     
 end


%% Model data analysis and prediction 01/02/2018
clear all
load('SyncPModel_3.mat')
% load('SyncNModel_3.mat')



ICI_list= [125 83.3333 62.5 50 41.6667 35.7143 31.25 27.7778 25 22.7273 20.8333]*1e-3;
for f = 1:12
firsthalf{f} = [];
secondhalf{f} = [];
end
for trial = 1:12
    if abs(UnitInfo.Info(trial).Rho) > 0.8 && UnitInfo.Info(trial).Pval <0.05
        output = UnitInfo.Info(trial).Output;
        for f = 1:11
            firsthalf{f} = [firsthalf{f} output.spiketime{f}(find(0.05 <= output.spiketime{f} & output.spiketime{f} <= 0.275))];           
            secondhalf{f} = [secondhalf{f} output.spiketime{f}(find(0.275 < output.spiketime{f} & output.spiketime{f} < 0.550))];
        end
        
    end
end

for f = 1:11
    if ~isempty(firsthalf{f})
        firsthalf{f} = 2*pi*mod(firsthalf{f},ICI_list(f))/ICI_list(f);
    end
    if ~isempty(secondhalf{f})
        secondhalf{f} = 2*pi*mod(secondhalf{f},ICI_list(f))/ICI_list(f);
    end
end

edges = 0: pi/64: 2*pi;
n = 1;
figure('position',[800 100 800 900])
for f = 1:11
    subplot(12,1,f)
    
    histogram(firsthalf{n,f},edges) %Blue
    %         pause
    hold on
    histogram(secondhalf{n,f},edges) %Orange
%     axis([0,2*pi,0,40])
    if f <12
        set(gca,'XTick',[]);
    else
        xticks([0 0.5*pi pi 1.5*pi 2*pi])
        xticklabels({'0' '0.5\pi' '\pi' '1.5\pi' '2\pi'})
    end
    
end

medians = zeros(2,11);
deviation = zeros(2,11);
for f = 1:11
    [p,h,stats] = ranksum(firsthalf{f},secondhalf{f});
    if h ==1
    medians(1,f) = median(firsthalf{f});
    medians(2,f) = median(secondhalf{f});
    deviation(1,f) = std(firsthalf{f});
    deviation(2,f) = std(secondhalf{f});
    end
end





%% New model data including facilitation. 

load('dataaaaa2.mat')

% UnitInfo.List : {f_dE, f_fI, tau_pE, tau_pI, E, I}
% E and I stays constant.

X1 = zeros(length(unique(UnitInfo.List(:,1))),length(unique(UnitInfo.List(:,2)))); % Matrix of Rho depending on adaptation values
X3 = X1;
for n = 1:length(UnitInfo.List)
    if UnitInfo.Info(n).Pval < 0.05
        
        X1(floor(UnitInfo.List(n,1)*10),floor(UnitInfo.List(n,2)*10)) = X1(floor(UnitInfo.List(n,1)*10),floor(UnitInfo.List(n,2)*10)) + UnitInfo.Info(n).Rho;
    end
    X3(floor(UnitInfo.List(n,1)*10),floor(UnitInfo.List(n,2)*10)) = X3(floor(UnitInfo.List(n,1)*10),floor(UnitInfo.List(n,2)*10)) + UnitInfo.Info(n).Output.mean_discharge_rate.mean(11);
end

X1 = X1/(8*8);
X3 = X3/64;

figure
imagesc(X1)
xticks([1:9])
xticklabels({'0.1','0.2','0.3','0.4','0.5','0.6','0.7','0.8','0.9'})
xlabel('f_DI')
yticks([1:3])
yticklabels({'0.1','0.2','0.3'})
ylabel('f_DE')
figure
imagesc(X3)
xticks([1:9])
xticklabels({'0.1','0.2','0.3','0.4','0.5','0.6','0.7','0.8','0.9'})
xlabel('f_DI')
yticks([1:3])
yticklabels({'0.1','0.2','0.3'})
ylabel('f_DE')

X2 = {};
figure
i = 1;
f_dE = 0.3;
tau_pE = unique(UnitInfo.List(:,3));
tau_pI = unique(UnitInfo.List(:,4));
for f_dI = 0.1:0.1:0.9
    subplot(2,5,i)
    X2{i} = zeros(8,8);
    
    x = find(UnitInfo.List(:,1)==f_dE & UnitInfo.List(:,2)==f_dI);
    for n = 1:length(x)
        if UnitInfo.Info(x(n)).Pval < 0.05
        X2{i}(find(tau_pE == UnitInfo.List(x(n),3)),tau_pI == UnitInfo.List(x(n),4)) =  X2{i}(find(tau_pE == UnitInfo.List(x(n),3)),tau_pI == UnitInfo.List(x(n),4))+ UnitInfo.Info(x(n)).Rho;
        end
    end
    imagesc(X2{i})
    colorbar
    caxis([-1 1])
    i = i+1;
end

 figure
subplot(1,3,2)
f_dE = 0.2;
f_dI = 0.9;
X2{1} = zeros(8,8);

x = find(round(UnitInfo.List(:,1)*10)==f_dE*10 & round(UnitInfo.List(:,2)*10)==f_dI*10);
for n = 1:length(x)
    if UnitInfo.Info(x(n)).Pval < 0.05
        
        X2{1}(find(tau_pE == UnitInfo.List(x(n),3)),tau_pI == UnitInfo.List(x(n),4)) =  X2{1}(find(tau_pE == UnitInfo.List(x(n),3)),tau_pI == UnitInfo.List(x(n),4))+ UnitInfo.Info(x(n)).Rho;
    end
end
imagesc(X2{1})
colorbar
caxis([-1 1])




%% spikes per click in model neurons. 

load('SyncPfMdl.mat')
ICI_list = [125 83.3333 62.5 50 41.6667 35.7143 31.25 27.7778 25 22.7273 20.8333];
Hz_list = [];
for i = 1:length(ICI_list)
    Hz_list = [Hz_list round(1000/ICI_list(i))];
end

figure
cmap = colormap(jet(length(ICI_list)+1));

for n =1:2:length(ICI_list)
    SEM = UnitInfo.Info.Output.spikes_per_click{n}.std/sqrt(20*600);
    ts = tinv([0.025  0.975],20*(600)-1);   
    error = ts(2)*SEM;
    errorbar(UnitInfo.Info.Output.spikes_per_click{n}.xaxis,UnitInfo.Info.Output.spikes_per_click{n}.mean,...
        error, 'linewidth', 1.7,'color',cmap(n+1,:),'DisplayName', ...
        [num2str(Hz_list(n)) 'Hz'])
    hold on
    %         norm_mean = [norm_mean UnitInfo.Info(p).Output.mean_discharge_rate.mean(n)/Hz_list(n)];
end
axis([500,1000,0,90]);
legend('show')





















