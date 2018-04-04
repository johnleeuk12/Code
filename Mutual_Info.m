%% Mutual information scripts

function Mutual_info()


clear all

bw = 1; %0.6:0.2:4;
% Moyen = zeros(2,length(bw));
% Erreurs =zeros(2,length(bw));
% parfor i = 1:length(bw)
%     
%     Out = MI_joint(bw(i))
%     Moyen(:,i) = Out.mean; 
%     Erreurs(:,i) = Out.std; 
% end
% 
% Output.mean = Moyen;
% Output.errors = Erreurs;

% save('MI_data.mat', 'Output')
% 
% 
% figure
% errorbar(bw,Output.mean(1,:),Output.errors(1,:))
% hold on
% errorbar(bw,Output.mean(2,:),Output.errors(2,:))


Ind_mean = zeros(3,length(bw));
Ind_error = zeros(3,length(bw));
for i = 1:length(bw)
    for d = 1:6
        Out = MI_ind(bw(i),d,10);
        Ind_mean(d,i) = Out.mean;
        Ind_error(d,i) = Out.error;
    end
end

figure
for d = 1:6
errorbar(bw,Ind_mean(d,:),Ind_error(d,:));
hold on
end




function Out = MI_joint(bw)
%%  Studying FR joint distribution for SyncP SyncN and SyncM neurons. 03/19/2018

% clear all
% bw = 1;
MI_NM  = [];
MI_PN = [];
% loading and reorganizing data
for trial = 1:100
    disp(trial)
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

% Kernel density estimation
[x1,x2] = meshgrid(0:1:ceil(max(max(X))),0:1:ceil(max(max(X))));
x1 = x1(:);
x2 = x2(:);
xi = [x1 x2];

% bw = 5;
vqKS = ksdensity(X,xi,'Bandwidth',bw);

pdfall = vqKS/sum(vqKS);

% Calculating Information content

InfoPN = 0;
Vq = {};
pdfs = {};


for f = 1:10
    % kernel density estimation
    Vq{f} = ksdensity([DR_trialN(:,f) DR_trialP(:,f)],xi,'Bandwidth',bw);
    Vq{f}(find(Vq{f}<0)) = 0;
    pdfs{f} = Vq{f}/sum(sum(Vq{f}));
    for r1 = 1:length(pdfall)
        if pdfs{f}(r1) ~= 0 && pdfall(r1) ~=0
            InfoPN = InfoPN + pdfs{f}(r1)*log2(pdfs{f}(r1)/pdfall(r1));
        end
        
    end
end

InfoPN = InfoPN/10;

MI_PN = [MI_PN InfoPN];

load('SyncNM_new.mat')
NNM = size(output.rates_stim{1, 1},1);
DR_trialNM1 = zeros(300,10);
DR_trialNM2 = zeros(300,10);
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
DR_trialNM1(find(DR_trialNM1>150)) = 0;
DR_trialNM2(find(DR_trialNM2>150)) = 0;



SyncNM1all = reshape(DR_trialNM1,[],1);
SyncNM2all = reshape(DR_trialNM2,[],1);

XNM = [SyncNM1all,SyncNM2all];


[x1,x2] = meshgrid(0:1:ceil(max(max(XNM))),0:1:ceil(max(max(XNM))));
x1 = x1(:);
x2 = x2(:);
xi = [x1 x2];
vqKS = ksdensity(XNM,xi,'Bandwidth',bw);
pdfallNM = vqKS/sum(vqKS);

VqNM = {};
pdfsNM = {};
InfoNM = 0;

% kernel density
for f = 1:10
    VqNM{f} = ksdensity([DR_trialNM1(:,f) DR_trialNM2(:,f)],xi,'Bandwidth',bw);
    VqNM{f}(find(VqNM{f}<0)) = 0;
    pdfsNM{f} = VqNM{f}/sum(sum(VqNM{f}));
    for r1 = 1:length(pdfallNM)
        if pdfsNM{f}(r1) ~= 0 && pdfallNM(r1) ~=0
            InfoNM = InfoNM + pdfsNM{f}(r1)*log2(pdfsNM{f}(r1)/pdfallNM(r1));
        end
    end
end
InfoNM = InfoNM/10;


MI_NM = [MI_NM InfoNM];
end



Out.mean = [mean(MI_PN); mean(MI_NM)];
Out.std = [std(MI_PN); std(MI_NM)];



% previous code 28/03/2018
% 
% clear all
% 
% % loading and reorganizing data
% load('SyncP_new.mat')
% NP = size(output.rates_stim{1, 1},1);
% DR_trialP = [];
% for f = 1:10
%     for n = 1:NP
%         DR_trialP(n,f) = mean(output.rates_stim{f+1}(n,:));
%     end
% end
% 
% load('SyncN_new.mat')
% NN = size(output.rates_stim{1, 1},1);
% DR_trialN = [];
% for f = 1:10
%     for n = 1:NN
%         DR_trialN(n,f) = mean(output.rates_stim{f+1}(n,:));
%     end
% end
% 
% DR_trialN(randsample(1:305,5),:) = []; %randomly cutting data size to match SyncP neurons. 
% 
% DR_trialN(find(DR_trialN>150)) = 0;
% DR_trialP(find(DR_trialP>150)) = 0;
% % DR_trialN(find(DR_trialN == 0)) = NaN;
% % DR_trialP(find(DR_trialP == 0)) = NaN;
% 
% SyncPall = reshape(DR_trialP,[],1);
% SyncNall = reshape(DR_trialN,[],1);
% 
% X = [SyncNall,SyncPall];
% % edges = {0:2:ceil(max(max(X))) 0:2:ceil(max(max(X)))};
% % hist3(X,edges)
% % [CountData, Xedges,Yedges] = histcounts2(X(:,1),X(:,2),edges{1},edges{2});
% % linear interpolarization of Distribution
% % Vqall = interp2(interp2(CountData,'cubic'),'cubic');
% % Vqall = interp2(CountData,'cubic');
% % Vqall(find(Vqall<0)) = 0;
% % pdfall = Vqall/sum(sum(Vqall));
% % figure
% % mesh(Vqall)
% 
% 
% % Kernel density estimation
% [x1,x2] = meshgrid(0:1:ceil(max(max(X))),0:1:ceil(max(max(X))));
% x1 = x1(:);
% x2 = x2(:);
% xi = [x1 x2];
% 
% bw = 5;
% vqKS = ksdensity(X,xi,'Bandwidth',bw);
% 
% pdfall = vqKS/sum(vqKS);
% 
% 
% % %Gaussian regression
% % [x1, x2] = meshgrid(Xedges(1:end-1),Yedges(1:end-1));
% % x = [x1(:).'; x2(:).'];
% % % x1 = x1(:).';
% % % x2 = x2(:).';
% % xi = [x1; x2];
% % y = CountData(:).';
% % for i = 1:size(x1,1)
% %     for j = 1:size(x2,2)
% %         xs = [x1(i,j); x2(i,j)];
% %         Vqall(i,j) = gaussian_kern_reg2(xs,x,y,[1; 1]);
% %     end
% % end
% % Vqall(find(Vqall<0)) = 0;
% % pdfall = Vqall/sum(sum(Vqall)) ;
% 
% 
% 
% 
% % Do the same thing for all frequencies
% 
% Countdatas = {};
% Vq = {};
% pdfs = {};
% 
% % Calculating Information content
% 
% InfoPN = 0;
% % 
% % for f = 1:10
% %     Countdatas{f} = histcounts2(DR_trialN(:,f),DR_trialP(:,f),edges{1},edges{2});
% %     Vq{f} = interp2(interp2(Countdatas{f},'cubic'),'cubic');
% %     Vq{f} = interp2(Countdatas{f},'cubic');
% %     Vq{f}(find(Vq{f}<0)) = 0;
% %     pdfs{f} = Vq{f}/sum(sum(Vq{f}));
% %     
% %     for r1 = 1:length(pdfall)
% %         for r2 = 1:length(pdfall)
% %             if pdfs{f}(r1,r2) ~= 0 && pdfall(r1,r2) ~=0
% %                 InfoPN = InfoPN + pdfs{f}(r1,r2)*log2(pdfs{f}(r1,r2)/pdfall(r1,r2));
% %             end
% %         end
% %     end
% % end
% % 
% % for f = 1:10
% %     Countdatas{f} = histcounts2(DR_trialN(:,f),DR_trialP(:,f),edges{1},edges{2});
% %     y = Countdatas{f}(:).';
% %     for i = 1:size(x1,1)
% %         for j = 1:size(x2,2)
% %             xs = [x1(i,j); x2(i,j)];
% %             Vq{f}(i,j) = gaussian_kern_reg2(xs,x,y,[1; 1]);
% %         end
% %     end
% %     
% %     VqNM{f}(find(Vq{f}<0)) = 0;
% %     pdfs{f} = Vq{f}/sum(sum(Vq{f}));
% %     for r1 = 1:length(pdfall)
% %         for r2 = 1:length(pdfall)
% %             if pdfs{f}(r1,r2) ~= 0 && pdfall(r1,r2) ~=0
% %                 InfoPN = InfoPN + pdfs{f}(r1,r2)*log2(pdfs{f}(r1,r2)/pdfall(r1,r2));
% %             end
% %         end
% %     end
% %     
% % end
% 
% 
% 
% for f = 1:10
%     % kernel density estimation
%     Vq{f} = ksdensity([DR_trialN(:,f) DR_trialP(:,f)],xi,'Bandwidth',bw);
%     Vq{f}(find(Vq{f}<0)) = 0;
%     pdfs{f} = Vq{f}/sum(sum(Vq{f}));
%     for r1 = 1:length(pdfall)
%         if pdfs{f}(r1) ~= 0 && pdfall(r1) ~=0
%             InfoPN = InfoPN + pdfs{f}(r1)*log2(pdfs{f}(r1)/pdfall(r1));
%         end
%         
%     end
% end
% 
% InfoPN = InfoPN/10;
% log2(10)
% 
% 
% load('SyncNM_new.mat')
% NNM = size(output.rates_stim{1, 1},1);
% DR_trialNM1 = zeros(300,10);
% DR_trialNM2 = zeros(300,10);
% NM1 = randsample(NNM,600);
% NM2 = NM1(301:600);
% NM1 = NM1(1:300);
% 
% for f = 1:10
%     for n = 1:300
%         DR_trialNM1(n,f) = mean(output.rates_stim{f+1}(NM1(n),:));
%         DR_trialNM2(n,f) = mean(output.rates_stim{f+1}(NM2(n),:));
%     end
% end
% 
% 
% % DR_trialNM1(find(DR_trialNM1 == 0)) = NaN;
% % DR_trialNM2(find(DR_trialNM2 == 0)) = NaN;
% DR_trialNM1(find(DR_trialNM1>150)) = 0;
% DR_trialNM2(find(DR_trialNM2>150)) = 0;
% 
% 
% 
% SyncNM1all = reshape(DR_trialNM1,[],1);
% SyncNM2all = reshape(DR_trialNM2,[],1);
% 
% XNM = [SyncNM1all,SyncNM2all];
% edges = {0:2:ceil(max(max(XNM))) 0:2:ceil(max(max(XNM)))};
% % hist3(XNM,edges)
% % [CountDataNM, Xedges,Yedges] = histcounts2(XNM(:,1),XNM(:,2),edges{1},edges{2});
% 
% % linear interpolarization of Distribution
% % Vqall = interp2(interp2(CountData,'cubic'),'cubic');
% % VqallNM = interp2(CountDataNM,'cubic');
% % VqallNM(find(VqallNM<0)) = 0;
% % pdfallNM = VqallNM/sum(sum(VqallNM));
% % figure
% % mesh(VqallNM)
% 
% % %kernel density
% 
% [x1,x2] = meshgrid(0:1:ceil(max(max(XNM))),0:1:ceil(max(max(XNM))));
% x1 = x1(:);
% x2 = x2(:);
% xi = [x1 x2];
% vqKS = ksdensity(XNM,xi,'Bandwidth',bw);
% pdfallNM = vqKS/sum(vqKS);
% 
% % %Gaussian regression
% % [x1, x2] = meshgrid(Xedges(1:end-1),Yedges(1:end-1));
% % x = [x1(:).'; x2(:).'];
% % % x1 = x1(:).';
% % % x2 = x2(:).';
% % xi = [x1; x2];
% % y = CountDataNM(:).';
% % for i = 1:size(x1,1)
% %     for j = 1:size(x2,2)
% %         xs = [x1(i,j); x2(i,j)];
% %         VqallNM(i,j) = gaussian_kern_reg2(xs,x,y,[1; 1]);
% %     end
% % end
% % 
% % VqallNM(find(VqallNM<0)) = 0;
% % 
% % pdfallNM = VqallNM/sum(sum((VqallNM))) ;
% % figure
% % mesh(x1,x2,ys)
% %     
% % Do the same thing for all frequencies
% 
% CountdatasNM = {};
% VqNM = {};
% pdfsNM = {};
% InfoNM = 0;
% 
% %Linear interpolation
% % for f = 1 :10
% %     CountdatasNM{f} = histcounts2(DR_trialNM1(:,f),DR_trialNM2(:,f),edges{1},edges{2});
% %     %     Vq{f} = interp2(interp2(Countdatas{f},'cubic'),'cubic');
% %     VqNM{f} = interp2(CountdatasNM{f},'cubic');
% %     VqNM{f}(find(VqNM{f}<0)) = 0;
% %     pdfsNM{f} = VqNM{f}/sum(sum(VqNM{f}));
% % end
% 
% % InfoNM = 0;
% % 
% % for f = 1:10
% %     for r1 = 1:length(pdfallNM)
% %         for r2 = 1:length(pdfallNM)
% %             if pdfsNM{f}(r1,r2) ~= 0 && pdfallNM(r1,r2) ~=0
% %                 InfoNM = InfoNM + pdfsNM{f}(r1,r2)*log2(pdfsNM{f}(r1,r2)/pdfallNM(r1,r2));
% %             end
% %         end
% %     end
% % end
% 
% % kernel density
% for f = 1:10
%     VqNM{f} = ksdensity([DR_trialNM1(:,f) DR_trialNM2(:,f)],xi,'Bandwidth',bw);
%     VqNM{f}(find(VqNM{f}<0)) = 0;
%     pdfsNM{f} = VqNM{f}/sum(sum(VqNM{f}));
%     for r1 = 1:length(pdfallNM)
%         if pdfsNM{f}(r1) ~= 0 && pdfallNM(r1) ~=0
%             InfoNM = InfoNM + pdfsNM{f}(r1)*log2(pdfsNM{f}(r1)/pdfallNM(r1));
%         end
%     end
% end
% 
% 
% % %Gaussian regression
% % 
% % for f = 1:10
% %     CountdatasNM{f} = histcounts2(DR_trialNM1(:,f),DR_trialNM2(:,f),edges{1},edges{2});
% %     y = CountdatasNM{f}(:).';
% %     for i = 1:size(x1,1)
% %         for j = 1:size(x2,2)
% %             xs = [x1(i,j); x2(i,j)];
% %             VqNM{f}(i,j) = gaussian_kern_reg2(xs,x,y,[1; 1]);
% %         end
% %     end
% %     
% %     VqNM{f}(find(VqNM{f}<0)) = 0;
% %     pdfsNM{f} = VqNM{f}/sum(sum(VqNM{f}));
% %     for r1 = 1:length(pdfallNM)
% %         for r2 = 1:length(pdfallNM)
% %             if pdfsNM{f}(r1,r2) ~= 0 && pdfallNM(r1,r2) ~=0
% %                 InfoNM = InfoNM + pdfsNM{f}(r1,r2)*log2(pdfsNM{f}(r1,r2)/pdfallNM(r1,r2));
% %             end
% %         end
% %     end
% %     
% % end
% 
% % Calculating Information content
% 
% 
% 
% InfoNM = InfoNM/10;

function Out = MI_ind(bw,d,Ntrial)
%% Calculating MI 07/11/2017

% Recalculating MI using a more simple method, directly comparable to joint
% MI

if d ==1
    load('SyncP_new.mat')
elseif d==2
    load('SyncN_new.mat')
elseif d==3
    load('SyncNM_new.mat')
elseif d==4
    load('SyncPMdlInfo.mat')
    output = UnitInfo.Info.Output;
elseif d==5
    load('SyncNMdlInfo.mat')
    output = UnitInfo.Info.Output;
else 
    load('SyncNMMdlInfo.mat')
    output = UnitInfo.Info.Output;
end
% NP = size(output.rates_stim{1, 1},1);
MI_P = [];

for trial = 1:Ntrial
    disp(trial)
    Sample_ind = randsample(1:size(output.rates_stim{1, 1},1),250);
    NP = 250;
    
    
    DR_trialP = [];
    for f = 1:10
        for n = 1:NP
            DR_trialP(n,f) = mean(output.rates_stim{f+1}(Sample_ind(n),:));
        end
    end
    
    DR_trialP(find(DR_trialP>100)) = 0;
    
    % 08/11/2017 attempt to linearly extrapolate a pdf.

    %Kernel density estimation
    SyncPall = reshape(DR_trialP,[],1);
    [vqKS,xi] = ksdensity(SyncPall,0:0.5: ceil(max(max(DR_trialP)))-1,'Bandwidth',bw);

    pdfPall = vqKS/sum(vqKS);
    
    InfoP = 0;
%     figure
    Hz_list = 8:4:48;
    pdfP = [];
    for f= 1:10

        % %     kernel density
        vq2N = ksdensity(DR_trialP(:,f),xi,'Bandwidth',bw);
        %below is common
        vq2N(find(vq2N <0)) = 0;
        vq2N= [vq2N(2:end) 0];
        pdfP(f,:) = vq2N/sum(vq2N);
        if mod(f,2) ==1
%         plot(pdfP(f,:),'DisplayName', [num2str(Hz_list(f)) 'Hz'])
        hold on
        end
        for r = 1:length(pdfPall)
            if pdfP(f,r) ~=0 && pdfPall(r) ~=0
                InfoP = InfoP + pdfP(f,r)*log2(pdfP(f,r)/pdfPall(r));
            end
        end
    end
    InfoP = InfoP/10;
    MI_P = [MI_P InfoP];
end

Out.mean = mean(MI_P);
Out.error = std(MI_P);




% for trial = 1:Ntrial
%     disp(trial)
%     Sample_ind = randsample(1:size(output.rates_stim{1, 1},1),250);
%     NP = 250;
%     
%     
%     DR_trialP = [];
%     for f = 1:10
%         for n = 1:NP
%             DR_trialP(n,f) = mean(output.rates_stim{f+1}(Sample_ind(n),:));
%         end
%     end
%     
%     DR_trialP(find(DR_trialP>100)) = 0;
%     
%     % 08/11/2017 attempt to linearly extrapolate a pdf.
%     edgesTotalP = 0:2: ceil(max(max(DR_trialP)));
%     [NtotalP,edgesTotalP] = histcounts(DR_trialP,edgesTotalP);
%     % for i = 1:length(NtotalP)
%     %     if NtotalP(i) == 0
%     %         NtotalP(i) = NaN;
%     %     end
%     % end
%     
%     % histogram(DR_trialP,edgesTotalP)
%     
%     % plot(NtotalP)
%     xq = 0:1: ceil(max(max(DR_trialP)));
%     % linear interpolation into infinite dataset
%     vq1 = interp1(edgesTotalP(2:end),NtotalP,xq(2:end),'pchip');
%     vq1(find(vq1 <0)) = 0;
%     vq1 = [vq1(2:end) 0];
%     pdfPall = vq1/sum(vq1);
%     % vq2 = interp1(xData, yData, xaxis,'spline');
%     % plot(xaxis,vq1)
%     % hold on
%     % plot(xaxis,vq2)
%     % vq2 = interp1(edgesTotalP(2:end),NtotalP,xq(2:end),'pchip');
%     % plot(vq2)
%     
%     
%     %Gaussian kernel based estimation.
%     edgesTotalPgauss = 0:1: ceil(max(max(DR_trialP)));
%     % [NtotalPgauss,edgesTotalPgauss] = histcounts(DR_trialP,edgesTotalPgauss);
%     xs = edgesTotalPgauss(1:end-1)+0.5;
%     h = 1.;
%     
%     
%     %Kernel density estimation
%     SyncPall = reshape(DR_trialP,[],1);
%     [vqKS,xi] = ksdensity(SyncPall,0:0.5: ceil(max(max(DR_trialP)))-1,'Bandwidth',bw);
%     % Offcial kernel density estimation, but I don't know how to use it yet.
%     
%     
%     %
%     % for i = 1:length(xs)
%     %     ys(i) = gaussian_kern_reg2(xs(i),xs,NtotalPgauss,h);
%     % end
%     % vqgauss = ys/sum(ys);
%     
%     %Compare
%     % figure
%     % plot(xs,vqgauss)
%     % hold on
%     % plot(xq(2:end),pdfPall)
%     % plot(xi,vqKS)
%     
%     
%     % pdfPall = vqgauss;
%     pdfPall = vqKS/sum(vqKS);
%     
%     InfoP = 0;
% %     figure
%     Hz_list = 8:4:48;
%     for f= 1:10
%         %linear interpolation
%         %     [NcountP,edgesTotalP] = histcounts(DR_trialP(:,f),edgesTotalP);
%         % %     NcountP(find(NcountP==0)) = NaN;
%         %     vq2N = interp1(edgesTotalP(2:end),NcountP,xq(2:end),'pchip');
%         
%         % %     gaussian kernel
%         %     [NcountPgauss,edgesTotalPgauss] = histcounts(DR_trialP(:,f),edgesTotalPgauss);
%         %     for i = 1:length(xs)
%         %         vq2N(i) = gaussian_kern_reg2(xs(i),xs,NcountPgauss,h);
%         %     end
%         %
%         % %     kernel density
%         vq2N = ksdensity(DR_trialP(:,f),xi,'Bandwidth',bw);
%         %below is common
%         vq2N(find(vq2N <0)) = 0;
%         vq2N= [vq2N(2:end) 0];
%         pdfP(f,:) = vq2N/sum(vq2N);
%         if mod(f,2) ==1
% %         plot(pdfP(f,:),'DisplayName', [num2str(Hz_list(f)) 'Hz'])
%         hold on
%         end
%         for r = 1:length(pdfPall)
%             if pdfP(f,r) ~=0 && pdfPall(r) ~=0
%                 InfoP = InfoP + pdfP(f,r)*log2(pdfP(f,r)/pdfPall(r));
%             end
%         end
%     end
%     InfoP = InfoP/10;
%     MI_P = [MI_P InfoP];
% end
% 
% Out.mean = mean(MI_P);
% Out.error = std(MI_P);


function ISI_analy()
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