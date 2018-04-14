%% Figure 1
% load('datafig1.mat')
% DR = zeros(1,3);
% DR(1,3) = outputmdl.mean_discharge_rate.mean(1)-outputmdl.mean_discharge_rate.mean(end);
% DR(1,1) = outputSyncP.meanDR(end)-outputSyncP.meanDR(2);
% DR(1,2) = outputSyncN.meanDR(end)-outputSyncN.meanDR(2);
% 
% 
% 
% 
% 
% DRerror = zeros(1,3);
% DRerror(1,3) = mean(outputmdl.mean_discharge_rate.error);
% DRerror(1,2) = mean(outputSyncN.errorDR);
% DRerror(1,1) = mean(outputSyncP.errorDR);
% 
% cmapp = [[0.1 0.7 0.1]; [0.9 0.6 0.1]; [0 0 0]   ];%    [0.26 0.5 0.9]   ]; %; [0.9 0.3 0.26]; ]; 
% figure
% for ind = 1:3
%     bar(ind, DR(1,ind),0.5,'FaceColor', 'none', 'EdgeColor',cmapp(ind,:),'LineWidth',3)
%     hold on
%     errorbar(ind, DR(1,ind),DRerror(1,ind),'.','Capsize', 50,'Color',cmapp(ind,:),'LineWidth',3)
%     set(gca, 'FontSize', 16)
% end
% 
% VS = zeros(6,3);
% VS(:,1) =  fliplr(outputSyncP.meanVS(1,2:2:end));
% VS(:,2) =  fliplr(outputSyncN.meanVS(1,2:2:end));
% VS(:,3) = outputmdl.VS.';
% errorVS = zeros(6,3);
% errorVS(:,1) = fliplr(outputSyncP.errorVS(1,2:2:end)).';
% errorVS(:,2) = fliplr(outputSyncN.errorVS(1,2:2:end)).';
% % errorVS(:,3) = fliplr(outputSyncN.errorVS(1,2:2:end)).';
% 
% for n = 1:length(UnitInfo.Info)
%     for f = 1:length(ICI_list)
%         totalVS(n,f) = UnitInfo.Info(n).Output.VS(f);
%     end
% end
% 
% meanVS = mean(totalVS,1);
% stdVS = std(totalVS,1);
% 
% 
% figure
% for ind = 1:2
%     hold on
%     errorbar(Hz_list,VS(:,ind).',errorVS(:,ind).','Capsize', 15,'LineWidth',3,'Color', cmapp(ind,:))
%     
% end
% 
% % plot(Hz_list,VS(:,3).','LineWidth',3,'Color', cmapp(3,:))
% errorbar(fliplr(Hz_list),meanVS(1,1:2:end),stdVS(1,1:2:end),'Capsize', 15,'LineWidth',3,'Color', cmapp(3,:))
% set(gca, 'FontSize', 16)
% 




%% Figure 5 on the previous model, and model parameters 
% % 
% load('model_datasync2.mat')
tau_PE = unique(UnitInfo.List(:,3)).'; %x
tau_PI = unique(UnitInfo.List(:,4)).';%y
% E =  unique(UnitInfo.List(:,5)).';
% I =  unique(UnitInfo.List(:,6)).';
f_DE = 0.5:0.1:1. ; %a
f_DI = 0.5:0.1:1. ; %b 

Z_positive1 = zeros(length(tau_PE),length(tau_PI));
Z_negative1 = zeros(length(tau_PE),length(tau_PI));
Z = zeros(length(f_DE),length(f_DI));
Z_onset = zeros(length(f_DE),length(f_DI));


for x = 1:length(tau_PE)
    for y = 1:length(tau_PI)
%         tic
        idx = find(UnitInfo.List(:,3) == tau_PE(x) & UnitInfo.List(:,4) == tau_PI(y)).';
        for z = idx 
            if UnitInfo.Info(z).Positive ==1 %1 is positive
                Z_positive1(x,y) = Z_positive1(x,y) - UnitInfo.Info(z).Rho;
            elseif UnitInfo.Info(z).Positive == -1
                Z_negative1(x,y) = Z_negative1(x,y) - UnitInfo.Info(z).Rho;
            end
        end
%         toc
    end
end



for a = 1:length(f_DE)
    for b = 1:length(f_DI)
        z = find(UnitInfo.List(:,1) == f_DE(a) & UnitInfo.List(:,2) == f_DI(b) &...
            UnitInfo.List(:,3) == tau_PE(7) & UnitInfo.List(:,4) == tau_PI(6));
        Z(a,b) = -UnitInfo.Info(z).Rho;
        Z_onset(a,b) = UnitInfo.Info(z).Output.spikes_per_click{1}.mean(1);
    end
end


Z_positive1 = Z_positive1/36; 
Z_negative1 = Z_negative1/36;

A_DE = 1-f_DE;
A_DI = 1-f_DI;

figure
imagesc(tau_PI,tau_PE,Z_positive1,[-1, 1])
% surf(tau_PI,tau_PE,Z_positive1)
ylabel('\tau_pE')
xlabel('\tau_pI')
zlabel('Rho')
title('Sync+')
colormap(flipud(jet))
set(gca, 'FontSize', 16)
colorbar
set(colorbar, 'ylim', [0 1])

figure
imagesc(tau_PI,tau_PE,Z_negative1,[-1, 1])
% surf(tau_PI,tau_PE,Z_negative1)
ylabel('\tau_pE')
xlabel('\tau_pI')
zlabel('Rho')
title('Sync-')
colormap(flipud(jet))
set(gca, 'FontSize', 16)
colorbar
set(colorbar, 'ylim', [-1 0])

% figure
% imagesc(A_DI,A_DE,Z)
% % surf(f_DI,f_DE,Z)
% ylabel('A_DE')
% xlabel('A_DI')
% zlabel('Rho')
% colormap(flipud(jet))
% title('with optimal \tau_pE and \tau_pI')
%     set(gca, 'FontSize', 16)
% 
% figure
% imagesc(A_DI,A_DE,Z_onset)
% % surf(f_DI,f_DE,Z_onset)
% ylabel('A_DE')
% xlabel('A_DI')
% zlabel('rate (spikes/s')
% title('onset response')
% % colormap hot 
% % size(Z_positive)
% set(gca, 'FontSize', 16)
% test =1;
% 
% 
% % %%% from here, EI 
% idx = find([UnitInfo.Info.Positive] == 1 & [UnitInfo.Info.Significant_rate] == 1);
% idx2 = find([UnitInfo.Info.Positive] == 0 & [UnitInfo.Info.Significant_rate] == 1);
% idx3 = find([UnitInfo.Info.Positive] == -1 & [UnitInfo.Info.Significant_rate] == 1);
% for i = 1:length(idx)
%     x1(i) = UnitInfo.List(idx(i),5);
%     y1(i) = UnitInfo.List(idx(i),6)/x1(i);
%     z1(i) = UnitInfo.Info(idx(i)).Rho;
% end
% 
% 
% 
% for i = 1:length(idx2)
%     x2(i) = UnitInfo.List(idx2(i),5);
%     y2(i) = UnitInfo.List(idx2(i),6)/x2(i);
%     z2(i) = UnitInfo.Info(idx2(i)).Rho;
% end
% 
% for i = 1:length(idx3)
%     x3(i) = UnitInfo.List(idx3(i),5);
%     y3(i) = UnitInfo.List(idx3(i),6)/x3(i);
%     z3(i) = UnitInfo.Info(idx3(i)).Rho;
% end
%  
% figure 
% cmapp = [[0.1 0.7 0.1]; [0.9 0.6 0.1]; [0 0 0]   ];%    [0.26 0.5 0.9]   ]; %; [0.9 0.3 0.26]; ]; 
% 
% scatter(x1,y1,'DisplayName','Sync+','LineWidth',1.5,'MarkerEdgeColor',cmapp(1,:))
% hold on 
% % pause
% scatter(x2,y2,'x','DisplayName','Non Significant','LineWidth',1.5,'MarkerEdgeColor',cmapp(3,:))
% % pause
% scatter(x3,y3,'DisplayName','Sync-','LineWidth',1.5,'MarkerEdgeColor',cmapp(2,:))
% axis([0, 7, 0, 7])
% set(gca, 'FontSize', 16)
% grid on
% legend show
% test=1 ;
% scatter3(x1,y1,z1)
% hold on 
% scatter3(x2,y2,z2)






%% Figure model

% cmapp = [[0.1 0.7 0.1]; [0.9 0.6 0.1]; [0 0 0]   ];
% n = 220;
% Hz_list = [];
% for i = 1:length(ICI_list)
% Hz_list = [Hz_list round(1000/ICI_list(i))];
% end
% figure
% % figure('units','normalized','outerposition',[0 0 1 1])
% cmap = colormap(jet(length(ICI_list)));
% subplot(2,2,1)
% 
% hold off
% for f = 1:length(ICI_list)
% plot(output.rate{n,f},'linewidth', 2.0, 'DisplayName', ...
%         [num2str(ceil(1000/ICI_list(f))) 'Hz'] );
% lgd{f} = [num2str(ceil(1000/ICI_list(f))) 'Hz'];
% hold on
% % % % x = 510:ICI_list(f):1110;
% % % % x = [x; ones(1,length(x))*(f-1)];
% % % % scatter(x(1,:),x(2,:),'k','filled')
% end
% axis([300,1300,0 80])
% set(gca, 'FontSize', 16)
% legend('show')
% 
% % subplot(2,2,2)
% % uitable('Data',d,'ColumnFormat', columnformat)
% subplot(2,2,3)
% plot(Hz_list,output.VS{n})
% hold on
% subplot(2,2,4)
% shadedErrorBar(Hz_list,output.DRmean{n},output.DRstd{n},{'--','Color',cmapp(2,:)})
% test = 1;
% hold on


%% Figure 4 neurons

% animal_list = {'m36n','m2p','m41o','m32q'};
% 
% % createFileInfo('m2p') %animal namecode
% %
% 
% TotalRho = [];
% VsNb = [];
% RayMean = [];
% categ = []; %1 Sync+ 2 Sync- 3 Nsync+ 4 Nsync-
% for ani = 1:4
%     animal = animal_list{ani};
%     load([animal '_List3']);
% 
%     indxList = find( [UnitInfo.Info.Stimuli_Nb] == 12); % & [UnitInfo.Info.Positive] == Positive ...        & [UnitInfo.Info.Stimuli_Nb] ==StimType);
%     % Criteria, nb of VS>0.1
%     for n = indxList
%         count = 0;
%         TotalRho = [TotalRho -UnitInfo.Info(n).Rho];
%         %         for f = 1:12
%         %             if UnitInfo.Info(n).output.vector(f)>0.1
%         %                 count = count+1;
%         %             end
%         %
%         %         end
%         %         VsNb = [VsNb count];
%         A = sort(UnitInfo.Info(n).output.rayleigh,'descend');
%         RayMean = [RayMean median(UnitInfo.Info(n).output.rayleigh(3:5),2)];
%         if UnitInfo.Info(n).Positive == -1 && UnitInfo.Info(n).Sync == 1
%             categ = [categ 1];
%         elseif UnitInfo.Info(n).Positive == 1 && UnitInfo.Info(n).Sync == 1
%             categ = [categ 2];
%         elseif UnitInfo.Info(n).Positive == -1 && UnitInfo.Info(n).Sync == 0
%             categ = [categ 3];
%         elseif UnitInfo.Info(n).Positive == 1 && UnitInfo.Info(n).Sync == 0
%             categ = [categ 4];
%         else
%             categ = [categ 5];
%         end
%     end
% end
% 
% figure
% test = 1;
% 
% scatter(TotalRho,RayMean,[],categ,'LineWidth',1.5)
% 
% % test = 1;
% 
% %% Rasterplot for Nsync
% i=0
% while i<11
% NN = randsample(NsyncIdx,1);
% if NN<100 || NN>900
%     NN = randsample(NsyncIdx,1);
% else
%     figure
%     idx = find([output.raster.rep{NN}] <11);
%     plot(output.raster.spikes{NN}(idx),10*(output.raster.stim{NN}(idx)-1)+output.raster.rep{NN}(idx),'k.','MarkerSize',9);
%     i = i+1;
% end
% end
% 
% %% Stats with network\
% 
% NsyncIdx([1:50]) = [];
% 
% nonSig.vector = [];
% nonSig.Rho =[];
% Sig.vector = [];
% Sig.Rho = [];
% test = []
%     
% NsyncIdx2 = 1:1450;
% NsyncIdx2(NsyncIdx) = [];
% NsyncIdx2(1:50) = [];
% 
% for n = NsyncIdx2
%     nonSig.vector = [nonSig.vector; output.VS{n}];
%     nonSig.Rho = [nonSig.Rho; output.Analysis(n).Rho];
% end
% 
% for nn = NsyncIdx
%     Sig.vector = [Sig.vector; output.VS{nn}];
%     Sig.Rho = [Sig.Rho; output.Analysis(nn).Rho];
% 
% end
% figure
% 
% plot(median(Sig.vector,1))
% hold on 
% plot(median(nonSig.vector,1))
% 
% figure
% 
% boxplot([Sig.Rho;nonSig.Rho],[ones(1,length(NsyncIdx)) zeros(1,length(NsyncIdx2))])
% set(gca, 'FontSize', 16)


























