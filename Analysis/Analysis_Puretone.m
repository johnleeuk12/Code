%% Puretone  06/06/2018

% real data, median spike time comparisons. 
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

% model with SFA