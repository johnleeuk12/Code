load('dataforstat.mat')

tabels = {};

mdl = {};
for f = 1:5
tabels{f} = table(output.spikes_per_click{f}.mean.', output.net_positivePclick{f}.mean.',...
    output.net_negativePclick{f}.mean.','VariableNames',{'spikes','net_Positive','net_negative',});
% onset_spikes = [onset_spikes tabels{f}.spikes(1)];
% onset_excPlus = [onset_excPlus tabels{f}.net_Positive(1)];
% mdl{f} = stepwiselm(tabels{f},'spikes~net_Positive + net_negative');
%    
end
% 
% % lm = fitlm(tabel,'spikes~time_period+gE+time_period*gE');
% md3 = stepwiselm(tabel,'spikes~net_Positive + net_negative')
% % md2 = stepwiselm(tabel,'time_period~gE+gI+gI*gE');




clear all
load('modeldata.mat')

onset_spikes = [];
onset_excPlus = [];
spikes = [];
NetPosi = [];
NetNega = [];
Norma = [];


for n = 1 %:length(UnitInfo.List)
    %     onset_spikes = [onset_spikes UnitInfo.Info(n).Output.spikes_per_click{1}.mean(1)];
    %     onset_excPlus = [onset_excPlus UnitInfo.Info(n).Output.net_positivePclick{1}.mean(1) ];
    
    for f = 1:5
        for sp = 1:length(UnitInfo.Info(n).Output.spikes_per_click{f}.mean)
            spikes = [spikes; UnitInfo.Info(n).Output.spikes_per_click{f}.mean(sp)];
            NetPosi = [NetPosi; UnitInfo.Info(n).Output.net_positivePclick{f}.mean(sp)];
            NetNega = [NetNega; UnitInfo.Info(n).Output.net_negativePclick{f}.mean(sp)];
        end
    end
end
for n = 1:length(spikes)
    Norma(n,1) = spikes(n,1)-(NetPosi(n,1)*5.05*1e7-15.76);
end


tabel2 = table(NetNega,Norma,'VariableNames',{'net_negative','normalised_spikes'});
tabel_onset = table(onset_spikes.',onset_excPlus.','VariableNames',{'spikes','net_Positive'});
lm = fitlm(tabel_onset,'spikes ~ net_Positive')







tabel2 = {}
norma = [];
for n = 1:25
    norma(n,1) = tabels{1}.spikes(n) - (tabels{1}.net_Positive(n)*5.05*1e7 -15.76);
end
tabel2 = table(tabels{1}.net_negative, norma,'VariableNames',{'net_negative','spikes_norm'});
net_nega = tabels{1}.net_negative;
