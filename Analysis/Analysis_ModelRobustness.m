%% Model robustness 04/06/2018

clear all

load('modeldataSN_Jitter.mat') %load file, SP, SN, Noise, Jitter
% noise = 0:2:16;
noise = 0:10; % this is for jitter 
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
        collect1{i} = [collect1{i} UnitInfo.Info((i-1)*10+trial).Rho]; %pooling data
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
% axis([-2 18 -1.2 1.2]) %for noise
axis([-2 12 -1.2 1.2]) % for jitter

figure
errorbar(noise,mean_spont,error_spont);

figure
plot(UnitInfo.Info(30).Output.mean_discharge_rate.mean) %test for individual neuron. 

%% Vector strength robustness 04/06/2018


clear all
load('modeldataSN_Noise.mat')

noise = 0:2:16;
% noise = 0:10;
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
set(gca,'YTickLabel',{'0' '2' '4' '6' '8' '10' '12' '14' '16'});
% set(gca,'YTickLabel',{'0' '1' '2' '3' '4' '5' '6' '7' '8' '9' '10'});
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
