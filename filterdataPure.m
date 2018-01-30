function filterdata()

clear all

% will be needed when taking all data
%load('infoset.mat');
%animals = unique(cellfun(@char,{output.animal},'unif',0));

%% file info
animal_list = {'m36n','m2p','m41o','m32q'};

% createFileInfo('m2p') %animal namecode
%
animal = animal_list{2};

UnitInfo = createFileInfo(animal);

% load([animal '_List3']);
%
% ICI_list1 = [2 2.5 3 5 7.5 10 12.5 15 20 25 30 35 40 45 50 55 60 65 70 75]; %ms, ICI
% ICI_list2 = [250 125 83.3333 62.5 50 41.6667 35.7143 31.25 27.7778 25 22.7273 20.8333];
% Attn_List = [20 30 40 50 60 70 80 90 100];
% Attn_List = [10 20 30 40 50 60 70 80 90 100];
% %
global directory
directory = ['U:\Neural and Behavioural Data' '\marmoset\' animal];

% directory = ['C:\Users\John\Dropbox\Bendorlab\Code\Marmoset\' animal];

%% Spike extraction
for n = 1:length(UnitInfo.List)
    %     disp(n)
    
    spiketable = load([directory filesep UnitInfo.List{1,n}]);
    if ~isempty(spiketable)
        PREstim_duration = UnitInfo.Info(n).Pre_stim_Duration/1000; %seconds
        stim_duration = UnitInfo.Info(n).stim_duration;
        TrialLength = UnitInfo.Info(n).Pre_stim_Duration + UnitInfo.Info(n).Post_stim_Duration +stim_duration;
        StimNb = unique(spiketable(:,1)).';
        if max(StimNb) ==9
            UnitInfo.Info(n).Stimuli_Nb = 1;
            Attn_List = [20 30 40 50 60 70 80 90 100];
        elseif max(StimNb) == 10
            UnitInfo.Info(n).Stimuli_Nb = 2;
            Attn_List = [10 20 30 40 50 60 70 80 90 100];
        end
    end
end


for f = 1:9
    output.rates_stim{f} = [];
    output.rates_pre{f} = [];
    output.rates_post{f} = [];
    output.isi_total{f} = [];
end


indxList = find([UnitInfo.Info.Stimuli_Nb] == 2 & [UnitInfo.Info.stim_duration] == 200);

Attn_List = [20 30 40 50 60 70 80 90 100];
Attn_List = [10 20 30 40 50 60 70 80 90 100];
for f = 1:length(Attn_List)
    output.rates_stim{f} = [];
    output.rates_pre{f} = [];
    output.rates_post{f} = [];
    output.isi_total{f} = [];
end


%         & [UnitInfo.Info.Stimuli_Nb] ==StimType);
% figure
for n = indxList
    disp(UnitInfo.List{1,n})
    spiketable = load([directory filesep UnitInfo.List{1,n}]);
    
    if ~isempty(spiketable)
        
        PREstim_duration = UnitInfo.Info(n).Pre_stim_Duration/1000; %seconds
        stim_duration = UnitInfo.Info(n).stim_duration;
        TrialLength = UnitInfo.Info(n).Pre_stim_Duration + UnitInfo.Info(n).Post_stim_Duration +stim_duration;
        StimNb = unique(spiketable(:,1)).';
        
        
        raster.stim=[];  raster.rep=[];  raster.spikes=[];
        spikes_pooled_for_VS = [];
        all_mean_rate_spont = [];
        all_std_rate_spont = [];
        all_mean_rate_stim = [];
        all_std_rate_stim = [];
        %     all_rate_total = [];
        Fanofactor = [];
        emptycount = 0;
        
        
        
        for f = StimNb
            nbrep = spiketable(find(spiketable(:,1)==f),2);
            nreps = max(nbrep.');
            indstim = find(spiketable(:,1)==f);
            spikes_pooled = [];
            rate_total = [];
            if ~isempty(nbrep)
                for r = unique(nbrep.')
                    spikes1 = []; %channel 1
                    spikes2 = []; %channel 2
                    for i = indstim.'
                        if spiketable(i,2) == r && spiketable(i,3)== 1 && spiketable(i,4) > 0 % channel 1
                            spikes1 = [spikes1 (spiketable(i,4)*1e-6 - PREstim_duration)]; % spikes rounded up to 0.1ms cad 100microseconds. converted to seconds
                        end
                        if spiketable(i,2) == r && spiketable(i,3)== 2 && spiketable(i,4) > 0 % channel 2
                            spikes2 = [spikes2 (spiketable(i,4)*1e-6 - PREstim_duration)];
                        end
                    end
                    if length(spikes1) > length(spikes2)
                        spikes =  spikes1;
                    else
                        spikes = spikes2;
                    end
                    if isempty(spikes(find(spikes>PREstim_duration & spikes< PREstim_duration + 0.2)))
                        emptycount = emptycount+1;
                    end
                    
                    spikes_pooled = [spikes_pooled spikes];
                    
                    raster.stim=[raster.stim f*ones(size(spikes))];
                    raster.rep=[raster.rep r*ones(size(spikes))];
                    raster.spikes=[raster.spikes spikes];
                    
                    % rate
                    rate = zeros(1,TrialLength);
                    spikes4rate = spikes(find(spikes<PREstim_duration+0.2)) + PREstim_duration;
                    
                    
                    for st = spikes4rate
                        if ceil(st*1e3) <= length(rate)
                            rate(1,ceil(st*1e3)) = rate(1,ceil(st*1e3))+1;
                        end
                        
                    end
                    rate_total = [rate_total ; rate*1000];
                    
                    
                end
                %calculating Fano factor
                SpikeCount = sum(rate_total(:,PREstim_duration*1000+1:PREstim_duration*1000+stim_duration+100),2)/1000.;
                Fanofactor = [Fanofactor std(SpikeCount)^2/mean(SpikeCount)];
                
                %average rate + gaussian kernel
                %             rate_av = mean(rate_total,1);
                %             xs = 1:TrialLength;
                %             h = 15; %kernal bandwidth. determines the shape of the function
                %             for i = 1:TrialLength
                %                 ys(i)=gaussian_kern_reg(xs(i),xs,rate_av,h);
                %             end
                %
                output.rates_stim{f} = [output.rates_stim{f}; ...
                    rate_total(:,PREstim_duration*1000+1:PREstim_duration*1000+stim_duration+100)];
                output.rates_pre{f} = [output.rates_pre{f}; ...
                    rate_total(:,PREstim_duration*1000-99:PREstim_duration*1000)];
                
                %
                %             all_rate{f} = ys;
                %             %             all_rate_total = [all_rate_total ; rate_total];
                %             all_mean_rate_spont(f) = mean2(rate_total(:,1:PREstim_duration*1000));
                %             all_std_rate_spont(f) = std2(rate_total(:,1:PREstim_duration*1000))/sqrt(PREstim_duration*1000*nreps);
                %             all_mean_rate_stim(f) = mean2(rate_total(:,PREstim_duration*1000+1:PREstim_duration*1000+stim_duration+100));% -all_mean_rate_spont(f);
                %             all_std_rate_stim(f) = std2(rate_total(:,PREstim_duration*1000+1:PREstim_duration*1000+stim_duration+100))/sqrt((stim_duration+100)*nreps);
                % %             spikes_pooled_for_VS=spikes_pooled(find(spikes_pooled>0.05 & spikes_pooled<=(stim_duration+0.05)));
                
                
                %%% checking individual neurons
%                 output.mean_rate_stim{f} = mean(output.rates_stim{f},1);
%                 output.mean_rate_pre{f} = mean(output.rates_pre{f},1);
%                 xs = 1:400;
%                 h = 10;
%                 for i = 1:400
%                     ys(i) = gaussian_kern_reg(xs(i),xs,[output.mean_rate_pre{f} output.mean_rate_stim{f}],h);
%                 end
%                 smooth_rate{f} = ys;
            end
            
        end
        
        
        %         figure
%         hold off
%         cmap = colormap(jet(length(Attn_List)));
%         for f = 1: length(Attn_List)
%             plot(smooth_rate{f},'color',cmap(f,:),'LineWidth',1.7,'DisplayName', ...
%                 [num2str(Attn_List(f)) 'dB'] );
%             %     lgd{f} = [num2str(ceil(1000/ICI_list(f))) 'Hz'];
%             hold on
%         end
%         legend('show')
%         pause(0.3)
%         
%         for f = 1:length(Attn_List)
%             output.rates_stim{f} = [];
%             output.rates_pre{f} = [];
%             output.rates_post{f} = [];
%             output.isi_total{f} = [];
%         end
         %%% /checking individual neurons
        test = 1;
    end
end
%% mean activity and stuff

output.Fanofactor = [];
% vector = zeros(1,length(ICI_list));

% output.meanVS = mean(output.VS,1);
% output.errorVS = std(output.VS,1)/sqrt(nn);


output.meanDR = [];
output.errorDR = [];

output.var_isi = [];
for f = 1:length(Attn_List)
    output.mean_rate_stim{f} = mean(output.rates_stim{f},1);
    output.std_rate_stim{f} = std(output.rates_stim{f},0,1);
    output.mean_rate_pre{f} = mean(output.rates_pre{f},1);
    output.std_rate_pre{f} = std(output.rates_pre{f},0,1);
    %     output.mean_rate_post{f} = mean(output.rates_post{f},1);
    SpikeCount = sum(output.rates_stim{f},2)/1000.;
    output.Fanofactor = [output.Fanofactor std(SpikeCount)^2/mean(SpikeCount)];
    output.var_isi = [output.var_isi std(output.isi_total{f})];
    
    %Gaussian smoothing
    xs = 1:400;
    h = 10;
    for i = 1:400
        ys(i) = gaussian_kern_reg(xs(i),xs,[output.mean_rate_pre{f} output.mean_rate_stim{f}],h);
    end
    smooth_rate{f} = ys;
    output.meanDR = [ output.meanDR mean2(output.rates_stim{f})];
    output.errorDR = [output.errorDR std2(output.rates_stim{f})/sqrt(size(output.rates_stim{1},1)*size(output.rates_stim{1},2))];
end

figure
cmap = colormap(jet(length(Attn_List)));
for f = 1: length(Attn_List)
    plot(smooth_rate{f},'color',cmap(f,:),'LineWidth',1.7,'DisplayName', ...
        [num2str(Attn_List(f)) 'dB'] );
    %     lgd{f} = [num2str(ceil(1000/ICI_list(f))) 'Hz'];
    hold on
end
legend('show')



test = 1
%% Create fileinfo
function UnitInfo = createFileInfo(animal)

% animal = 'm2p';

global directory
directory = ['U:\Neural and Behavioural Data' '\marmoset\' animal];
% directory = ['C:\Users\John\Dropbox\Bendorlab\Code\Marmoset\' animal];
filelist = dir(directory);
filelist = filelist(3:end,:);
filelist = extractfield(filelist,'name');
UnitInfo = {};
UnitInfo.List = []; i = 1;
UnitInfo.Info = [];
UnitInfo.Info = struct('Channel_Nb',{},'Stimuli_Nb',{},'Pre_stim_Duration',{},'Post_stim_Duration',{},'Sync',{},'Significant_rate',{},'Positive',{},'Rho',{});% for click trains
% UnitInfo.Info = struct('Channel_Nb',{},'tone_Hz',{},'Attn_dB',{},'Pre_stim_Duration',{},'Post_stim_Duration',{},'Stim_duration',{});

% filelist{1} = 'M2p0010.dat';
for k = 1 :length(filelist)
    %     fid = fopen(['U:\Neural and Behavioural Data' filesep 'marmoset' filesep animal filesep filelist{k}]);
    fid = fopen([directory filesep filelist{k}]);
    frewind(fid)
    fileinfo = textscan(fid,'%s %s',30,'delimiter',{' = ','% '}...
        ,'MultipleDelimsAsOne',1,'headerlines', 6);
    fclose(fid);
%     A = fileread(['U:\Neural and Behavioural Data' filesep 'marmoset' filesep animal filesep filelist{k}]);
    A = fileread([directory filesep filelist{k}]);
    
    A1 = strfind(A,'Pre-stimulus record time (ms) = ');
    A2 = strfind(A,'Post-stimulus record time (ms) = ');
    %     A3 = strfind(A,'Number of Stimuli = ');
    A3 = strfind(A,'Number of Stimuli = ');
    l1 = length('Pre-stimulus record time (ms) = ');
    l2 = length('Post-stimulus record time (ms) = ');
    l3 = length('Number of Stimuli = ');
    
    Prestim = str2num(A(A1+l1:A1+l1+2));
    Poststim = str2num(A(A2+l2:A2+l2+2));
    NbOstimuli = str2num(A(A3+l3:A3+l3+1));
    %     if strcmp(fileinfo{1,2}{1,1}, '1')==1
    %         UnitInfo.List{1,i} = filelist{k};
    % %         UnitInfo.Info(i).Channel_Nb = str2num(fileinfo{1,2}{10,1});
    %         UnitInfo.Info(i).Stimuli_Nb = NbOstimuli;
    %         UnitInfo.Info(i).Pre_stim_Duration = Prestim;
    %         UnitInfo.Info(i).Post_stim_Duration = Poststim;
    %         i = i+1
    %     end
    % end
    
    if strcmp(fileinfo{1,2}{1,1}, '2')==1 %This is for CLICK TRAINS
        if k ~= 7174
            UnitInfo.List{1,i} = filelist{k};
            UnitInfo.Info(i).Channel_Nb = str2num(fileinfo{1,2}{10,1});
            %         UnitInfo.Info(i).Stimuli_Nb = str2num(fileinfo{1,2}{17,1});
            UnitInfo.Info(i).Stimuli_Nb = 0;
            UnitInfo.Info(i).Pre_stim_Duration = Prestim;
            UnitInfo.Info(i).Post_stim_Duration = Poststim;
            kk = strfind(A,'duration(ms)');
            l3 = length('duration(ms)');
            stim =  str2num(A(kk(1,1)+l3:kk(1,1)+l3+2));
            UnitInfo.Info(i).stim_duration = stim ;
            i = i+1;
            disp(filelist{k})
        end
    end
    
    %         if strcmp(fileinfo{1,2}{16,1}, '12')==1 || strcmp(fileinfo{1,2}{16,1}, '20')==1
    %             %         pause
    %             UnitInfo.List{1,i} = filelist{k};
    %             UnitInfo.Info(i).Channel_Nb = str2num(fileinfo{1,2}{10,1});
    %             UnitInfo.Info(i).Stimuli_Nb = str2num(fileinfo{1,2}{16,1});
    %             UnitInfo.Info(i).Pre_stim_Duration = Prestim;
    %             UnitInfo.Info(i).Post_stim_Duration = Poststim;
    %             i = i+1
    %         elseif strcmp(fileinfo{1,2}{17,1}, '12')==1 || strcmp(fileinfo{1,2}{17,1}, '20')==1
    %             UnitInfo.List{1,i} = filelist{k};
    %             UnitInfo.Info(i).Channel_Nb = str2num(fileinfo{1,2}{10,1});
    %             UnitInfo.Info(i).Stimuli_Nb = str2num(fileinfo{1,2}{17,1});
    %             UnitInfo.Info(i).Pre_stim_Duration = Prestim;
    %             UnitInfo.Info(i).Post_stim_Duration = Poststim;
    %             i = i+1
    %         end
    %     end
end

%% plots

function plotdata(ICI_list, vector, raster, nreps, stim_duration,total_rate)

%     if vector(end)>0;
subplot(2,2,1)
plot(ICI_list,vector)
axis([0 80 -0.1 1.1])
%     pause()


subplot(2,2,3)
cmap = colormap(hsv(length(ICI_list)));
x_axis = linspace(-0.5, 1.0, 1500);
for n = 1:length(ICI_list)
    plot(x_axis,total_rate{n},'color',cmap(n,:),'DisplayName',num2str(ICI_list(n)))
    hold on
    
end
legend('show')

axis([-0.2,0.8,0,100]);


subplot(2,2,[2,4])
xlabel('time (s)')
ylabel('IPI (ms)')

ytext = [];
ytext_position = [];
for s = 1:length(ICI_list)
    ytext_position = [ytext_position;s*nreps-nreps/2];
    ytext = char(ytext,num2str(ICI_list(s)));
end
ytext(1,:) = [];
area([0 stim_duration stim_duration 0],[length(ICI_list)*nreps+1 length(ICI_list)*nreps+1 0 0],'LineStyle','none','FaceColor',[.85 .85 1]);
hold on
plot(raster.spikes,nreps*(raster.stim-1)+raster.rep,'k.','MarkerSize',9);
set(gca,'yTick',ytext_position);
set(gca,'yTickLabel',ytext);
%     end



%% GUI










