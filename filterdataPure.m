function filterdata()

clear all

% will be needed when taking all data
%load('infoset.mat');
%animals = unique(cellfun(@char,{output.animal},'unif',0));

%% file info
animal_list = {'m36n','m2p','m41o','m32q'};

% createFileInfo('m2p') %animal namecode
%

nn = 1;
output.rates_stim{1} = [];
output.rates_pre{1} = [];
output.rates_post{1} = [];
output.isi_total{1} = [];
figure
for i = 2:4
    animal = animal_list{i};
    load([animal '_List3']);
    
    SyncInfo = UnitInfo;
    Sync = 1; %or 0 for Nsync
    Positive = 1; % or 1 for Negative
    StimType = 12; % or 20
    
    % for i = 2:4
    %     animal = animal_list{i};
    % UnitInfo = createFileInfo(animal);
    %
    % filename = [animal '_ListPure.mat'];
    % save(filename,'UnitInfo');
    % end
    % %
    load([animal '_ListPure']);
    
    % %
    global directory
    % directory = ['U:\Neural and Behavioural Data' '\marmoset\' animal];
    directory = ['C:\Users\John.Lee\OneDrive\Bendorlab\Marmoset\' animal];
    
    
    
    
    % Adding attenuation info on SyncInfo.
    for k = 1 :length(SyncInfo.List)
        A = fileread([directory filesep SyncInfo.List{k}]);
        A5 = strfind(A,'Attenuation(s) Used (dB) = ');
        l5 = length('Attenuation(s) Used (dB) = ');
        Attenuation = str2num(A(A5+l5:A5+l5+3));
        if isempty(Attenuation) == 1
            Attenuation = str2num(A(A5+l5:A5+l5+2));
        end
        SyncInfo.Info(k).attn =Attenuation;
    end
    
    %defining indxlist
    indxList = find([SyncInfo.Info.Stimuli_Nb] ==StimType &...
        [SyncInfo.Info.Positive] == Positive &...
        [SyncInfo.Info.Sync] == Sync &  ...
        [SyncInfo.Info.Pre_stim_Duration] == 500 &...
        [SyncInfo.Info.Post_stim_Duration] == 500 & [SyncInfo.Info.Significant_rate] ==1 );
    
    
    %initializing output vectors
    
    
    
    for n = indxList
        Channel = find( [UnitInfo.Info.Channel_Nb] == SyncInfo.Info(n).Channel_Nb); %Corresponding index for puretones matching Sync response Channel number.
        for k = Channel
            if UnitInfo.Info(k).stim_duration == 200
                spiketable = load([directory filesep UnitInfo.List{1,k}]);
                disp(UnitInfo.List{1,k})
                
                if ~isempty(spiketable)
                    
                    PREstim_duration = UnitInfo.Info(k).Pre_stim_Duration/1000; %seconds
                    stim_duration = UnitInfo.Info(k).stim_duration;
                    TrialLength = UnitInfo.Info(k).Pre_stim_Duration + UnitInfo.Info(k).Post_stim_Duration +stim_duration;
                    StimNb = unique(spiketable(:,1)).';
                    if max(StimNb) ==9
                        StimNb = 2:10;
                    end
                    
                    raster.stim=[];  raster.rep=[];  raster.spikes=[];
                    spikes_pooled_for_VS = [];
                    all_mean_rate_spont = [];
                    all_std_rate_spont = [];
                    all_mean_rate_stim = [];
                    all_std_rate_stim = [];
                    %     all_rate_total = [];
                    Fanofactor = [];
                    emptycount = 0;
                    
                    
                    f = SyncInfo.Info(n).attn*1e-1;
                    
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
                            
                            raster.stim=[raster.stim ones(size(spikes))];
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
                        
                        %Extracting individual neurons.
                        output.rates_stim{nn,1} = rate_total(:,PREstim_duration*1000+1:PREstim_duration*1000+stim_duration+100);
                        output.rates_pre{nn,1} = rate_total(:,PREstim_duration*1000-99:PREstim_duration*1000);
                        
                    end
                    
                    % plotting individual neurons
                    xs = 1:400;
                    h = 10;
                    for iii = 1:400
                        ys(iii) = gaussian_kern_reg(xs(iii),xs,mean([output.rates_pre{nn,1} output.rates_stim{nn,1}],1),h);
                    end
                    subplot(6,5,nn)
                    plot(xs,ys,'LineWidth',1.7)
                    xlabel(UnitInfo.List{1,k})
                    pause(0.1)
                end
                nn = nn+1;
                test = 1;
            end
        end
    end
    
end


output.meanDR = [];
output.errorDR = [];
output.mean_rate_stim = [];
output.std_rate_stim = [];
output.mean_rate_pre = [];
output.std_rate_pre = [];
output.var_isi = [];
rate_stim_total = [];
rate_pre_total = [];
for n = 1:nn-1
    rate_stim_total = [rate_stim_total; output.rates_stim{n}];
    rate_pre_total = [rate_pre_total; output.rates_pre{n}];
end





output.mean_rate_stim = mean(rate_stim_total,1);
output.std_rate_stim = std(rate_stim_total,0,1);
output.mean_rate_pre = mean(rate_pre_total,1);
output.std_rate_pre = std(rate_pre_total,0,1);
%     output.mean_rate_post{f} = mean(output.rates_post{f},1);
%     SpikeCount = sum(output.rates_stim{f},2)/1000.;
%     output.Fanofactor = [output.Fanofactor std(SpikeCount)^2/mean(SpikeCount)];
%     output.var_isi = [output.var_isi std(output.isi_total{f})];

%Gaussian smoothing
xs = 1:400;
h = 10;
for i = 1:400
    ys(i) = gaussian_kern_reg(xs(i),xs,[output.mean_rate_pre output.mean_rate_stim],h);
end
% smooth_rate = ys;
output.meanDR = [ output.meanDR mean2(output.mean_rate_stim)];
output.errorDR = [output.errorDR std2(output.mean_rate_stim)/sqrt(size(output.mean_rate_stim,1)*size(output.mean_rate_stim,2))];

% figure
plot(xs,ys,'LineWidth',1.7)

% figure
% cmap = colormap(jet(length(Attn_List)));
% for f = 1: length(Attn_List)
%     plot(smooth_rate,'color',cmap(f,:),'LineWidth',1.7,'DisplayName', ...
%         [num2str(Attn_List(f)) 'dB'] );
%     %     lgd{f} = [num2str(ceil(1000/ICI_list(f))) 'Hz'];
%     hold on
% end
% legend('show')

    
    
    
    
% directory = ['C:\Users\John\Dropbox\Bendorlab\Code\Marmoset\' animal];

%% Create fileinfo
function UnitInfo = createFileInfo(animal)

% animal = 'm2p';

global directory
% directory = ['U:\Neural and Behavioural Data' '\marmoset\' animal];
% directory = ['C:\Users\John\Dropbox\Bendorlab\Code\Marmoset\' animal];
    directory = ['C:\Users\John.Lee\OneDrive\Bendorlab\Marmoset\' animal];

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










