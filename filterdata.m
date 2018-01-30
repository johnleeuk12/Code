function filterdata(n)


% will be needed when taking all data
%load('infoset.mat');
%animals = unique(cellfun(@char,{output.animal},'unif',0));

%% file info
animal_list = {'m36n','m2p','m41o','m32q'};

% createFileInfo('m2p') %animal namecode
animal = animal_list{n};

% UnitInfo = createFileInfo(animal);

load([animal '_List3']);

ICI_list1 = [2 2.5 3 5 7.5 10 12.5 15 20 25 30 35 40 45 50 55 60 65 70 75]; %ms, ICI
ICI_list2 = [250 125 83.3333 62.5 50 41.6667 35.7143 31.25 27.7778 25 22.7273 20.8333];
%

directory = ['U:\Neural and Behavioural Data' '\marmoset\' animal];
% for n = 1:length(UnitInfo.List)
%     disp(n)
%     spiketable = load([directory filesep UnitInfo.List{1,n}]);
%
% %     load info on recording
%         if UnitInfo.Info(n).Stimuli_Nb == 12
%             ICI_list = ICI_list2;
%         else % UnitInfo.Info(2,n) == 20
%             ICI_list = ICI_list1;
%         end
%     PREstim_duration = UnitInfo.Info(n).Pre_stim_Duration/1000; %seconds
%     stim_duration = 0.5;
%         TrialLength = UnitInfo.Info(n).Pre_stim_Duration + UnitInfo.Info(n).Post_stim_Duration +500;
%     TrialLength = UnitInfo.Info(n).Pre_stim_Duration + UnitInfo.Info(n).Post_stim_Duration +200;
%     NbOstim = UnitInfo.Info(n).Stimuli_Nb;
%
%     %     vector = zeros(1,length(ICI_list));
%     raster.stim=[];  raster.rep=[];  raster.spikes=[];
%     spikes_pooled_for_VS = [];
%     all_mean_rate_spont = [];
%     all_std_rate_spont = [];
%     all_mean_rate_stim = [];
%     all_std_rate_stim = [];
%     %     all_rate_total = [];
%     Fanofactor = [];
%     emptycount = 0;
%     isi = [];
%     var_isi = [];
%         for f = 1:length(ICI_list)
%             freq = (1e3)./ICI_list(f);
%     for f = 1:NbOstim
%          indstim = find(spiketable(:,1)==f);
%         nbrep = spiketable(find(spiketable(:,1)==f),2);
%         nreps = max(nbrep.');
%         spikes_pooled = [];
%         rate_total = [];
%         if ~isempty(nbrep)
%             for r = unique(nbrep.')
%                 spikes1 = []; %channel 1
%                 spikes2 = []; %channel 2
%                 for i = indstim.'
%                     if spiketable(i,2) == r && spiketable(i,3)== 1 && spiketable(i,4) > 0 % channel 1
%                         spikes1 = [spikes1 (spiketable(i,4)*1e-6 - PREstim_duration)]; % spikes rounded up to 0.1ms cad 100microseconds. converted to seconds
%                     end
%                     if spiketable(i,2) == r && spiketable(i,3)== 2 && spiketable(i,4) > 0 % channel 2
%                         spikes2 = [spikes2 (spiketable(i,4)*1e-6 - PREstim_duration)];
%                     end
%                 end
%                 if length(spikes1) > length(spikes2)
%                     spikes =  spikes1;
%                 else
%                     spikes = spikes2;
%                 end
%                 if isempty(spikes(find(spikes>PREstim_duration & spikes< PREstim_duration + 0.2)))
%                     emptycount = emptycount+1;
%                 end
%
%                 spikes_pooled = [spikes_pooled spikes];
%
%                 raster.stim=[raster.stim f*ones(size(spikes))];
%                 raster.rep=[raster.rep r*ones(size(spikes))];
%                 raster.spikes=[raster.spikes spikes];
%
%                 % rate
%                 rate = zeros(1,TrialLength);
%                 spikes4rate = spikes(find(spikes<PREstim_duration+0.2)) + PREstim_duration;
%
%
%                 for st = spikes4rate
%                     if ceil(st*1e3) <= length(rate)
%                         rate(1,ceil(st*1e3)) = rate(1,ceil(st*1e3))+1;
%                     end
%
%                 end
%                 rate_total = [rate_total ; rate*1000];
%
%
%             end
%         end
%         rate_av = mean(rate_total,1);
%         xs = 1:TrialLength;
%         h = 15; %kernal bandwidth. determines the shape of the function
%         for i = 1:TrialLength
%             ys(i)=gaussian_kern_reg(xs(i),xs,rate_av,h);
%         end
% %         plot(ys)
%     end
% end
%





%% spikes, VS and raster
directory = ['U:\Neural and Behavioural Data' '\marmoset\' animal];
% directory = ['C:\Users\John\Dropbox\Bendorlab\Code' '\marmoset\' animal];

for n = 1:length(UnitInfo.List)
    UnitInfo.Info(n).Sync = 0;
    UnitInfo.Info(n).Significant_rate = 0;
    UnitInfo.Info(n).Positive = 0;
end
%
for n = 1:length(UnitInfo.List)
    disp(n)
    spiketable = load([directory filesep UnitInfo.List{1,n}]);
    
%     load info on recording
    if UnitInfo.Info(n).Stimuli_Nb == 12
        ICI_list = ICI_list2;
    else % UnitInfo.Info(2,n) == 20
        ICI_list = ICI_list1;
    end
    PREstim_duration = UnitInfo.Info(n).Pre_stim_Duration/1000; %seconds
    stim_duration = 0.5;
    TrialLength = UnitInfo.Info(n).Pre_stim_Duration + UnitInfo.Info(n).Post_stim_Duration +500;
    %     TrialLength = UnitInfo.Info(n).Pre_stim_Duration + UnitInfo.Info(n).Post_stim_Duration +UnitInfo.Info(n).Stim_duration;
    NbOstim = UnitInfo.Info(n).Stimuli_Nb;
    
    vector = zeros(1,length(ICI_list));
    raster.stim=[];  raster.rep=[];  raster.spikes=[];
    spikes_pooled_for_VS = [];
    all_mean_rate_spont = [];
    all_std_rate_spont = [];
    all_mean_rate_stim = [];
    all_std_rate_stim = [];
    all_rate_total = [];
    all_empty = [];
    Fanofactor = [];
    emptycount = 0;
    isi = [];
    var_isi = [];
    for f = 1:length(ICI_list)
        freq = (1e3)./ICI_list(f);
%     end
%     for f = 1:NbOstim
        
        
        indstim = find(spiketable(:,1)==f);
        nbrep = spiketable(find(spiketable(:,1)==f),2);
        nreps = max(nbrep.');
        spikes_pooled = [];
        rate_total = [];
        empty = 0;
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
                
                if isempty(spikes4rate)
                    empty = empty+1;
                end
                for st = spikes4rate
                    if ceil(st*1e3) <= length(rate)
                        rate(1,ceil(st*1e3)) = rate(1,ceil(st*1e3))+1;
                    end
                    
                end
                rate_total = [rate_total ; rate*1000];
                
                
            end
            if empty*2 >nreps
                all_empty(f) = 0;
            end    
            
            %calculating Fano factor
            SpikeCount = sum(rate_total(:,PREstim_duration*1000+1:PREstim_duration*1000+600),2)/1000.;
            Fanofactor = [Fanofactor std(SpikeCount)^2/mean(SpikeCount)];
            
            %average rate + gaussian kernel
            rate_av = mean(rate_total,1);
            xs = 1:TrialLength;
            h = 15; %kernal bandwidth. determines the shape of the function
            for i = 1:TrialLength
                ys(i)=gaussian_kern_reg(xs(i),xs,rate_av,h);
            end
            all_rate{f} = ys;
            %             all_rate_total = [all_rate_total ; rate_total];
            all_mean_rate_spont(f) = mean2(rate_total(:,1:PREstim_duration*1000));
            all_std_rate_spont(f) = std2(rate_total(:,1:PREstim_duration*1000))/sqrt(PREstim_duration*1000*nreps);
            all_mean_rate_stim(f) = mean2(rate_total(:,PREstim_duration*1000+1:PREstim_duration*1000+600));% -all_mean_rate_spont(f);
            all_std_rate_stim(f) = std2(rate_total(:,PREstim_duration*1000+1:PREstim_duration*1000+600))/sqrt(600*nreps);
            spikes_pooled_for_VS=spikes_pooled(find(spikes_pooled>0.05 & spikes_pooled<=(stim_duration+0.05)));
            
            %             %Interspike interval
            %             spikes1 = spikes_pooled_for_VS;
            %             spikes2 = [0 spikes1(1:end-1)];
            %             isi = spikes1-spikes2;
            %             var_isi = [var_isi std(isi)];
            
            %
            
            %for calculating vector strength, subtract the first 50 ms and
            %include 50 ms post stimulus
            
            
            if ~isempty(spikes_pooled_for_VS)
                total_spikes = length(spikes_pooled_for_VS);
                x = 0;
                y = 0;
                if total_spikes > 0
                    x = sum(cos(2*pi*(spikes_pooled_for_VS*freq)));
                    y = sum(sin(2*pi*(spikes_pooled_for_VS*freq)));
                    vector(f) = sqrt(x^2+y^2)/total_spikes;
                elseif total_spikes == 0
                    vector(f) = 0;
                end
                rayleigh(f) = 2*total_spikes*vector(f)^2;
            else
                vector(f)=0;
                rayleigh(f)=0;
            end
            
            if rayleigh(f)<13.8
                vector(f)=0;
            end
        end
        
    end
    
    %Testing for Sync counter1
    % Testing for significance rate counter2
    UnitInfo.Info(n).output.vector = vector;
    UnitInfo.Info(n).output.rayleigh = rayleigh;
    counter1 = 0;
    counter2 = 0;
    if ~isempty(nreps)
        if ~isempty(all_mean_rate_stim)
            if length(ICI_list) == 12
                new_ICI_list = ICI_list(2:end);
                new_all_mean_rate_stim = all_mean_rate_stim(2:end);
                indF = [2:length(ICI_list)];
                % This is to test monotonicity in flutter range. previously, 9
                % = 4
            else
                new_ICI_list = ICI_list(9:end);
                new_all_mean_rate_stim = all_mean_rate_stim(9:end);
                indF = [9:length(ICI_list)];
            end
             
            for f = indF
                if vector(f) >0.1
                    counter1 = counter1 + 1;
                else
                    counter1 = 0;
                end
                
                if all_mean_rate_stim(f) - 2*all_std_rate_stim(f)>0 && all_mean_rate_stim(f)-all_mean_rate_spont(f)>2
                    counter2 = counter2 + 1;
                else
                    counter2 = 0;
                end
                
                if ~isempty(all_empty)
                    counter2 = 0;
                end
                
                if counter2 == 3
                    UnitInfo.Info(n).Significant_rate = 1;% significant rate response
                end
                if counter1 == 3
                    UnitInfo.Info(n).Sync = 1; %Sync
                end
            end
            % testing for monotonicity take between 8 and 50Hz
            
            
            [RHO,PVAL] = corr(new_ICI_list.',new_all_mean_rate_stim.','Type','Spearman');
            
            %             if PVAL <0.05
            UnitInfo.Info(n).Rho = RHO;
            if RHO > 0.7 && PVAL <0.05
                if length(ICI_list) == 12
                    UnitInfo.Info(n).Positive = 1;
                else
                    UnitInfo.Info(n).Positive = -1;
                end
            elseif RHO < -0.7 && PVAL <0.05
                if length(ICI_list) == 12
                    UnitInfo.Info(n).Positive = -1;
                else
                    UnitInfo.Info(n).Positive = 1;
                end
            end
            %             else
            UnitInfo.Info(n).pval = PVAL;
            %             end
        end
        % Plot adaptation.
    end
    
    % UnitInfo.Info(n).FanoFactor
    % %
    %         plotdata(ICI_list, vector, raster, nreps, stim_duration,all_rate)
    %         figure
    %         shadedErrorBar(ICI_list,all_mean_rate_stim,all_std_rate_stim)
    
end
filename = [animal '_List3.mat'];
save(filename,'UnitInfo');


% test = 1
%% Create fileinfo
function UnitInfo = createFileInfo(animal)

% animal = 'm2p';
directory = ['U:\Neural and Behavioural Data' '\marmoset\' animal];
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
    fid = fopen(['U:\Neural and Behavioural Data' filesep 'marmoset' filesep animal filesep filelist{k}]);
    frewind(fid)
    fileinfo = textscan(fid,'%s %s',30,'delimiter',{' = ','% '}...
        ,'MultipleDelimsAsOne',1,'headerlines', 6);
    fclose(fid);
    A = fileread(['U:\Neural and Behavioural Data' filesep 'marmoset' filesep animal filesep filelist{k}]);
    A1 = strfind(A,'Pre-stimulus record time (ms) = ');
    A2 = strfind(A,'Post-stimulus record time (ms) = ');
    A3 = strfind(A,'Number of Stimuli = ');
    l1 = length('Pre-stimulus record time (ms) = ');
    l2 = length('Post-stimulus record time (ms) = ');
    l3 = length('Number of Stimuli = ');
    A4 = strfind(A,'pblaster_version = ');
    l4 = length('pblaster_version = ');
    
    
    Prestim = str2num(A(A1+l1:A1+l1+2));
    Poststim = str2num(A(A2+l2:A2+l2+2));
    NbOstimuli = str2num(A(A3+l3:A3+l3+1));
    pblaster =  str2num(A(A4+l4:A4+l4+2));
    %     if strcmp(fileinfo{1,2}{1,1}, '1')==1
    %         UnitInfo.List{1,i} = filelist{k};
    % %         UnitInfo.Info(i).Channel_Nb = str2num(fileinfo{1,2}{10,1});
    %         UnitInfo.Info(i).Stimuli_Nb = NbOstimuli;
    %         UnitInfo.Info(i).Pre_stim_Duration = Prestim;
    %         UnitInfo.Info(i).Post_stim_Duration = Poststim;
    %         i = i+1
    %     end
    % end
    
    if strcmp(fileinfo{1,2}{1,1}, '1111')==1 %This is for CLICK TRAINS
        if strcmp(fileinfo{1,2}{16,1}, '12')==1 % && pblaster
            if strcmp(fileinfo{1,2}{17,1}, '6.1')==1
                UnitInfo.List{1,i} = filelist{k};
                UnitInfo.Info(i).Channel_Nb = str2num(fileinfo{1,2}{10,1});
                UnitInfo.Info(i).Stimuli_Nb = str2num(fileinfo{1,2}{16,1});
                UnitInfo.Info(i).Pre_stim_Duration = Prestim;
                UnitInfo.Info(i).Post_stim_Duration = Poststim;
                i = i+1
            end
            
        elseif strcmp(fileinfo{1,2}{16,1}, '20')==1 %         pause
            UnitInfo.List{1,i} = filelist{k};
            UnitInfo.Info(i).Channel_Nb = str2num(fileinfo{1,2}{10,1});
            UnitInfo.Info(i).Stimuli_Nb = str2num(fileinfo{1,2}{16,1});
            UnitInfo.Info(i).Pre_stim_Duration = Prestim;
            UnitInfo.Info(i).Post_stim_Duration = Poststim;
            i = i+1
        elseif strcmp(fileinfo{1,2}{17,1}, '12')==1 % && pblaster
            if strcmp(fileinfo{1,2}{18,1}, '6.1')==1
                UnitInfo.List{1,i} = filelist{k};
                UnitInfo.Info(i).Channel_Nb = str2num(fileinfo{1,2}{10,1});
                UnitInfo.Info(i).Stimuli_Nb = str2num(fileinfo{1,2}{16,1});
                UnitInfo.Info(i).Pre_stim_Duration = Prestim;
                UnitInfo.Info(i).Post_stim_Duration = Poststim;
                i = i+1
            end
            
        elseif strcmp(fileinfo{1,2}{17,1}, '20')==1
            UnitInfo.List{1,i} = filelist{k};
            UnitInfo.Info(i).Channel_Nb = str2num(fileinfo{1,2}{10,1});
            UnitInfo.Info(i).Stimuli_Nb = str2num(fileinfo{1,2}{17,1});
            UnitInfo.Info(i).Pre_stim_Duration = Prestim;
            UnitInfo.Info(i).Post_stim_Duration = Poststim;
            i = i+1
        end
    end
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










