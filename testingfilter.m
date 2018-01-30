
animal = 'm2p';
directory = ['U:\Neural and Behavioural Data' '\marmoset\' animal];
filelist = dir(directory);
filelist = filelist(3:end,:);
filelist = extractfield(filelist,'name');
UnitInfo = {};
UnitInfo.List = []; i = 1;
UnitInfo.Info = [];
UnitInfo.Info = struct('Channel_Nb',{},'Stimuli_Nb',{},'Pre_stim_Duration',{},'Post_stim_Duration',{},'Sync',{},'Significant_rate',{},'Positive',{},'Positive_sig', {});

for k = 1:length(filelist)
    fid = fopen(['U:\Neural and Behavioural Data' filesep 'marmoset' filesep animal filesep filelist{k}]);
    frewind(fid)
    fileinfo = textscan(fid,'%s %s',20,'delimiter',{' = ','% '}...
        ,'MultipleDelimsAsOne',1,'headerlines', 6);
    fclose(fid);
    A = fileread(['U:\Neural and Behavioural Data' filesep 'marmoset' filesep 'm2p' filesep filelist{k}]);
    A1 = strfind(A,'Pre-stimulus record time (ms) = ');
    A2 = strfind(A,'Post-stimulus record time (ms) = ');
    l1 = length('Pre-stimulus record time (ms) = ');
    l2 = length('Post-stimulus record time (ms) = ');
    Prestim = str2num(A(A1+l1:A1+l1+2));
    Poststim = str2num(A(A2+l2:A2+l2+2));
    
    if strcmp(fileinfo{1,2}{1,1}, '1111')==1
        if strcmp(fileinfo{1,2}{16,1}, '12')==1 || strcmp(fileinfo{1,2}{16,1}, '20')==1
            %         pause
            UnitInfo.List{1,i} = filelist{k};
            UnitInfo.Info(i).Channel_Nb = str2num(fileinfo{1,2}{10,1});
            UnitInfo.Info(i).Stimuli_Nb = str2num(fileinfo{1,2}{16,1});
            UnitInfo.Info(i).Pre_stim_Duration = Prestim;
            UnitInfo.Info(i).Post_stim_Duration = Poststim;
            i = i+1
        elseif strcmp(fileinfo{1,2}{17,1}, '12')==1 || strcmp(fileinfo{1,2}{17,1}, '20')==1
            UnitInfo.List{1,i} = filelist{k};
            UnitInfo.Info(i).Channel_Nb = str2num(fileinfo{1,2}{10,1});
            UnitInfo.Info(i).Stimuli_Nb = str2num(fileinfo{1,2}{17,1});
            UnitInfo.Info(i).Pre_stim_Duration = Prestim;
            UnitInfo.Info(i).Post_stim_Duration = Poststim;
            i = i+1
        end
    end
end