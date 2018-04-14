animal = 'm2p';
directory = ['U:\Neural and Behavioural Data' '\marmoset\' animal];
filelist = dir(directory);
filelist = filelist(3:end,:);
filelist = extractfield(filelist,'name');
UnitInfo = {};
UnitInfo.List = []; i = 1;
UnitInfo.Info = [];
% UnitInfo.Info = struct('Channel_Nb',{},'Stimuli_Nb',{},'Pre_stim_Duration',{},'Post_stim_Duration',{},'Sync',{},'Significant_rate',{},'Positive',{},'Positive_sig',{}); for click trains 
UnitInfo.Info = struct('Channel_Nb',{},'tone_Hz',{},'Attn_dB',{},'Pre_stim_Duration',{},'Post_stim_Duration',{},'Stim_duration',{});

filelist{1} = 'M2p0010.dat';
for k = 1 %:length(filelist)
    fid = fopen(['U:\Neural and Behavioural Data' filesep 'marmoset' filesep animal filesep filelist{k}]);
    frewind(fid)
    fileinfo = textscan(fid,'%s %s',20,'delimiter',{' = ','% '}...
        ,'MultipleDelimsAsOne',1,'headerlines', 6);
    fclose(fid);
    A = fileread(['U:\Neural and Behavioural Data' filesep 'marmoset' filesep animal filesep filelist{k}]);
    A1 = strfind(A,'Pre-stimulus record time (ms) = ');
    A2 = strfind(A,'Post-stimulus record time (ms) = ');
    l1 = length('Pre-stimulus record time (ms) = ');
    l2 = length('Post-stimulus record time (ms) = ');
    Prestim = str2num(A(A1+l1:A1+l1+2));
    Poststim = str2num(A(A2+l2:A2+l2+2));
    
end