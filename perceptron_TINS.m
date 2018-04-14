function perceptron_TINS()

clear all

N = 100; % Nb of learning trials

total_error = {};
Accuracy = {};
final_data = {};
for state = 1:2;
    total_error{state} = zeros(1,320);
    Accuracy{state} = zeros(N,320);
    global_Accuracy{state} = zeros(1,32);
    global_error{state} = zeros(1,32);
    bin_Accuracy{state} = zeros(N,32);
    bin_y{state}{1} = zeros(N,32);
    bin_y{state}{2} = zeros(N,32);% y{1} = target
    final_data.weight{state} = [];
    %     final_data.y{state}{2} = [];
    %     final_data.y{state}{1} = [];
    %     y_error{state}{1} = zeros(1,32);
    %     y_error{state}{2} = zeros(1,32);
    final_data.z{state}{2} = [];
    final_data.z{state}{1} = [];
    z_error{state}{1} = zeros(1,32);
    z_error{state}{2} = zeros(1,32);
    
    weights{state} = [];
end

for tryout = 1:N
    disp(['trial ' int2str(tryout)])
    [output, data] = learning_TINS();
    
    for state = 1:2;
        final_data.weight{state} = [final_data.weight{state}; output.weight{state}];
        for ind = 1:2
            final_data.z{state}{ind} = [final_data.z{state}{ind}; output.z{state}{ind}];
        end
        %             total_error{state} = total_error{state} + error{state};
        Accuracy{state}(tryout,:) = 1-output.error{state}/2;
        weights{state}=[weights{state}; output.weight{state}];
    end
end


disp('Calculating accuracy...')


for state = 1:2;
    z_baseline{state} = (mean(mean(final_data.z{state}{1}(:,1:300))) + ...
        mean(mean(final_data.z{state}{2}(:,1:300))))/2;
end

% y_baseline{3} = (y_baseline{2} + y_baseline{1})/2;
for bin = 1:32
    for state = 1:2;
        for tryout = 1:N
            bin_Accuracy{state}(tryout,bin) = ...
                mean(Accuracy{state}(tryout,10*(bin-1)+1:10*bin));
            bin_y{state}{1}(tryout,bin) = ...
                mean(final_data.z{state}{1}(tryout,10*(bin-1)+1:10*bin)) ...
                - z_baseline{state};
            bin_y{state}{2}(tryout,bin) = ...
                mean(final_data.z{state}{2}(tryout,10*(bin-1)+1:10*bin)) ...
                - z_baseline{state};
        end
        global_Accuracy{state}(bin) = mean(bin_Accuracy{state}(:,bin));
        y_mean{state}{1}(:,bin) = mean(bin_y{state}{1}(:,bin));
        y_mean{state}{2}(:,bin) = mean(bin_y{state}{2}(:,bin));
        
        SEM_y_target = std(bin_y{state}{1}(:,bin))/sqrt(N);
        SEM_y_ref = std(bin_y{state}{2}(:,bin))/sqrt(N);
        SEM_acc = std(bin_Accuracy{state}(:,bin))/sqrt(N);
        ts = tinv([0.025, 0.975],N-1);
        global_error{state}(bin) = ts(2)*SEM_acc;
        y_error{state}{1}(bin) = ts(2)*SEM_y_target;
        y_error{state}{2}(bin) = ts(2)*SEM_y_ref;
    end
end


for state = 1:2 
    weights_mean{state} = mean(weights{state});
    ts = tinv([0.025, 0.975],N-1);
    for n = 1:359
    weights_error{state}(1,n) = ts(2)*sqrt(N)/std(weights{state}(:,n));
    end
end


final_data.Accuracy_mean = global_Accuracy;
final_data.Accuracy_error = global_error;
final_data.y_mean = y_mean;
final_data.y_error = y_error;


% analysis for alpha. data_analysis{1} = difference between Target and
% Reference, and data-analysis{2} = distance of Reference to baseline (0)
% mean and std for timebin = 24 : 28
ts = tinv([0.025, 0.975],3);
data_analysis{1}.mean = mean(abs(y_mean{2}{1}(1,24:28)-y_mean{2}{2}(1,24:28)));
data_analysis{1}.error = ts(2)*sqrt(3)/std(abs(y_mean{2}{1}(1,24:28)-y_mean{2}{2}(1,24:28)));
data_analysis{2}.mean = mean(abs(y_mean{2}{2}(1,24:28)));
data_analysis{2}.error = ts(2)*sqrt(3)/std(abs(y_mean{2}{2}(1,24:28)));
% Drawing figure

figure

shadedErrorBar(linspace(10,320,32),global_Accuracy{2},global_error{2},'-b') % active
hold on
shadedErrorBar(linspace(10,320,32),global_Accuracy{1},global_error{1},'-r') % passive
xlabel('ms');
ylabel('Accuracy')
axis([0,330,0.4,0.9])
save('perceptron_data.mat')

% plot(global_Accuracy{2},'b')
for x1 = [ 40 165 240 280 ]
    hold on
    line([x1 x1], get(gca, 'ylim'),'LineStyle','--');
end

hold on 
line(get(gca,'xlim'),[-100 100],'LineStyle','--');

figure 


shadedErrorBar(linspace(10,320,32),y_mean{1}{1},y_error{1}{1},'-r') % target
hold on
shadedErrorBar(linspace(10,320,32),y_mean{1}{2},y_error{1}{2},'-b') % ref

for x1 = [ 40 165 240 280 ]
    hold on
    line([x1 x1],[-150 150],'LineStyle','--');
end
xlabel('ms');
hold on 
line(get(gca,'xlim'),[0. 0.],'LineStyle','--');
axis([0,330,-150,150])

figure

shadedErrorBar(linspace(10,320,32),y_mean{2}{1},y_error{1}{1},'-r') % target
hold on
shadedErrorBar(linspace(10,320,32),y_mean{2}{2},y_error{1}{2},'-b') % ref

for x1 = [ 40 165 240 280 ]
    hold on
    line([x1 x1], [-150 150],'LineStyle','--');
end
xlabel('ms');
hold on 
line(get(gca,'xlim'),[0. 0.],'LineStyle','--');
axis([0,330,-150,150])

end

%% learning, individual learning trials

function [output, data] = learning_TINS()


Rootname = '/home/johnleeuk21/Documents/MATLAB/TINS/SortedPreActPost/';

% Rootname = 'C:\Users\John\Documents\MATLAB\TINS\SortedPreActPost\';

% Load SESSION NAMES
sessioninf = [];
load([Rootname filesep 'SessionInfo.mat'])
sessioninf.sessofint=[1:length(Session)];

data = [];
counter = 0;
for i = 1:2
    for j = 1:6
        for k = 1:7
            data{i}{j}{k} =[];
            % 1 = Training_target, 2 = Training_reference
            % 3 = Validation_target, 4 = validation_ref
            % 5 = Test_target, 6 = Test_ref
        end
    end
end

% Construct data
disp('Compiling Data...')
blacklist = '';
for SessionNum = sessioninf.sessofint;
    behave_state = Session{SessionNum}(end-6);
    SessionName=Session{SessionNum}(1:end-8);
    %     disp(SessionName)
    Animal = Session{SessionNum}(1:end-11);
    Filename = [Rootname Session{SessionNum}(1:end-8) '_' behave_state '_CLK' 'PSTHDSortedAt10ms.mat'];
    %     disp(SessionName)
    load(Filename)
    if sum(SpikeDat.Type == 1)>=15 && strcmp(blacklist,Animal) == 0;
        if behave_state =='a';
            [data, data_session] = construct_data(data,'', Filename, behave_state);
        else % second iteration, passive
            counter = counter +1;
            if mod(counter,2) == 0;
                [data, data_session] = construct_data(data, data_session, Filename, behave_state);
                if size(data{1}{1}{1},1) ~= size(data{2}{1}{1} ,1)
                    disp(Animal)
                end
            end
        end
    elseif sum(SpikeDat.Type == 1)<=15 
        blacklist = Animal;
    else counter = counter+1;
    end
end

% Learning phase

[output, data] = perceptron1(data);







end

% save variable called current animal trial (avo021) s

%% create database for perceptron training
function [data, data_session] = construct_data(data, data_session, Filename, behave_state)
%
% Rootname = '/home/johnleeuk21/Documents/MATLAB/TINS/SortedPreActPost/';
% Filename = 'avo021d18_p_CLKPSTHDSortedAt10ms.mat';
% behave_state = 'p';

load(Filename)
Nb_Neurons = size(SpikeDat.Average,2);
% disp(Filename)
% disp(Nb_Neurons);

nb_training_sample = 7;
indice.target = find(SpikeDat.Type==1); % find session numbers for correct target sessions
indice.ref = find(SpikeDat.Type==0); % find reference sessions
Training_ind.target = randsample(indice.target,nb_training_sample);
Training_ind.ref = randsample(indice.ref,nb_training_sample);

validation_ind.target = indice.target;
validation_ind.ref = indice.ref;

for ind = Training_ind.target;
    validation_ind.target = validation_ind.target(validation_ind.target ~=ind);
end
for ind1 = Training_ind.ref;
    validation_ind.ref = validation_ind.ref(validation_ind.ref ~=ind1);
end

Test_ind.target = randsample(validation_ind.target,1);
Test_ind.ref = randsample(validation_ind.ref,1);

validation_ind.target = validation_ind.target(validation_ind.target ~=Test_ind.target);
validation_ind.ref = validation_ind.ref(validation_ind.ref ~=Test_ind.target);
validation_ind.target = randsample(validation_ind.target,6);
validation_ind.ref = randsample(validation_ind.ref,6);


if behave_state == 'a' % initialise data to feed into perceptron
    state = 1;
    data_session = {};
    for i = 1:2
        for j = 1:4
            for k = 1:nb_training_sample
                data_session{i}{j}{k} = zeros(Nb_Neurons,320);
                % 1 = Training_target, 2 = Training_reference
                % 3 = Test_target,     4 = Test_target
            end
        end
        for j = 5:6
            data_session{i}{j}{1} = zeros(Nb_Neurons,320);
        end
        
    end
else state = 2; % passive state
end

% cross validation, separate training and test   set here

% take all sessions

for n = 1:Nb_Neurons
    for k = 1:nb_training_sample
        data_session{state}{1}{k}(n,:) = SpikeDat.Summed{Training_ind.target(k),n};
        data_session{state}{2}{k}(n,:) = SpikeDat.Summed{Training_ind.ref(k),n};
    end
    for k = 1:nb_training_sample-1
        data_session{state}{3}{k}(n,:) = SpikeDat.Summed{validation_ind.target(k),n};
        data_session{state}{4}{k}(n,:) = SpikeDat.Summed{validation_ind.ref(k),n};
    end
    data_session{state}{5}{1}(n,:) = SpikeDat.Summed{Test_ind.target,n};
    data_session{state}{6}{1}(n,:) = SpikeDat.Summed{Test_ind.ref,n};
end
% 
% for i = 1:2
%     for j = 1:2
%         for k = 2:6
%             data_session{i}{j}{1} = data_session{i}{j}{1} + data_session{i}{j}{k};
%         end
%         data_session{i}{j}{1} = data_session{i}{j}{1}/7;
%     end
% end
% 
% 
for j = 1:4
    for k = 1:nb_training_sample
    data{state}{j}{k} = [data{state}{j}{k} ; data_session{state}{j}{k}];
    end
end
for j = 5:6
    data{state}{j}{1} = [data{state}{j}{1} ; data_session{state}{j}{1}];
end

% mean of neuron activity over sessions



end


%% perceptron learning
function [output, data] = perceptron1(data)
disp('Learning...')
nb_iter = 1:30; % right now, completely useless
nb_training_sample = 7;
bias = 1;
period = 241:280; %165:240 sound, 241:280 silence after sound
% test_period = setdiff([165:240],training_period);
% active or passive
% First we will classify them with two different weights

weight = {};
coeff = 2; % learning rate

%training

% for state = [1 2]
%     Nb_Neurons{state} = size(data{state}{1}{1},1);
%     weight{state} = rand(1,Nb_Neurons{state}+1);
% end
% 
% c = linspace(1,10);
% for j = nb_iter
%     for state = [1 2] % here i = state
%         %         Nb_Neurons{state} = size(data{state}{1}{1},1);
% %         weight{state} = rand(1,Nb_Neurons{state}+1);
%         z = {};
%         estim = {};
%         error = {};
%         error{state} = zeros(1,320);
%         y{state}{1} = zeros(1,320);
%         y{state}{2} = zeros(1,320);
%         % start learning
%         
%         for t = 165:240
% 
%             for k = 1: nb_training_sample
%                 for ind = [1 2] % test_sets
%                     z1{ind} = dot(weight{state},[bias data{state}{ind}{k}(:,t).']);
%                     estim{ind} = 1/(1+exp(-z1{ind}));
%                     delta = abs(ind-2) - estim{ind}; % this way, target = 1, ref= 0
%                     weight{state} = weight{state} + coeff*delta*[bias data{state}{ind}{k}(:,t).'];
%                 end
% 
%             end
%             for ind = [5 6] % Test sets
%                 z{state}{ind-4}(1,t) = dot(weight{1},[bias data{state}{ind}{1}(:,t).']);
%                 y{state}{ind-4}(1,t) = 1/(1+exp(-z{state}{ind-4}(1,t)));
%                 if ind == 5 && y{state}{ind-4}(1,t) <0.5 %target but estim = ref
%                     error{state}(1,t) = error{state}(1,t) +1;
%                 elseif ind == 6 && y{state}{ind-4}(1,t) >0.5 % ref but estim = target
%                     error{state}(1,t) = error{state}(1,t) +1;
%                 end
%             end
%             if state == 1
%             disp(t)
%             axis([164,245,0,1])
%             scatter(t,sum(error{1})/(2*(t-164)),30,c(j),'filled')
%             hold on
%             end
%         end
%     end
% end

% for state = 1:2;
%     error{state} = zeros(1,320);
%     y{state}{1} = zeros(1,320);
%     y{state}{2} = zeros(1,320);
% end
% 
% 
% y = {}; % real estimation for testing set.
% % % here y{1} = target, y{2} = ref
% % z{3} = [];
% % z{4} = [];
% for state = [1 2]
%     
%     for t = 1:320
%         for ind = [5 6] % Test sets
%             z{state}{ind-4}(1,t) = dot(weight{1},[bias data{state}{ind}{1}(:,t).']);
%             y{state}{ind-4}(1,t) = 1/(1+exp(-z{state}{ind-4}(1,t)));
%             if ind == 5 && y{state}{ind-4}(1,t) <0.5 %target but estim = ref
%                 error{state}(1,t) = error{state}(1,t) +1;
%             elseif ind == 6 && y{state}{ind-4}(1,t) >0.5 % ref but estim = target
%                 error{state}(1,t) = error{state}(1,t) +1;
%             end
%         end
%     end
% end
% 
%     
% 
% 
% for state = 1:2;
%     error{state} = zeros(1,320);
%     y{state}{1} = zeros(1,320);
%     y{state}{2} = zeros(1,320);
% end
% 
%     
% y = {}; % real estimation for testing set. 
% % % here y{1} = target, y{2} = ref
% % z{3} = [];
% % z{4} = [];
% for state = [1 2]
% 
%     for t = 1:320
%         for ind = [5 6] % Test sets
%             z{state}{ind-4}(1,t) = dot(weight{1},[bias data{state}{ind}{1}(:,t).']);
%             y{state}{ind-4}(1,t) = 1/(1+exp(-z{state}{ind-4}(1,t)));
%             if ind == 5 && y{state}{ind-4}(1,t) <0.5 %target but estim = ref
%                 error{state}(1,t) = error{state}(1,t) +1;
%             elseif ind == 6 && y{state}{ind-4}(1,t) >0.5 % ref but estim = target
%                 error{state}(1,t) = error{state}(1,t) +1;
%             end
%         end
%     end
% end



for state = [1 2] % here i = state
    Nb_Neurons{state} = size(data{state}{1}{1},1);
    weight{state} = rand(1,Nb_Neurons{state}+1);
    z = {};
    estim = {};
% %     start learning
    for j = nb_iter
        for t = period
            for k = 1: nb_training_sample
                for ind = [1 2] % test_sets
                    z1{ind} = dot(weight{state},[bias data{state}{ind}{k}(:,t).']);
                    estim{ind} = 1/(1+exp(-z1{ind}));
                    delta = abs(ind-2) - estim{ind}; % this way, target = 1, ref= 0
                    weight{state} = weight{state} + coeff*delta*[bias data{state}{ind}{k}(:,t).'];
                end
            end
        end
    end
end



    
% calculating baseline
for state = [1 2]
        for t = 1:40
            for k = 1: nb_training_sample
                for ind = [1 2] % test_sets
                    z_base{state}{ind}(1,t) = dot(weight{state},[bias data{state}{ind}{k}(:,t).']);
                    y_base{state}{ind}(1,t) = 1/(1+exp(-z_base{state}{ind}(1,t)));
                end
            end
        end
end

for state = [1 2]
    for ind = [1 2]
        z_baseline{state}{ind} = mean(z_base{state}{ind});
        y_baseline{state}{ind} = mean(y_base{state}{ind});
    end
end

% Testing accuracy % valdation + updating weight

coeff1 = 2.; %active modifier
coeff2 = 0.25 ; %passive modifier

for j = 1:30
    error = {};
    for state = 1:2;
        error{state} = zeros(1,320);
        y{state}{1} = zeros(1,320);
        y{state}{2} = zeros(1,320);
    end
    
    for t = 1:320
        for state = [1 2]
            for ind = [3 4]
                for k = nb_training_sample-1
                    z{state}{ind-2}(1,t) = dot(weight{1},[bias data{state}{ind}{k}(:,t).']);
                    y{state}{ind-2}(1,t) = 1/(1+exp(-z{state}{ind-2}(1,t)));
                    if state ==1
                        delta = abs(ind-4)-y{state}{ind-2}(1,t);
                        weight{1} = weight{1} + coeff1*delta*[bias data{state}{ind}{k}(:,t).'];
                    else
                        delta = y_baseline{state}{ind-2}-y{state}{ind-2}(1,t);
                        weight{1} = weight{1} + coeff2*delta*[bias data{state}{ind}{k}(:,t).'];
                    end
%                     
                end
                if ind == 3 && y{state}{ind-2}(1,t) <0.5 %target but estim = ref
                    error{state}(1,t) = error{state}(1,t) +1;
                elseif ind == 4 && y{state}{ind-2}(1,t) >0.5 % ref but estim = target
                    error{state}(1,t) = error{state}(1,t) +1;
                end
            end
        end
    end
end

% Update accuracy and error 

error = {};
for state = 1:2;
    error{state} = zeros(1,320);
    y{state}{1} = zeros(1,320);
    y{state}{2} = zeros(1,320);
end

    
y = {}; % real estimation for testing set. 
% % here y{1} = target, y{2} = ref
% z{3} = [];
% z{4} = [];
for state = [1 2]

    for t = 1:320
        for ind = [5 6] % Test sets
            z{state}{ind-4}(1,t) = dot(weight{1},[bias data{state}{ind}{1}(:,t).']);
            y{state}{ind-4}(1,t) = 1/(1+exp(-z{state}{ind-4}(1,t)));
            if ind == 5 && y{state}{ind-4}(1,t) <0.5 %target but estim = ref
                error{state}(1,t) = error{state}(1,t) +1;
            elseif ind == 6 && y{state}{ind-4}(1,t) >0.5 % ref but estim = target
                error{state}(1,t) = error{state}(1,t) +1;
            end
        end
    end
end




output.error = error;
output.weight = weight;
output.y = y;
output.z = z;





end

