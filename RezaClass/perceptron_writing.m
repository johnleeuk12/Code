function perceptron_writing()

clear all

%% data extraction. 
% Here I save each training example for each number (from 0 to 9) in
% separate cells
for n = 0:9
    X{n+1} = {};
    fid = fopen(['data' num2str(n)],'r');
    for tr = 1:1000
    [dataset{n+1}{tr},~] = fread(fid,[28,28],'uchar');
    
    end
    X{n+1} = zeros(3,1000); %initializing feature matrix
end
fclose('all')

%% data 

for n = 0:9
    for tr = 1:1000
        % Here I convert pixels with values higher than one to '1' and all the
        % other pixels to '0'
        dataset{n+1}{tr}([find(dataset{n+1}{tr}<100)]) =0;
        dataset{n+1}{tr}([find(dataset{n+1}{tr}>=100)]) =1;
        
        %Here I extract three features from each training example
        X{n+1}(1,tr) = mean2(dataset{n+1}{tr});                                                      % total average pixel count
        X{n+1}(2,tr) = sum(sum((dataset{n+1}{tr}(11:20,11:20))))/100;                                % average pixel count at center
        X{n+1}(3,tr) = (length(find(dataset{n+1}{tr}(1:14,:) == dataset{n+1}{tr}(15:28,:))) + ...
            length(find(dataset{n+1}{tr}(:,1:14) == dataset{n+1}{tr}(:,15:28))))/(28*28);            % percentage of symetry over x and over y axis.
    end
end

Inputdata = [X{1} X{2}];
Outputdata = [zeros(1,1000) ones(1,1000)];
figure

scatter3(X{1}(1,:),X{1}(2,:),X{1}(3,:)) 
hold on 
scatter3(X{2}(1,:),X{2}(2,:),X{2}(3,:))
hold off
legend('0','1')
title('distribution of sample images in feature space')
xlabel('feature_1')
ylabel('feature_2')
zlabel('feature_3')

%% learning
w = zeros(2000,3);
error = zeros(1,500);
for step = 1:500
    indlist = randsample(1:2000,200);              %randomly select 200 image samples
    for tr = randsample(1:2000,2000)
        yhat(tr) = sign(w(tr,:)*Inputdata(:,tr));  %calculating estimation of y
        if ismember(tr,indlist)==1                 %for 200 randomly selected images, update estimation of y if the estimation is wrong
            if yhat(tr) ~=Outputdata(tr)
           w(tr,:) = w(tr,:) +  Outputdata(tr)*Inputdata(:,tr).';
           error(step) = error(step)+1;            %calculate number of wrong predictions for each training step.
            end
        end
    end
end
    
figure

plot(error)
title('the number of wrong predictions per each step')
ylabel('number of errors')
xlabel('steps')

            
        




