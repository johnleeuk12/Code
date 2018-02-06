function Perceptron()

clear all
%% initializing data
load('input.dat')
load('output.dat')
w = randn(size(input,2),1);

%% learning
error = zeros(1,1000);
for step = 1:1000
    for n = 1:length(output)
        yhat = sign(input(n,:)*w);                      %calculate the estimation of y 
        if output(n,1) ~= 0
            if yhat ~= output(n,1)
                w = w + output(n,1)*input(n,:).';       %if the estimation of y is incorrect, update the weight vector w
                error(step) = error(step)+1;            %calculate number of wrong predictions for each training step.
            end
        else
            error(step) = error(step)+1;                %In case of prediction of 0, count that as an error.
        end
    end
end

figure
plot(error)
title('the number of wrong predictions per each step')
ylabel('number of errors')
xlabel('steps')