% running TIMMA test for different cell lines.
% this is the function to find the optimal find by k_max, S ,y

function[minInd] = TIMMATest(S, y, k_max)
[d,k] = size(S);
error = zeros(1,k);
for i= 0:k-1
    disp(i)
    initial_list = [zeros(1,i),1,zeros(1,k-1-i)];
    [~,err,~] = TIMMA_floating(k,d, [], initial_list, S, y, k_max,[],[],[]);
    error(i+1) = mean(err);
   
end
[~, minInd] = min(error);