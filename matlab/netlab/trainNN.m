%VideoType: string containing 'smoke' or 'fire'
%iterations: number of iterations to be performed for training, ~ 3000 is
%            the experimentally derived stopping criteria
%step: number of iterations of scg between showing user current error, must
%            be large enough(~30 is sufficient)
%nhidden: number of hidden units in the 2 layers neural network
%act_fn: string containing 'linear', 'logistic', or 'softmax' (usually
%        'logistic' is choosen


function trainNN(VideoType, iterations,step, nhidden,act_fn,result) 

num_samples=length(result.target);
index_array=round(rand(round(num_samples/5),1)*num_samples);  %pick random 20% of data
index_array(find(index_array==0))=1;        %check index is in bounds
index_array(find(index_array>num_samples))=num_samples;  %check index is in bounds
help_array=zeros(num_samples,1);
help_array(index_array)=1;
train_ind=find(help_array==0);
test_ind=find(help_array~=0);
all_data=result.inputdata;
all_target=result.target;
x=all_data(train_ind,:);
t=all_target(train_ind,:);
    
x_train=all_data(test_ind,:);
t_train=all_target(test_ind,:);



%setup neural network architecture
net = mlp(size(x, 2), nhidden, size(t, 2), act_fn); % Create net and find initial error
net = mlpinit(net, 10); % Initialise network with inverse variance of 10

fprintf('---------Training started--------- \n')

h=1;
cross_validate_error=zeros(ceil(iterations/step),1);
for i = [step-1:step:iterations iterations]
    [net, error] = mlptrain(net, x, t, step);
    y_train = mlpfwd(net, x_train);
    y_train =round(y_train);
    cross_validate_error(h)=sum(abs(y_train-t_train))/length(t_train)*100;   
    h=h+1;
    fprintf('Iterations completed= %d/%d, Current error= %d \n',i,iterations, error);
end;

save cross_validate_error.mat cross_validate_error

%figure out not having which input makes the largest difference
h=1;
for c=1:size(x_train,2)
    x_train_nulled=x_train;
    x_train_nulled(:,c)=0;
    y_train = mlpfwd(net, x_train_nulled);
    y_train =round(y_train);
    error_min_info(h)=sum(abs(y_train-t_train))/length(t_train)*100; 
    h=h+1;
end
save error_min_info.mat error_min_info
%above block optional

if strcmpi(VideoType, 'smoke')
    save mlptrainSmoke.net net%%save the network weights to a "net" data structure
elseif strcmpi(VideoType, 'fire')
    save mlptrainFire.net net%%save the network weights to a "net" data structure
end

fprintf('---------Training Completed--------- \n')

end
