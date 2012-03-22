function classifier = trainNetlabNN(featureVec, targetVec, params,net0) 
%train Netlab Neural Network
%
%INPUT:
% - featureVec: double array
%       + cols: features
%       + rows: pixels (first all pos pixels, then all neg pixels)
% - targetVec: double column vector [1;1;....;1;0;0;....;0] with 1
%              belonging to pos pixels and 0 to neg pixels
% - params: cell with entries
%   +iterations: number of random ski slopes to start gradient descent from
%                default 500
%   +step: scaled conjugate gradient iterations ( steps down the slope from
%          starting point), default 5
%   +nhidden: number of hidden units in the 2 layers neural network
%             (optimization is ~quadratic or worse cost in this, but might capture
%             more structure with higher value. ), default 30
%   +act_fn: string containing 'linear', 'logistic', or 'softmax',
%            default 'logistic'
%
%OUTPUT: classifier, struct with fields
%   -net: Netlab neural net
%   -cross_validate_error

%default params
iterations = 1000;
step = 10;
nhidden = 30;
act_fn = 'logistic';
if ~isempty(params)
    iterations = params{1};
end
if length(params) > 1
    step = params{2};
end
if length(params) > 2
    nhidden = params{3};
end
if length(params) > 3
    act_fn = params{4};
end
if length(params) > 4
    error('too many parameters for trainNetlabNN, check u.TrainClassifierParams')
end

num_samples = length(targetVec);
%pick randomly about 20% of data for validation
index_array = round(rand(round(num_samples/5),1)*num_samples);  
%make sure index is between 1 and num_samples
index_array(index_array==0) = 1;        
index_array(index_array>num_samples)=num_samples;
%transform to logical map
val_ind = false(num_samples,1);
val_ind(index_array) = true;
train_ind = ~val_ind;
%define train and validation data vector
x_train=featureVec(train_ind,:);
t_train=targetVec(train_ind,:);
x_val=featureVec(val_ind,:);
t_val=targetVec(val_ind,:);

if( (nargin < 4) || (rand(1,1)<0.5 ) )
%setup neural network architecture
  net = mlp(size(x_train, 2), nhidden, size(t_train, 2), act_fn); % Create net and find initial error
  net = mlpinit(net, 10); % Initialise network with inverse variance of 10
else
  net  = net0;
  net.b1 = net.b1 .* (0.5+1e-1*randn( size( net.b1 ) ) );
  net.b2 = net.b2 .* (0.5+1e-1*randn( size( net.b2 ) ) );
  net.w1 = net.w1 .* (0.5+1e-1*randn( size( net.w1 ) ) );
  net.w2 = net.w2 .* (0.5+1e-1*randn( size( net.w2 ) ) );
end

%initialize loop
h=1;
cross_validate_error=zeros(ceil(iterations/step),1);
%sfigure(2); clf; 
% (Could use parfor)
for i = [step-1:step:iterations iterations]
    [net, error_] = mlptrain(net, x_train, t_train, step);
    y_val = mlpfwd(net, x_val);
    cross_validate_error(h)=sum(abs(y_val-t_val))/length(t_val)*100;   
    h=h+1;
    %sfigure(2); hold on; plot( i, error_, 'ro' );  
     if ~mod(h,5) || i==iterations
%         fprintf('Training iterations completed= %d/%d, Current error= %d \n',i,iterations, error_);
     %     drawnow;
     end
     
end
hold off; 
%assign to output
classifier.net=net;
classifier.cross_validate_error = cross_validate_error;



