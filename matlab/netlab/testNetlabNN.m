function prob = testNetlabNN(classifier, Feature)
%classify based on feature vectors 'Feature' using netlab NN 'net'
%
%INPUT:
% - classifier: resTrain.classifier
% - Feature: feature vectors. One column for each feature
%
% OUTPUT:
% - prob: vector containing classifier output

net = classifier.net;
%classify
prob = mlpfwd(net,Feature);
