spgl1_path = genpath('..');
addpath(spgl1_path);

data = load('gps_data_1.mat');
xhat = imresize(data.dlon,1/8,'bilinear');
yhat = imresize(data.dlat,1/8,'bilinear');
npts = numel(xhat);

n = npts;
integrator_mtrx = gallery('triw',npts,1,npts)';
H    = integrator_mtrx * integrator_mtrx;
W    = eye(npts,npts);
D    =  [ H , 0*H ; 0*H, H ];
xhat_noisy(1) = 0;
yhat_noisy(1) = 0;
b    =  [ xhat;  yhat];
a_LS = D \ b;

group= [1:npts , 1:npts]' ;

foo = spgSetParms();
foo.verbosity  = 2;
foo.iterations = 10000;
foo.iter_skip  = 20;
sval           = 5e-2*norm(b);
[X,R,G,INFO]   = spg_group(D,b,group,sval,foo);
disp(INFO);

L1_x = X(1:end/2); L1_y = X(end/2+1:end);
abs_acc = abs(L1_x)+abs(L1_y);
sfigure(3); stem( abs_acc); hold on;
legend('L_1 reconstructed |acceleration|'); 
hold off;
x_recons = H * L1_x;
y_recons = H * L1_y;
sfigure(1); clf; hold on; 
plot( xhat, yhat, '-');   
plot( x_recons, y_recons, 'r--' );  legend('measured','reconstructed');
% TODO: need a 'moving max' function, with 
% "minimum signal level" to avoid the near-zero stretches

