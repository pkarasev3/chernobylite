npts    = 100;

e_sigma = 10e-1;
u = randn(1,1);
v = randn(1,1);
normuv = sqrt( u^2 + v^2 );
u = u / normuv;
v = v / normuv;
sum2switch = 0;
idx_switch = [];

% TODO: for a fixed set of accelerations,
% run in batch over parameter settings! 

% TODO: make the function take npts, lambda, and return the accelerations


% A micro-test! 
  integrator_mtrx = gallery('triw',10,1,10)';
  aa = [0;0;1;1;0;0;1;-1;0;0];
  vv = integrator_mtrx * aa;
  xx = integrator_mtrx * vv;
%

n = npts;
integrator_mtrx = gallery('triw',npts,1,npts)';
ax = zeros(npts,1);
ay = zeros(npts,1);
delta_uv_prv = [0;0];
for k = 2:npts
  
  lambda = 1;
  sum2switch = sum2switch + poissrnd(lambda,1,1);
  
  if( sum2switch > npts/5 )
    idx_switch = [idx_switch, k-1];
    sum2switch = 0;
    delta_uv=zeros(2,1);
    
    % enforce a minimum angle turn
    while( norm(delta_uv)<1e-1 || abs( delta_uv'*delta_uv_prv ) > 0.3 )
      delta_uv = randn(2,1); 
      delta_uv = delta_uv / norm(delta_uv);
    end
    ax(k) = delta_uv(1);
    ay(k) = delta_uv(2);
    delta_uv_prv = delta_uv + delta_uv_prv;
  end
  
end

vx     = integrator_mtrx * ax;
vy     = integrator_mtrx * ay;
xhat   = integrator_mtrx * vx;
yhat   = integrator_mtrx * vy;

% TODO: implement H as row * RHS operations, not full dense matrix

% The 'realistic' observed values
H     = integrator_mtrx*integrator_mtrx;
xhat0 = H * (ax ); assert( norm(xhat0 - xhat ) < 1e-1 );
yhat0 = H * (ay ); assert( norm(yhat0 - yhat ) < 1e-1 );

h          = fspecial('gaussian',[3 1],1.0);
xhat_noisy = imfilter(xhat0 + e_sigma * randn(size(xhat0)),h,'replicate');
yhat_noisy = imfilter(yhat0 + e_sigma * randn(size(xhat0)),h,'replicate');
xhat_obs = xhat_noisy;
yhat_obs = yhat_noisy;


sfigure(1); clf; hold on;
plot( xhat, yhat, '-');   hold on;
for k = 1:npts
  
  cval = [1.0-k/npts, abs(0.5-k/npts), k/npts];
  plot( xhat_obs(k), yhat_obs(k), '.','color',cval);  
  legend('true path (line)','measurements (points)'); hold on;
  if( ~isempty(find(k==idx_switch, 1) ) )
    plot( xhat(k), yhat(k), '*','color',cval);   hold on;
    plot( xhat_obs(k), yhat_obs(k), '*','color',cval);   hold on;
  end
end
sfigure(1); axis equal; hold off;

W    = eye(npts,npts);
D    =  [ H , 0*H ; 0*H, H ];
xhat_noisy(1) = 0;
yhat_noisy(1) = 0;
b    =  [ xhat_noisy;  yhat_noisy];
a_LS = D \ b;

group= [1:npts , 1:npts]' ;

foo = spgSetParms();
foo.verbosity  = 2;
foo.iterations = 5000;
foo.iter_skip  = 20;
sval           = 10;
[X,R,G,INFO]   = spg_group(D,b,group,sval,foo);
disp(INFO);

L1_x = X(1:end/2); L1_y = X(end/2+1:end);
sfigure(3); stem(abs(L1_x)+abs(L1_y)); hold on; stem(abs(ax)+abs(ay),'r'); 
legend('L_1 reconstructed |acceleration|', 'true |acceleration|'); hold off;









