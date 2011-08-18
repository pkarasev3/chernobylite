npts    = 40;
xhat    = zeros( npts , 1 );
yhat    = zeros( npts , 1 );
e_sigma = 1e-2;
u = randn(1,1);
v = randn(1,1);
normuv = sqrt( u^2 + v^2 );
u = u / normuv;
v = v / normuv;
sum2switch = 0;
idx_switch = [];
for k = 2:npts
  
  lambda = 1;
  sum2switch = sum2switch + poissrnd(lambda,1,1);
  
  if( sum2switch > npts/5 )
    idx_switch = [idx_switch, k-1];
    sum2switch = 0;
    u = 0.5 * u + 0.5 * randn(1,1);
    v = 0.5 * v + 0.5 * randn(1,1);
    normuv = sqrt( u^2 + v^2 );
    u = u / normuv;
    v = v / normuv;
  end
  xhat(k) = xhat(k-1) + u;
  yhat(k) = yhat(k-1) + v;
  
end


xhat_obs = xhat + randn(size(xhat))*e_sigma;
yhat_obs = yhat + randn(size(yhat))*e_sigma;
sfigure(1); clf; hold on;
plot( xhat, yhat, '-');   hold on;
for k = 1:npts
  
  cval = [1.0-k/npts, abs(0.5-k/npts), k/npts];
  plot( xhat_obs(k), yhat_obs(k), '.','color',cval);   hold on;
  if( ~isempty(find(k==idx_switch, 1) ) )
    plot( xhat(k), yhat(k), '*','color',cval);   hold on;
    plot( xhat_obs(k), yhat_obs(k), '*','color',cval);   hold on;
  end
end
sfigure(1); axis equal; hold off;

n = npts;
e = ones(n,1);
A = spdiags([e -2*e e], -1:1, n, n);
A(1,:)         = 0*A(1,:); 
A(end,:)       = 0*A(end,:); 



ax= A * xhat_obs;
ay= A * yhat_obs;
sfigure(2); clf; hold on; 
plot(ax,'r.'); plot(ay,'b.');   % plot all
plot(idx_switch,ax(idx_switch),'mx','MarkerSize',12); 
plot(idx_switch,ay(idx_switch),'cx','MarkerSize',12);  % plot switch points
hold off; 
title('A * x, A * y'); hold off;

