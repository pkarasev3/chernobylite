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
A(1,:)          = 0*A(1,:);   A(1,1) = 1;
A(end,:)        = 0*A(end,:); A(end,end) = 1;

ax = A * xhat - [zeros(npts-1,1); xhat(end) ];
ay = A * yhat - [zeros(npts-1,1); yhat(end) ];
assert( 1e-6 > abs( ax(end) ) );
assert( 1e-6 > abs( ay(end) ) );

% The 'realistic' observed values
H     = inv(A);  % funky:  imagesc( H );
xhat0 = A \ (ax + [zeros(npts-1,1); xhat(end) ] ); assert( norm(xhat0 - xhat ) < 1e-6 );
yhat0 = A \ (ay + [zeros(npts-1,1); yhat(end) ] ); assert( norm(yhat0 - yhat ) < 1e-6 );




sfigure(2); clf; hold on; 
plot(sqrt(ax.^2+ay.^2),'r.');   % plot all
plot(idx_switch,sqrt(ax(idx_switch).^2+ay(idx_switch).^2),'mx','MarkerSize',12); 
hold off; 
title('ax, ay'); hold off;

