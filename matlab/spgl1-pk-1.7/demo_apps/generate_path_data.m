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
num_turns = (randn(1,1) > 0 ) * (-1)  + (randn(1,1) > 0 ) * (-1) +  4
for k = 2:npts
  
  lambda = 1;
  sum2switch = sum2switch + poissrnd(lambda,1,1);
  
  if( sum2switch > npts/num_turns )
    idx_switch = [idx_switch, k-1]; %#ok<AGROW>
    sum2switch = 0;
    delta_uv=zeros(2,1);
    
    % enforce a minimum angle turn
    while( norm(delta_uv)<5e-2 || abs( delta_uv'*delta_uv_prv ) > 0.3 )
      delta_uv = randn(2,1); 
      delta_uv = delta_uv / norm(delta_uv);
    end
    ax(k) = delta_uv(1);
    ay(k) = delta_uv(2);
    delta_uv_prv = delta_uv + delta_uv_prv;
  end
  
end

ax     = upsample(ax,5);
ay     = upsample(ay,5);

npts   = numel(ax);
integrator_mtrx = gallery('triw',npts,1,npts)';


vx     = integrator_mtrx * ax;
vy     = integrator_mtrx * ay;
xhat   = integrator_mtrx * vx;
yhat   = integrator_mtrx * vy;

% The 'realistic' observed values
H     = integrator_mtrx*integrator_mtrx;
xhatT = H * (ax ); assert( norm(xhatT - xhat ) < 1e-1 );
yhatT = H * (ay ); assert( norm(yhatT - yhat ) < 1e-1 );

idx0  = find( xhatT == 0 );
xhatT(idx0(1:end-1)) = [];
yhatT(idx0(1:end-1)) = [];
npts = numel(xhatT);

h          = fspecial('gaussian',[3 1],1.0);
xhat_noisy = imfilter(xhatT + e_sigma * randn(size(xhatT)),h,'replicate');
yhat_noisy = imfilter(yhatT + e_sigma * randn(size(xhatT)),h,'replicate');
xhat_obs = xhat_noisy;
yhat_obs = yhat_noisy;


sfigure(1); clf; hold on;
plot( xhat, yhat, '-');   hold on;
for k = 1:10:npts
  
  cval = [1.0-k/npts, abs(0.5-k/npts), k/npts];
  plot( xhat_obs(k), yhat_obs(k), '.','color',cval);  
  legend('true path (line)','measurements (points)'); hold on;
  if( ~isempty(find(k==idx_switch, 1) ) )
    plot( xhat(k), yhat(k), '*','color',cval);   hold on;
    plot( xhat_obs(k), yhat_obs(k), '*','color',cval);   hold on;
  end
end
sfigure(1); axis equal; hold off;


xhat0 = smooth(xhat_obs,npts/20);
yhat0 = smooth(yhat_obs,npts/20);

N = npts ;
oversamp = N / npts;
A=eye(N)-diag(ones(N-1,1),1);
B = (A^-1)';
bUseIntegrator = true;
A=A(1:end-1,1:end);
B2 = B*B;
B=B(1:end-1,1:end);
B2=B2(1:end-1,1:end);
dt=12/60 * npts / N;
D=A(1:end-1,1:end-1);
A = sparse(A);
H = zeros(npts,N);
for k = 1:npts
  H(k, oversamp*(k-1) + 1 ) = 1;
end
H = sparse(H);


ell_zero_norm = 100;
ell_zero_max  = 5;
ell_zero_min  = 2;
max_iters     = 10;
iter          = 0;
thresh_sigma = 0.05;

while( (ell_zero_norm > ell_zero_max   || ell_zero_norm < ell_zero_min ) && iter < max_iters)

  iter = iter+1;
  if( ell_zero_norm > ell_zero_max )
    thresh_sigma = 1.2 * thresh_sigma
  else
    thresh_sigma = 0.9 * thresh_sigma
  end

  max_err_y   = thresh_sigma * norm(yhat0)
  max_err_x   = thresh_sigma * norm(xhat0)

  cvx_begin
          variables  ax(N-1) ay(N-1) vx(N) vy(N) Ex(npts) Ey(npts) x(N) y(N)

          minimize(norm(ax+ay,1)+norm(ax-ay,1) )   % group sparse
          %minimize(norm(ay,1)+norm(ax,1))

          subject to
          ones(1,npts)*(Ey.^2) <= max_err_y
          ones(1,npts)*(Ex.^2) <= max_err_x  

          A*vx+dt*(ax) == 0
          A*vy+dt*(ay) == 0
          A*x+dt*vx(1:N-1)+(dt^2)/2*(ax) == 0
          A*y+dt*vy(1:N-1)+(dt^2)/2*(ay) == 0

          y(1)  == yhat0(1);
          y(N)  == yhat0(npts);
          x(1)  == xhat0(1);
          x(N)  == xhat0(npts);

          Ex == xhat0 - H * x
          Ey == yhat0 - H * y
  cvx_end

  ell_zero_norm = sum( abs(ax)+abs(ay) >5e-2 )

end

x_  = H*x;
y_  = H*y;
vx_ = H*vx;
vy_ = H*vy;

sfigure(1); hold on; plot( x_, y_, 'g--','LineWidth',2); hold off;
legend('meas.','recons.'); hold off;
    ylabel('Rotated Y Position [NM]')
            xlabel('Rotated X Position [NM]')
vnorm    = sqrt( vy_.^2 + vx_.^2 )+1e-99;
heading1 = 180/pi * (atan2( vy_(:)./vnorm, vx_(:)./vnorm ));
sfigure(2); clf; hold on; plot( heading1, 'r.','LineWidth',2);       ylabel(sprintf('Heading\n[degrees]'))
hold off







