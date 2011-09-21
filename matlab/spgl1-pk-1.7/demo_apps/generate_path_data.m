cvx_path = genpath('~/source/cvx/');
addpath(cvx_path);

npts    = 200;

e_sigma = 100e-1;
u = randn(1,1);
v = randn(1,1);
normuv = sqrt( u^2 + v^2 );
u = u / normuv;
v = v / normuv;
sum2switch = 20;
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
%num_turns = (randn(1,1) > 0 ) * (-1)  + (randn(1,1) > 0 ) * (-1) + 5
num_turns = 3
for k = 2:npts
  
  lambda = 1;
  sum2switch = sum2switch + 1*poissrnd(lambda,1,1);
  
  if( (k < npts*0.95) && sum2switch > npts/num_turns )
    idx_switch = [idx_switch, k-1]; %#ok<AGROW>
    sum2switch = 0;
    delta_uv=zeros(2,1);
    
    % enforce a minimum angle turn
    while( norm(delta_uv)<1e-1 || abs( delta_uv'*delta_uv_prv ) > 0.4 )
      delta_uv = randn(2,1); 
      delta_uv = delta_uv / norm(delta_uv);
    end
    ax(k) = delta_uv(1);
    ay(k) = delta_uv(2);
    delta_uv_prv = delta_uv + delta_uv_prv;
  end
  
end

ax0     = smooth(upsample(ax,5),3);
ay0     = smooth(upsample(ay,5),3);

npts   = numel(ax0);
integrator_mtrx = gallery('triw',npts,1,npts)';


vx0     = integrator_mtrx * ax0;
vy0     = integrator_mtrx * ay0;
xhat   = integrator_mtrx * vx0;
yhat   = integrator_mtrx * vy0;

% The 'realistic' observed values
H     = integrator_mtrx*integrator_mtrx;
xhatT = H * (ax0 ); assert( norm(xhatT - xhat ) < 1e-1 );
yhatT = H * (ay0 ); assert( norm(yhatT - yhat ) < 1e-1 );

idx0  = find( xhatT == 0 ); 
xhatT(idx0) = []; 
yhatT(idx0) = [];
ax0(idx0)   = []; ax0(1:10) = 0;
ay0(idx0)   = []; ay0(1:10) = 0;
vx0(idx0)   = [];
vy0(idx0)   = [];
npts = numel(xhatT);

h          = fspecial('gaussian',[3 1],1.0);
xhat_noisy = imfilter(xhatT + e_sigma * randn(size(xhatT)),h,'replicate');
yhat_noisy = imfilter(yhatT + e_sigma * randn(size(xhatT)),h,'replicate');
xhat_obs = xhat_noisy;
yhat_obs = yhat_noisy;

under_samp_ratio = 8;
xhat_obs = downsample(xhat_obs,under_samp_ratio);
yhat_obs = downsample(yhat_obs,under_samp_ratio);
npts     = numel(xhat_obs);

clear integrator_mtrx

sfigure(1); clf; hold on;
plot( xhat, yhat, '-');   hold on;
for k = 1:floor((npts/20)):npts
  
  cval = [1.0-k/npts, abs(0.5-k/npts), k/npts];
  plot( xhat_obs(k), yhat_obs(k), '.','color',cval);  
  legend('true path (line)','measurements (points)'); hold on;
  if( ~isempty(find(k==idx_switch, 1) ) )
    plot( xhat(k), yhat(k), '*','color',cval);   hold on;
    plot( xhat_obs(k), yhat_obs(k), '*','color',cval);   hold on;
  end
end
sfigure(1); axis equal; hold off;


xhat0 = smooth(xhat_obs,npts/20);  xhat0 = xhat0 - xhat0(1);
yhat0 = smooth(yhat_obs,npts/20);  yhat0 = yhat0 - yhat0(1);

N = 1 * npts ;
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
ell_zero_max  = 8;
ell_zero_min  = 1;
max_iters     = 20;
iter          = 0;
thresh_sigma  = 0.5;

% differential tolerance
dTol          = 0.1*dt * ( mean( abs( diff( xhat0 ) ) ) + mean( abs( diff( yhat0 ) ) ) );

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
          variables  Ax(N-1) Ay(N-1) ax(N-1) ay(N-1) vx(N) vy(N) Ex(npts) Ey(npts) x(N) y(N)

          % group sparse L1/L_infty
          minimize( ( norm(ax+ay,1)+norm(ax-ay,1) ) )  
          
        
          %minimize(norm(ay,1)+norm(ax,1))

          subject to
          ones(1,npts)*(Ey.^2) <= max_err_y
          ones(1,npts)*(Ex.^2) <= max_err_x  

          % Extra 'slack' goes into velocities
          Ax == A*vx+dt*(ax)
          Ay == A*vy+dt*(ay)
          
          % Error in 2nd order taylor series
          0 == A*x+dt*vx(1:N-1)+(dt^2)/2*(ax)
          0 == A*y+dt*vy(1:N-1)+(dt^2)/2*(ay) 
          
          norm(Ay,2) <= dTol
          norm(Ax,2) <= dTol
                    

          abs( y(1)  - yhat0(1) ) <= e_sigma
          abs( y(N)  - yhat0(npts) ) <= e_sigma
          abs( x(1)  - xhat0(1) ) <= e_sigma
          abs( x(N)  - xhat0(npts) ) <= e_sigma
          
          Ex == xhat0 - H * x
          Ey == yhat0 - H * y
  cvx_end

  ell_zero_norm = sum( abs(ax)+abs(ay) >5e-2 )

end

final_differential_error_max = norm(Ay,2) + norm(Ax,2)

x_  = x;
y_  = y;
vx_ = smooth([vx(1);vx(1)+cumsum(ax*dt)]);
vy_ = smooth([vy(1);vy(1)+cumsum(ay*dt)]);

tt   = dt * (1:N);
tt0  = linspace(tt(1),tt(end),numel(ax0));

sfigure(1); clf; hold on; 
plot( xhat0(1:1:end), yhat0(1:1:end), 'r.','LineWidth',3);   hold on;
plot( x_, y_, 'b-','LineWidth',2); axis equal; hold off
legend('measured data points','sparse-input path estimate'); hold off;
ylabel('Y Position [m]','FontSize',14);            xlabel('X Position [m]','FontSize',14);
vnorm    = sqrt( vy_.^2 + vx_.^2 )+1e-99;
heading1 = 180/pi * (atan2( vy_(:)./vnorm, vx_(:)./vnorm ));

sfigure(2); clf; hold on; plot( tt,heading1, 'r.','LineWidth',2);    
axis([0 tt(end) -180 180 ] );
ylabel('Heading [degrees]','FontSize',14); 
xlabel('time [min]','FontSize',14);
hold off

sfigure(3); clf; hold on; 
idx_nz = find( abs(ax)+abs(ay) > 1e-1 );
idx_z = setdiff(1:(N-1),idx_nz);
plot( [tt],[0;ax], 'b--','LineWidth',1);    
plot( [tt],[0;ay], 'r-.','LineWidth',1);
plot( tt0,ax0 * under_samp_ratio,'c*','LineWidth',1);
plot( tt0,ay0 * under_samp_ratio,'m*','LineWidth',1);
ylabel('Acceleration [m/min^2]','FontSize',14); 
xlabel('time [min]','FontSize',14);
legend('X-acceleration','Y-acceleration');
ax_show = 0.5*dt^2 * ax(idx_nz);
ay_show = 0.5*dt^2 * ay(idx_nz);
plot( tt(idx_nz),ax_show, 'bo','LineWidth',3);    
plot( tt(idx_nz),0.5*dt^2 * ay(idx_nz), 'rx','LineWidth',2);
axis([0, tt(end), min( [ax_show ; ay_show] )*1.1, max(  [ax_show ; ay_show] )*1.1 ]);
hold off






