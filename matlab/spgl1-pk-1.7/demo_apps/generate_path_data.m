cvx_path = genpath('~/source/cvx/');
addpath(cvx_path);
addpath('../../display_helpers/');
addpath('../../util/');

nturns_goal = 4;
npts    = 300;

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
num_turns = nturns_goal
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
ell_zero_max  = 6;
ell_zero_min  = 1;
max_iters     = 20;
iter          = 0;
thresh_sigma  = 0.3;

% differential tolerance
dTol          = 0.1*dt * ( mean( abs( diff( xhat0 ) ) ) + mean( abs( diff( yhat0 ) ) ) );

while( (ell_zero_norm > ell_zero_max   || ell_zero_norm < ell_zero_min ) && iter < max_iters)

  iter = iter+1;
  if( ell_zero_norm > ell_zero_max )
    thresh_sigma = 1.1 * thresh_sigma
  else
    thresh_sigma = 0.9 * thresh_sigma
  end

  max_err_y   = thresh_sigma * norm(yhat0)
  max_err_x   = thresh_sigma * norm(xhat0)
  fprintf('nnz = %f \n',ell_zero_norm);
  
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
                    

          abs( y(1)  - yhat0(1) ) <= e_sigma*5
          abs( y(N)  - yhat0(npts) ) <= e_sigma*5
          abs( x(1)  - xhat0(1) ) <= e_sigma*5
          abs( x(N)  - xhat0(npts) ) <= e_sigma*5
          
          Ex == xhat0 - H * x
          Ey == yhat0 - H * y
  cvx_end

  t_input_recons = find( abs(ax) + abs(ay) > 0.01 * max(abs(ax)+abs(ay)) );
  ell_zero_norm  = numel(t_input_recons);
  fprintf('nnz = %f \n',ell_zero_norm);

end

t_input_recons = find( abs(ax) + abs(ay) > 0.01 * max(abs(ax)+abs(ay)) );
t_input_truth  = find( abs(ax0) + abs(ay0) > 0.01 * max(abs(ax0)+abs(ay0)) );
final_differential_error_max = norm(Ay,2) + norm(Ax,2);
ell_zero_norm  = numel(t_input_recons);
fprintf('nnz = %f , ||d_x||_2 + ||d_y||_2 = %f \n',ell_zero_norm,final_differential_error_max);


x_  = x;
y_  = y;
vx_ = medfilt1(smooth([vx(1);vx(1)+cumsum(ax*dt)],5),3);
vy_ = medfilt1(smooth([vy(1);vy(1)+cumsum(ay*dt)],5),3);

tt   = dt * (1:N);
tt0  = linspace(tt(1),tt(end),numel(ax0));

sfigure(1); clf; hold on; 
plot( x_, y_, 'b-','LineWidth',3); %axis equal;
plot(xhat_obs-xhat_obs(1), yhat_obs-yhat_obs(1), 'r.','LineWidth',1); 
plot( x_(1), y_(1), 'mo','LineWidth',3,'MarkerSize',8,'MarkerFaceColor',[0 1 0]); plot( x_, y_, 'b-','LineWidth',1);
hold off
sh=legend('sparse-input path estimate','measured data points','initial position'); hold off;
set(sh,'FontSize',16);
ylabel('Y Position [m]','FontSize',16);            xlabel('X Position [m]','FontSize',16);
vnorm    = sqrt( vy_.^2 + vx_.^2 )+1e-99;
heading1 = 180/pi * (atan2( vy_(:)./vnorm, vx_(:)./vnorm ));

sfigure(2); clf; hold on; plot( tt,heading1, 'r.','LineWidth',2);
plot(-dt+tt0( t_input_truth ),0*tt0( t_input_truth ),'b^','LineWidth',3,'MarkerSize',12);
axis([0 tt(end) -180 180 ] );
sh=legend('reconstructed vehicle heading','true heading change marker'); set(sh,'FontSize',16); hold off;
ylabel('Heading [degrees]','FontSize',16); 
xlabel('time [min]','FontSize',16);
hold off

sfigure(3); clf; hold on; 
idx_nz = find( abs(ax)+abs(ay) > 5e-2 * max( abs(ax)+abs(ay) ) );
idx_z = setdiff(1:(N-1),idx_nz);
plot( [tt],[0;ax], 'b-.','LineWidth',2);    
plot( [tt],[0;ay], 'r-','LineWidth',1);
%plot( tt0,ax0 * under_samp_ratio,'c*','LineWidth',1);
%plot( tt0,ay0 * under_samp_ratio,'m*','LineWidth',1);
ylabel('Acceleration [m/min^2]','FontSize',16); 
xlabel('time [min]','FontSize',16);
sh=legend('X-acceleration','Y-acceleration'); set(sh,'FontSize',16 );
ax_show = ax(idx_nz); 
ay_show = ay(idx_nz); 
plot( tt(idx_nz+1),ax_show, 'bo','LineWidth',3);    
plot( tt(idx_nz+1),ay_show, 'rx','LineWidth',2);
axis([0, tt(end), min( -abs([ax_show ; ay_show] ))*1.1, max(  abs([ax_show ; ay_show] ))*1.1 ]);
hold off






