cvx_path = genpath('~/source/cvx/');
addpath(cvx_path);
addpath('../../display_helpers/');
addpath('../../util/');
load ~/source/visioncontrol/visctrl-papers/input_recovery_adan/trajectory_change_detection/flight_data/IFF_ZMP_20070521_060000_86185_high.mat
  %idx  = ceil( rand(1,1) * length( plane ) );
  %idxf  = 3413; % 3413 is nice and straight, 1 bend
  idx1 = 10;
  % index 1815 is crazy !
  %idxf  = 1815
  idxf = 3328  %  3328 is moderate
  plane_t = plane(idxf);
  [ times_rs, xpRs, ypRs, alts, v_nom] = preprocess_latlon_data( plane_t );
  
  if( idxf == 1815 )
    xpRs = xpRs(400:end);
    ypRs = ypRs(400:end);
  end
  xhat0 = [xpRs - xpRs(1)];
  yhat0 = [ypRs - ypRs(1)];
  
  
  
  vxnom = (xhat0(end)-xhat0(1))/(numel(xhat0)-1);
  vynom = (yhat0(end)-yhat0(1))/(numel(yhat0)-1);
  vy0   = yhat0(2)-yhat0(1);
  xhat_predict = [0; cumsum( vxnom * ones((numel(yhat0)-1),1) ) ] ;
  yhat_predict = [0; cumsum( vynom * ones((numel(yhat0)-1),1) ) ] ;
  xhat  = xhat0 - xhat_predict;
  yhat  = yhat0 - yhat_predict;
  npts = numel(xhat);
  
  
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
ell_zero_max  = 5;
ell_zero_min  = 1;
max_iters     = 20;
iter          = 0;
thresh_sigma  = 0.05 ;

% differential tolerance
dTol          = 0.1*dt * ( mean( abs( diff( xhat0 ) ) ) + mean( abs( diff( yhat0 ) ) ) );
e_sigmax       = 0.05 * mean(abs( xhat0) );
e_sigmay       = 0.05 * mean(abs( yhat0) );

while( (ell_zero_norm > ell_zero_max   || ell_zero_norm < ell_zero_min ) && iter < max_iters)

  iter = iter+1;
  if( ell_zero_norm > ell_zero_max )
    thresh_sigma = 1.1 * thresh_sigma
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
                    

          abs( y(1)  - yhat0(1) ) <= e_sigmay*10
          abs( y(N)  - yhat0(npts) ) <= e_sigmay*10
          abs( x(1)  - xhat0(1) ) <= e_sigmax*10
          abs( x(N)  - xhat0(npts) ) <= e_sigmax*10
          
          Ex == xhat0 - H * x
          Ey == yhat0 - H * y
  cvx_end

  ell_zero_norm = sum( abs(ax)+abs(ay) >5e-2 )

end

final_differential_error_max = norm(Ay,2) + norm(Ax,2)

x_  = x;
y_  = y;
vx_ = medfilt1(smooth([vx(1);vx(1)+cumsum(ax*dt)],1),1);
vy_ = medfilt1(smooth([vy(1);vy(1)+cumsum(ay*dt)],1),1);

tt  = dt * (1:N);
sfigure(1); clf; hold on; 
plot( x_, y_, 'b-','LineWidth',3); %axis equal;
plot(xhat0, yhat0, 'r.','LineWidth',1); 
plot( x_(1), y_(1), 'mo','LineWidth',3,'MarkerSize',8,'MarkerFaceColor',[0 1 0]); plot( x_, y_, 'b-','LineWidth',1);
hold off
sh=legend('sparse-input path estimate','measured data points','initial position'); hold off;
set(sh,'FontSize',16);
ylabel('Y Position [NM]','FontSize',16);            xlabel('X Position [NM]','FontSize',16);
vnorm    = sqrt( vy_.^2 + vx_.^2 )+1e-99;
heading1 = 180/pi * (atan2( vy_(:)./vnorm, vx_(:)./vnorm ));

sfigure(2); clf; hold on; plot( tt,heading1, 'r.','LineWidth',2);
axis([0 tt(end) min(heading1)-10  max(heading1)+10 ] );
sh=legend('reconstructed vehicle heading'); set(sh,'FontSize',16); hold off;
ylabel('Heading [degrees]','FontSize',16); 
xlabel('time [min]','FontSize',16);
hold off

sfigure(3); clf; hold on; 
idx_nz = find( abs(ax)+abs(ay) > 5e-2 * max( abs(ax)+abs(ay) ) );
idx_z = setdiff(1:(N-1),idx_nz);
plot( [tt],[0;ax], 'b-.','LineWidth',2);    
plot( [tt],[0;ay], 'r-','LineWidth',1);

ylabel('Acceleration [NM/min^2]','FontSize',16); 
xlabel('time [min]','FontSize',16);
sh=legend('X-acceleration','Y-acceleration'); set(sh,'FontSize',16 );
ax_show = ax(idx_nz); 
ay_show = ay(idx_nz); 
plot( tt(idx_nz+1),ax_show, 'bo','LineWidth',3);    
plot( tt(idx_nz+1),ay_show, 'rx','LineWidth',2);
axis([0, tt(end), min( -abs([ax_show ; ay_show] ))*1.1, max(  abs([ax_show ; ay_show] ))*1.1 ]);
hold off
