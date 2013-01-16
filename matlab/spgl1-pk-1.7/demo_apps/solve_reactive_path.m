cvx_path = genpath('~/source/cvx/');
addpath(cvx_path);
addpath('../../display_helpers/');
addpath('../../util/');

load reactive_data1
%load reactive_data2
xhat0 = xy2(1,:)'; 
yhat0 = xy2(2,:)'; 

preproc = @(v)  smooth(smooth(v,2) + smooth(randn(size(v))*1e-2,1),2);

xhat0=preproc(xhat0); 
yhat0=preproc(yhat0); 


disp( ['measurements l2 error: ' ... 
         num2str([norm(xhat0 - xy2(1,:)'), norm(yhat0 - xy2(2,:)') ]) ] );
meas_err = 1e-1*max(1e-3*var(xhat0), max( [abs(xhat0 - xy2(1,:)') ; abs(yhat0 - xy2(2,:)') ] ));
npts  = numel(xhat0);


N = 1 * npts ;
oversamp = N / npts;
A=eye(N)-diag(ones(N-1,1),1);
A=A(1:end-1,1:end);
D=A(1:end-1,1:end-1);
A = sparse(A); D = sparse(D);
H = zeros(npts,N);
for k = 1:npts
  H(k, oversamp*(k-1) + 1 ) = 1;
end
H = sparse(H);

e_sigma       = 0.1;
ell_zero_norm = numel(xhat0)-1;
ell_zero_max  = 20;
ell_zero_min  = 1;
max_iters     = 30;
iter          = 0;
thresh_sigma  = 0.3;


% differential tolerance
dTol          = 5*( mean( abs( diff( xhat0 ) ) ) + mean( abs( diff( yhat0 ) ) ) );
max_err_xy    = 1e-2 + meas_err;  

nnz_all       = [];
fval_all      = [];
dTol_all      = [];
bRaise_dTol   = false();
erxy_all      = [];
while( (ell_zero_norm > ell_zero_max   || ell_zero_norm < ell_zero_min ) && iter < max_iters)

  iter = iter+1;
  if( ell_zero_norm > ell_zero_max )
      
    dTol         = 1.05 * dTol
  else
    dTol         = dTol * 0.95
   % max_err_xy   = max_err_xy * 0.99
  end
  
  dTol_all       = [dTol_all, dTol]
  
  fprintf('nnz = %f \n',ell_zero_norm);
  
  cvx_begin
          variables  Ax(N-1) Ay(N-1) ax(N-1) ay(N-1) vx(N) vy(N) Ex(npts) Ey(npts) x(N) y(N)

          % group sparse L1/L_infty
          
          minimize( ( norm(ax+ay,1)+norm(ax-ay,1) ) )  
          
          %minimize(norm(ay,1)+norm(ax,1))
          
          subject to
          
          % doh ...
          %-ax(2:N-1)'*eye(N-2,N-2)*ax(1:N-2) <= 0
          %-ay(2:N-1)'*eye(N-2,N-2)*ay(1:N-2) <= 0
          %-ax(3:N-1).*ax(1:N-3) <= 0
          %-ay(3:N-1).*ay(1:N-3) <= 0
          
          ones(1,npts)*(Ex.^2+Ey.^2)           <=     10*max_err_xy
          (Ex.^2+Ey.^2)                        <=       max_err_xy
          
          % ones(1,npts)*(Ex.^2+Ey.^2) <= max_err_x+max_err_y
          
          % Extra 'slack' goes into velocities
          Ax == A*vx+dt*(ax)
          Ay == A*vy+dt*(ay)
          
          % Error in 2nd order taylor series
          0 == A*x+dt*vx(1:N-1)+(dt^2)/2*(ax)
          0 == A*y+dt*vy(1:N-1)+(dt^2)/2*(ay) 
          
          ones(1,npts-1)*(Ax.^2+Ay.^2) <= 2*dTol
          
          
          
%           abs( y(1)  - yhat0(1) ) <= e_sigma*5
%           abs( y(N)  - yhat0(npts) ) <= e_sigma*5
%           abs( x(1)  - xhat0(1) ) <= e_sigma*5
%           abs( x(N)  - xhat0(npts) ) <= e_sigma*5
          
         
          
          Ex == xhat0 - H * x
          Ey == yhat0 - H * y
  cvx_end

  t_input_recons = find( abs(ax) + abs(ay) > 1e-4 );%max([1e-2,0.01 * max(abs(ax)+abs(ay)) ]) );
  ell_zero_norm  = numel(t_input_recons);
  
  C1 = [max( (Ex.^2+Ey.^2)   )  , max_err_xy               ]
  C2 = [ones(1,npts)*(Ex.^2+Ey.^2)   ,    10*max_err_xy]
  C3 = [ones(1,npts-1)*(Ax.^2+Ay.^2) , 2*dTol ]
  
  fval = ( norm(ax+ay,1)+norm(ax-ay,1) )
  
  bRaise_dTol = false();
  if( C3(1) > C3(2) )
    if( (C1(1) < C1(2)) || (C2(1) < C2(2) ) )
      bRaise_dTol   = true();
    end
  else
    
  end
  
  
  
  fprintf('nnz = %f \n',ell_zero_norm);
  nnz_all = [nnz_all, ell_zero_norm]; disp(nnz_all);
  fval_all= [fval_all, fval]; disp(fval_all);
  erxy_all= [erxy_all, max_err_xy];
end

xhat_obs = xhat0;
yhat_obs = yhat0;
ax0 = sparse_ak(1,:)'; %downsample(smooth(sparse_ak(1,:)',20),20);
ay0 = sparse_ak(2,:)'; %downsample(smooth(sparse_ak(2,:)',20),20);
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
idx_nzx = find( abs(ax) > 1e-1 * max( abs(ax) ) );
idx_nzy = find( abs(ay) > 1e-1 * max( abs(ay) ) );
idx_zx  = setdiff(1:(N-1),idx_nzx);
idx_zy  = setdiff(1:(N-1),idx_nzy);
sfigure(1); clf; hold on; 
x_1 = xy1(1,:); y_1 = xy1(2,:);
%plot( x_, y_, 'b-','LineWidth',3); 
plot(xhat_obs, yhat_obs, 'r-','LineWidth',2); 
plot( x_(1), y_(1), 'mo','LineWidth',3,'MarkerSize',8,'MarkerFaceColor',[0 1 0]); 
plot( x_1, y_1, 'b-','LineWidth',2); 
axis equal; 
sh=legend('maneuvering vehicle path','initial position','clutter object path'); 

%sh=legend('sparse-input path estimate','measured data points','initial position'); hold off;

set(sh,'FontSize',14);
ylabel('Y Position [m]','FontSize',16);            xlabel('X Position [m]','FontSize',16);
hold off;
vnorm    = sqrt( vy_.^2 + vx_.^2 )+1e-99;
heading1 = 180/pi * (atan2( vy_(:)./vnorm, vx_(:)./vnorm ));

sfigure(2); clf; hold on; plot( tt,heading1, 'r.','LineWidth',2);
plot(-dt+tt0( t_input_truth ),0*tt0( t_input_truth ),'b^','LineWidth',3,'MarkerSize',12);
axis([0 tt(end) -180 180 ] );
sh=legend('reconstructed vehicle heading','true heading change marker'); set(sh,'FontSize',16); hold off;
ylabel('Heading [degrees]','FontSize',16); 
xlabel('time [sec]','FontSize',16);
hold off

sfigure(3); clf; hold on; 

ax_ = dt*ax; ax_(idx_zx) = 0;
ay_ = dt*ay; ay_(idx_zy) = 0;
plot( [tt],[0;ax_], 'b-.','LineWidth',2);    
plot( [tt],[0;ay_], 'r-','LineWidth',1);
%plot( tt0,ax0 * under_samp_ratio,'c*','LineWidth',1);
%plot( tt0,ay0 * under_samp_ratio,'m*','LineWidth',1);
ylabel('Acceleration [m/sec^2]','FontSize',16); 
xlabel('time [sec]','FontSize',16);
sh=legend('X-acceleration','Y-acceleration'); set(sh,'FontSize',16 );
ax_show = dt*ax(idx_nzx); 
ay_show = dt*ay(idx_nzy); 
plot( tt(idx_nzx+1),ax_show, 'bo','LineWidth',3);    
plot( tt(idx_nzy+1),ay_show, 'rx','LineWidth',2);
axis([0, tt(end), min( -abs([ax_show ; ay_show] ))*1.1, max(  abs([ax_show ; ay_show] ))*1.1 ]);
hold off


dTol_all
nnz_all

set(0,'defaultaxesfontweight','bold');
set(0,'defaultlinelinewidth',2);
set(0,'defaultlinemarkersize',8);
sfigure(5); clf; 
plot( log(fval_all), '-b','Marker','o','MarkerFaceColor',[1 0.5 0]);
grid on; 
xlabel('iterations','FontSize',16); 
ylabel('log f(a_x,a_y)','FontSize',16); 
axis([1,numel(fval_all),min(log(fval_all)-0.5),max(log(fval_all))]);

sfigure(6); clf; 
plot( nnz_all/numel(ax), '--r','Marker','o','MarkerFaceColor',[0 1.0 0.5]);
xlabel('iterations','FontSize',16); 
ylabel('proportion of non-zero (a_x,a_y)','FontSize',16); 
axis([1,numel(fval_all),-0.01,1.01]);

sfigure(7); clf; 
plot( dTol_all, '--g','Marker','s','MarkerFaceColor',[0.4 0.4 0.2],'LineWidth',3,'MarkerSize',8);
xlabel('iterations','FontSize',16); 
ylabel('\delta tolerance on (d_x,d_y)','FontSize',16); 
axis([1,numel(fval_all),0,max(dTol_all)*1.01]);

sfigure(8); clf; 
plot( dist_12(1:end-2), abs(diff(diff(xhat0)))+abs(diff(diff(yhat0))),'r.' )
xlabel('distance [m], clutter to vehicle','FontSize',16);
ylabel('forward-difference |x"|+|y"|, [m/s^2]','FontSize',16);
axis([1 20.5 0 0.3]);

sfigure(9); clf; 
plot(tt(1:end-1),dist_12(1:end-1),'m-');                                                                                    
hold on; 
axP=ax; axP( abs(axP)<abs(circshift(axP,[1,0])) ) = 0; axP( abs(axP)<abs(circshift(axP,[-1,0])) ) = 0;
ayP=ay; ayP( abs(ayP)<abs(circshift(ayP,[1,0])) ) = 0; ayP( abs(ayP)<abs(circshift(ayP,[-1,0])) ) = 0;
plot(tt(abs(axP)+abs(ayP)>2),dist_12(abs(axP)+abs(ayP)>2), ... 
                 'go','MarkerFaceColor',[0.25 0.25 0.5] );  hold off;
xlabel('time (sec)','FontSize',16);
ylabel('distance between clutter and vehicle','FontSize',16);
sh=legend('clutter to vehicle distance','marker: time of detected maneuver');
 set(sh,'FontSize',16 ); axis([0 tt(end) 1 20.0]);
