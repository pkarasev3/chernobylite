spgl1_path = genpath('..');
addpath(spgl1_path);
addpath('../../util/');
addpath('../../display_helpers/');
cvxpath = genpath('~/source/cvx');
addpath(cvxpath);

%filename = 'pk_to_starbux_loop_raw_sensor_data.txt';
filename = 'starbux_to_pk_raw_sensor_data.txt';
%filename = 'pk_to_thai_raw_sensor_data.txt';
bUsePhoneData = false;
if( bUsePhoneData )
  [xhat yhat tvals] = read_sensor_data(filename);
  idx1  = round(numel(xhat)*0.05);
  xhat  = [xhat( idx1:end)];
  yhat  = [yhat( idx1:end)];
  tvals = tvals(idx1:end);
  xhat  = (imresize(xhat,1/16,'bilinear')); xhat = smooth(xhat-xhat(1),10);
  yhat  = (imresize(yhat,1/16,'bilinear')); yhat = smooth(yhat-yhat(1),10);
  tvals = (imresize(tvals,1/16,'bilinear'));
  xysum = cumsum( abs(yhat)+abs(xhat) );
  idx   = find( xysum > 0 );
  xhat(1:(idx(1)-1)) = [];
  yhat(1:(idx(1)-1)) = [];
  tvals(1:(idx(1)-1))= []; %#ok<NASGU>
  tvals = linspace(0,numel(xhat)-1,numel(xhat))'; %tvals-tvals(1);
  xhat0 = (xhat - xhat(1))*500;
  yhat0 = (yhat - yhat(1))*500;
else
  load ~/source/visioncontrol/visctrl-papers/input_recovery_adan/trajectory_change_detection/flight_data/IFF_ZMP_20070521_060000_86185_high.mat
  %idx  = ceil( rand(1,1) * length( plane ) );
  idxf  = 3413; % 3413 is nice and straight, 1 bend
  idx1 = 10;
  % index 1815 is crazy !
  %idxf  = 1815
  %idxf = 3328  %  3328 is moderate
  plane_t = plane(idxf);
  [ times_rs, xpRs, ypRs, alts, v_nom] = preprocess_latlon_data( plane_t );
  xhat0 = [xpRs - xpRs(1)];
  yhat0 = [ypRs - ypRs(1)];
  
  vxnom = (xhat0(end)-xhat0(1))/(numel(xhat0)-1);
  vynom = (yhat0(end)-yhat0(1))/(numel(yhat0)-1);
  vy0   = yhat0(2)-yhat0(1);
  xhat_predict = [0; cumsum( vxnom * ones((numel(yhat0)-1),1) ) ] ;
  yhat_predict = [0; cumsum( vynom * ones((numel(yhat0)-1),1) ) ] ;
  xhat  = xhat0 - xhat_predict;
  yhat  = yhat0 - yhat_predict;
  
end


%max_val = max( abs( [xhat ; yhat] ) ); xhat = xhat / max_val ; yhat = yhat / max_val;
npts = numel(xhat);


sfigure(1); clf; hold on;
plot( xhat0, yhat0, '-');


N = npts;
A=eye(N)-diag(ones(N-1,1),1);
B = (A^-1)';
bUseIntegrator = true;
A=A(1:end-1,1:end);
B2 = B*B;
B=B(1:end-1,1:end);
B2=B2(1:end-1,1:end);
dt=12/60;
D=A(1:end-1,1:end-1);

max_err_y = 0.15 * norm(yhat0)
max_err_x = 0.15 * norm(xhat0)

cvx_begin
        variables  ax(N-1) ay(N-1) vx(N) vy(N) Ex(N) Ey(N) x(N) y(N)
    
        minimize(norm(ax+ay,1)+norm(ax-ay,1))   % group sparse
        %minimize(norm(ay,1)+norm(ax,1))
        
        subject to
        ones(1,N)*(Ey.^2) <= max_err_y
        ones(1,N)*(Ex.^2) <= max_err_x
        
        A*vx+dt*ax == 0
        A*vy+dt*ay == 0
        A*x+dt*vx(1:N-1)+(dt^2)/2*ax == 0
        A*y+dt*vy(1:N-1)+(dt^2)/2*ay == 0
                       
        y(1)  == yhat0(1);
        y(N)  == yhat0(N);
        x(1)  == xhat0(1);
        x(N)  == xhat0(N);
        
        Ex == xhat0 - x
        Ey == yhat0 - y
cvx_end

sfigure(1); hold on; plot( x, y, 'g--'); hold off;
legend('meas.','recons.'); hold off;
    ylabel('Rotated Y Position [NM]')
            xlabel('Rotated X Position [NM]')
vnorm    = sqrt( vy.^2 + vx.^2 )+1e-99;
heading1 = 180/pi * (-atan2( vy(:)./vnorm, vx(:)./vnorm ));
sfigure(2); clf; hold on; plot( heading1, 'r.','LineWidth',2);       ylabel(sprintf('Heading\n[degrees]'))
%title(['heading for flight # ' num2str( idxf )]);
hold off
