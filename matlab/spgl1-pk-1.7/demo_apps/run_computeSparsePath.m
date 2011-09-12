spgl1_path = genpath('..');
addpath(spgl1_path);
addpath('../../util/');
addpath('../../display_helpers/');

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
else
  load ~/source/visioncontrol/visctrl-papers/input_recovery_adan/trajectory_change_detection/flight_data/IFF_ZMP_20070521_060000_86185_high.mat
  %idx  = ceil( rand(1,1) * length( plane ) );
  idx1 = 10;
  % index 1815 is crazy !
  idx = 3328  %  3328 is moderate
  plane_t = plane(idx);
  [ times_rs, xpRs, ypRs, alts, v_nom] = preprocess_latlon_data( plane_t );
  xhat = [zeros(5,1); xpRs - xpRs(1)];
  yhat = [zeros(5,1); ypRs-ypRs(1)];
  
end


max_val = max( abs( [xhat ; yhat] ) ); xhat = xhat / max_val ; yhat = yhat / max_val;
npts = numel(xhat);
tvals= linspace(0,npts-1,npts)';

Kgroups = npts;
[D b group tau_cen H] = setup_matrices( xhat, yhat, tvals, Kgroups);

opts = spgSetParms();
opts.verbosity  = 2;
opts.iterations = 10000;
opts.iter_skip  = 60;
sval           = 5e-2;
[X,R,G,INFO]   = spg_group(D,b,group,sval,opts);
disp(INFO);

L1_x = X(1:end/2); L1_y = X(end/2+1:end);

sfigure(3); 
plot( L1_x,'b-.'); hold on; plot( L1_y,'r--');

accel_xy = [L1_x(:)'; L1_y(:)'];
x_recons = [0; H * L1_x];
y_recons = [0; H * L1_y];
dy       = smooth( diff([0; y_recons]), 3 );
dx       = smooth( diff([0; x_recons]), 3 );
vnorm    = sqrt( dy.^2 + dx.^2 )+1e-99;
heading0 = 180/pi * (-atan2( dy(:), dx(:) ));
heading  = [dx(:)'; dy(:)' ] ./ [ vnorm' ; vnorm' ];
normal   = [0,1;1,0] * heading;

% TODO: rotate these into the current heading! not "global xy" referenced! 
lin_acc  = sum( heading .* accel_xy, 1 );%.* heading;
tan_acc  = sum( normal  .* accel_xy, 1 );%.* normal;
sfigure(3);  clf; hold on;
plot( tan_acc,'r-.'); plot( lin_acc,'b-'); legend('normal acc','linear acc');
hold off;
sfigure(2); clf; hold on; plot( heading0 ); legend('heading (deg)'); hold off;

%peaks_y = ( imfilter( L1_y,fspecial('log',[5 1],3),'replicate' ) );
%peaks_x = ( imfilter( L1_x,fspecial('log',[5 1],3),'replicate' ) );
%local_max_x = ( abs(peaks_x) > 0.1*max( abs(peaks_x) ) );
%local_max_y = ( abs(peaks_y) > 0.1*max( abs(peaks_y) ) );

[maxtabx, mintabx]=peakdet((smooth(L1_x,11)), 1e-10 );
[maxtaby, mintaby]=peakdet((smooth(L1_y,11)), 1e-10 );
local_max_x = 0*xhat; local_max_x( maxtabx(:,1) ) = 1;
local_max_y = 0*yhat; local_max_y( maxtaby(:,1) ) = 1;
local_min_x = 0*xhat; local_min_x( mintabx(:,1) ) = 1;
local_min_y = 0*yhat; local_min_y( mintaby(:,1) ) = 1;



sfigure(1); clf; hold on;
plot( xhat, yhat, '-');
plot( x_recons, y_recons, 'r--' );
%plot( xhat(local_max_x>0), yhat(local_max_x>0), 'mx','MarkerSize',12);
%%plot( xhat(local_max_y>0), yhat(local_max_y>0), 'cx','MarkerSize',12);
%plot( xhat(local_min_x>0), yhat(local_min_x>0), 'mo','MarkerSize',12);
%plot( xhat(local_min_y>0), yhat(local_min_y>0), 'co','MarkerSize',12);
legend('meas.','recons.'); hold off;

