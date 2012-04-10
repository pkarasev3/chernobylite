function run_U_xi_and_phi_versus_phistar_demo()
set(0,'defaultaxesfontsize',16);
set(0,'defaulttextfontsize',18);
set(0,'defaulttextfontname','Arial');
set(0,'defaultaxesfontweight','bold');
set(0,'defaultlinelinewidth',2);
set(0,'defaultlinemarkersize',4);

dt_init         = 0.7;

%[Dval_alla t_alla] = run_core( sqrt(1/(2)) , dt_init);
%sfigure(2); semilogy(t_alla,Dval_alla,'--','color',[0 0 0.8]); hold on;

%[Dval_allb t_allb] = run_core(     (1/(2)) ,dt_init);
% sfigure(2); semilogy(t_allb,Dval_allb,'-.','color',[0 0.4 .6]); hold on;

[Dval_allc t_allc] = run_core(     1/(2) ,dt_init);
%    sfigure(2); semilogy(t_allc,Dval_allc,'--','color',[0 0.8 .2]); hold on;

%    [Dval_alld t_alld] = run_core(     (1/(16)) ,dt_init);
%    sfigure(2); semilogy(t_alld,Dval_alld,'-.','color',[0.6 0.2 .2]); hold on;
%
%   [Dval_alle t_alle] = run_core(     (1/(256)) ,dt_init);
%    sfigure(2); semilogy(t_alle,Dval_alle,'--','color',[0.9 0.4 .2]);
%
%    legend('\rho=(1/2)^{1/2}','\rho=(1/2)','\rho=(1/4)','\rho=(1/16)','\rho=(1/256)');
%    xlabel('time (sec)');
%    ylabel('labeling error');
%    title('Labeling Error: D(\phi,\phi^*)'); grid on;
%
%    hold off;
%
%    save run_phi_versus_phistar_demo_BAK

end


function [img phi_star phi2] = get_img_phantom_bridge()
img      = phantom();
img(img==0) = 0.1;
[star_i star_j] = find( ( (img < 1e-1) - (img == 0) ) > 1e-3  );
keep_idx  = find( star_j < 128 );
star_i    = star_i(keep_idx);
star_j    = star_j(keep_idx);
idx_left  = sub2ind( size(img), star_i, star_j);
left_stub = img*0-1; left_stub( idx_left ) = 1;
img(145:156,120:145) =  img(150,124); % Make a 'bridge' connecting the two chunks
img = img + (randn(size(img))*5e-2);
img = abs(img+0.1).^(1.5);
img(img>1)=1;

[m n] = size(img);
[xx yy] = meshgrid(1:m,1:n);
xy2     = [100,120]; xy1          = [120,92];
RadInit = 25;%15;
d1      = RadInit - sqrt( ((xx-xy1(1))).^2 + ((yy-xy1(2))).^2 );
d2      = RadInit - sqrt( ((xx-xy2(1))).^2 + ((yy-xy2(2))).^2 );

phi_star  = 1e2*tanh(imfilter( left_stub, fspecial('gaussian',[3 3],1.5),'replicate') );
phi2      = 1e2*tanh(d2);

% phi, "ghost">=1,dX,iters,spatial order, time order
phi_init = phi_star;
phi_star = reinitializeLevelSetFunction(phi_init, 1, 1.0, 500, 3, 2);
phi_init = phi2;
phi2     = reinitializeLevelSetFunction(phi_init, 1, 1.0, 500, 3, 2);

save('initial_data_phantom_bridge.mat','phi_star','phi2','img');

end

function [img phi_star phi2] = get_img_phantom_white()
  img      = phantom();
  img(img==0) = 0.1;
  [star_i star_j] = find( (img>0.295).*(img<3.050) > 0 );
  keep_idx  = find( (star_i <= 114).*(star_i >= 45).*(star_j>=94).*(star_j<=174) > 0 );
  star_i    = star_i(keep_idx);
  star_j    = star_j(keep_idx);
  idx_circ  = sub2ind( size(img), star_i, star_j);
  circ_stub = img*0-1; circ_stub( idx_circ ) = 1;
  
  img = img + (randn(size(img))*5e-2);
  img = abs(img+0.1).^(1.5);
  img(img>1)=1;

  [m n] = size(img);
  [xx yy] = meshgrid(1:m,1:n);
  xy2     = [130,100]; xy1          = [120,92];
  RadInit = 25;%15;
  d1      = RadInit - sqrt( ((xx-xy1(1))).^2 + ((yy-xy1(2))).^2 );
  d2      = RadInit - sqrt( ((xx-xy2(1))).^2 + ((yy-xy2(2))).^2 );

  phi_star  = 1e2*tanh(imfilter( circ_stub, fspecial('gaussian',[3 3],1.5),'replicate') );
  phi2      = 1e2*tanh(d2);

  % phi, "ghost">=1,dX,iters,spatial order, time order
  phi_init = phi_star;
  phi_star = reinitializeLevelSetFunction(phi_init, 1, 1.0, 500, 3, 2);
  phi_init = phi2;
  phi2     = reinitializeLevelSetFunction(phi_init, 1, 1.0, 500, 3, 2);

  save('initial_data_phantom_white.mat','phi_star','phi2','img');

end

function [Dval_all t_all psi1 phi2 img_show U tt xx yy] = run_core( rho_argin, dt_init )
% run demo func in-place:
% [psi1 phi2 img_show] = run_lskk_demo();

dbstop if error;
addpath('~/source/chernobylite/matlab/util/');
addpath('~/source/chernobylite/matlab/display_helpers/');
addpath('~/source/chernobylite/matlab/LevelSetMethods/');
addpath('~/source/chernobylite/matlab/LSMlibPK/');

%[img phi_star phi2] = get_img_phantom_bridge();
[img phi_star phi2] = get_img_phantom_white();


psi1 = phi_star;




sfigure(1); clf;

epsilon   = sqrt(2)*sqrt(2); %1.1;%0.8;
Heavi     = @(z)  1 * (z >= epsilon) + (abs(z) < epsilon).*(1+z/epsilon+1/pi * sin(pi*z/epsilon))/2.0;
delta     = @(z)  (abs(z) < epsilon).*(1 + cos(pi*z/epsilon))/(epsilon*2.0);

% the cost of overlap that we want to shrink
overlap   = @(p1,p2) trapz(trapz( (Heavi(p1*1e2).*Heavi(p2*1e2)).^2 ) );

% the cost of empty-gap to shrink. note: not symmetric!
gap_pointwise_cost  = @(p1,p2,qE)  ( (Heavi(p1+qE) - Heavi(p1)).*(1-Heavi(p2)) ).^2;
emptygap            = @(p1,p2,qE)  trapz(trapz(  gap_pointwise_cost(p1,p2,qE) ) );

x=linspace(-1,1,100);
sfigure(1); subplot(2,1,2);  plot(x,Heavi(x),'r--'); hold on; plot(x,delta(x),'b-.'); hold off;



tt = 0;
img0 = img;
img  = imfilter( img * 255, fspecial('gaussian',[5 5],1),'replicate');
img  = img - min(img(:));
img  = 10*img / max(img(:)) ;

img_show_mid  = img * 0;
img_show_init = img * 0;
phi2_init     = 0 * img;
phi2_mid      = 0 * img;

lambda     = mean(abs(img(:)))^2;
kappa_phi  = 0*psi1;
kappa_psi  = 0*psi1;
kappa_xi   = 0*psi1;
delta_rel1 = [1];
delta_rel2 = [1];
delta_abs1 = [1];
delta_abs2 = [1];

t_all      = [0];
relTol     = 1e-4;
absTol     = 1e2;
phi_show_thresh = max([0.95,epsilon/2.0]);
tsum            = 0;
U               = 0 * psi1;
eps_u           = 1e-1;
steps           = 0;
MaxSteps        = 600;
psi1            = imfilter(phi2,fspecial('gaussian',[5,5],2),'replicate'); % with observer, we start out equal ... no "ground truth"
Dval            = eval_label_dist(psi1,phi2);
Dval_all        = [Dval];
Fval_all        = [0.5];
Norm_U_all      = [0.5];
Norm_du_all     = [0.0];
deltasqr_vol_all= sqrt([trapz(trapz((delta(psi1)).^2))]);

xi_integral    = 0 * phi2;

Gmax            = (min(img(:))-mean(img(:))).^2 + (max(img(:))-mean(img(:))).^2; % maximum the G(\phi,I) term can ever be

dt0             = dt_init;
MaxTime         = 0.25;

rho             =  rho_argin; %(1/2); % 1/16 %(1/4);
alphaPsi        =  1e-1;
Umax            =  sqrt(Gmax/rho + 1)


[~, ~, ~, gval ]  = update_phi( img, psi1, phi2, 0*psi1, 0*psi1, 0);
f_phi=0*phi2;
phi0 = phi2;

imgForU = imfilter(img,ones(3,3)/9,'replicate');

while( (steps < MaxSteps) )
  
  
  
  % Generate and accumulate user inputs
  num_inputs = 5;
  if( steps > 150 )
    num_inputs = 1;
    %lambda          =  (Gmax + 1);
  end
  k = 1;
  U_   = U;
  while(  (steps > 50) && (k <= num_inputs) )  % User is the only place that reference phi_star exists !
    
    idx_u = find( abs( (phi_star > 0).*(  0 > phi2 ) - ...
      (phi_star < 0).*(  0 < phi2 ) ) > 0 );
    %idx_dont_u = find(  (abs(U)<1e-9).*(phi_star>epsilon) > 0 );
    %idx_u = setdiff( idx_u, idx_dont_u );
    
    if( (steps > 400) || (numel(idx_u) < k ) )
      px = 1; py = 1;
    else
      idx_u   = idx_u( randperm(numel(idx_u)) );
      [py px] = ind2sub( size( phi2 ),idx_u(k) );
    end
    U_ = updateU( U_, phi_star,phi2,px,py,imgForU,Umax);
    diffU=norm( U(:)-U_(:) );
    
    if( k<= 1)
      fprintf('diffU = ');
    end
    fprintf(' %6.2g ,  ',diffU);
    k  = k+1;
    if( k==num_inputs)
      fprintf('\n');
    end
    
  end
  
  % Update U
  %U_( U_.*U < 0 ) = U( U_.*U < 0 );
  U               = U_;
  deltaU          = U-U_;
  
  
  
  % Update error signals
  xi                   = Heavi(phi2)-Heavi(psi1);
  eU                   = Heavi(psi1)-Heavi(U);
  
  % Update phi
  f1                   = -(U.^2).*(xi);
  %kappa_xi(:)          = 0*phi2;%
  kappa_xi(:)          = kappa(delta(phi2).^2 .* xi,1:numel(xi(:))); %#ok<NASGU>
  f2                   = lambda * (kappa_xi);
  f_phi                = f1 + f2;     assert( sum(isnan(f_phi(:))) == 0 );
  
  
  redist_iters         = 2;
  phi2_prev            = phi2;
  [psi1 phi2 tb1 dphi ]= update_phi( img, psi1, phi2, 0*psi1, f_phi, redist_iters );
  
  % Update psi
  xi                     = Heavi(phi2)-Heavi(psi1);
  g1                     = xi;
  alphaPsi               = 0.5 / Umax;  % D(t) to zero  F bounded
  %alphaPsi              = 2.0 / Umax;  % F(t) to zero, D bounded
  
  g2                   = -eU .* ( alphaPsi * U).^2;
  f_psi                = g1+g2;
  psi1_prev            = psi1;
  [psi1 phi2 tb2 ~ ]   = update_psi( img, psi1, phi2, f_psi, dphi, redist_iters );
  
  tb = min([tb1 tb2]);
  mf1=max(abs(f_phi(:))); mf2 = max(abs(f_psi(:)));
  fprintf('max f_phi = %f, max f_psi = %f \n',mf1,mf2 );
  
  %phi2 = phi2_prev * ( mf2/(mf1+mf2) ) + phi2 * (mf1/(mf2+mf1));
  %H_U = Heavi(U.^2-epsilon);
  
  % % % % Evaluate whether we're really shrinking D(\phi,\phi^*) % % % %
  Dval        = eval_label_dist(psi1,phi2);
  Fval        = eval_label_dist(psi1,U,alphaPsi*U);
  Fval_all    = [Fval_all, (Fval)];
  Norm_U_all  = [Norm_U_all,  sqrt( trapz(trapz( U.^2 ) ) )];
  Norm_du_all  = [Norm_du_all, sqrt(trapz(trapz(deltaU.^2)))];
  deltasqr_vol= (trapz(trapz( (delta(phi2) ).^2 ) ))+(trapz(trapz( (delta(psi1) ).^2 ) ));
  deltasqr_vol_all = [deltasqr_vol_all, deltasqr_vol];
  Dval_all    = [Dval_all, Dval];                                      %#ok<AGROW>
  
  xi                   = Heavi(phi2)-Heavi(psi1);
  eU                   = Heavi(psi1)-Heavi(U);
  Fhat       = trapz(trapz( delta(psi1).^2 .* (alphaPsi*U).^2 ...
    .* eU.^2 ) );
  FprimeLHS  = trapz(trapz( delta(psi1).^2 .* U.^2 ...
    .* xi .* eU ) );
  FprimeRHS  = trapz(trapz( delta(psi1).^2 .* U.^2 ...
    .* eU.^2 .* (alphaPsi*U).^2 ) );
  fprintf('\nDval = %g ,  Fval = %g ,  Fp_lhs = %g ,  Fp_rhs = %g, Fhat = %g \n',...
    Dval, Fval, FprimeLHS, FprimeRHS, Fhat );
  
  tt         = tt + tb;
  tsum       = tsum + tb;
  t_all      = [t_all, tt]; %#ok<*AGROW>
  
  
  
  % setup display image
  if( steps >= 3 )
    displayLevelSets();
  end
  fprintf('');
  steps = steps+1;
end

result = save_all( );
fprintf('result = %f \n',result);

  function res = save_all( )
    fprintf('done! saving .... \n');
    save run_whole_shebang_demo t_all Fval_all rho_argin...
      phi2_init phi2_mid img_show_mid ...
      img_show_init psi1 phi2 img img_show U ...
      eU Umax alphaPsi tt  steps Dval_all Norm_du_all Norm_U_all
    setenv('rhoval',num2str(rho_argin))
    !cp -v run_whole_shebang_demo.mat  "bridge_demo_rho=${rhoval}_`date +%d%b%Y-%H-%M`.mat"
    res = 1;
  end
  function c = eval_label_dist( phiA, phiB, W )
    
    if( nargin < 3 )
      c = 0.5 * trapz(trapz( (Heavi(phiA)-Heavi(phiB)).^2 ) );
    else
      c = 0.5 * trapz(trapz( W.^2 .* (Heavi(phiA)-Heavi(phiB)).^2 ) );
    end
    
  end

  function  [psi phi dt_a g_source] = update_psi( Img, psi, phi, f_psi, dphi,...
      redist_iters)
    g_source = 0*phi;
    kappa_psi(1:numel(psi)) = kappa(psi,1:numel(psi));
    lambda_now = 0*lambda;
    dpsi  = delta(psi) .* ( f_psi   + lambda_now * kappa_psi) ;
    
    
    both_maxes = max(abs([dpsi(:); dphi(:)])); % max of dphi and dpsi
    dt_a  = dt0 / both_maxes;
    psi    = psi + dt_a * dpsi;
    
    if( redist_iters > 0 )
      dX = 1/sqrt(2);
      psi   =  reinitializeLevelSetFunction(psi,1,dX,redist_iters,3,3,false() );
    end
    fprintf('');
  end

  function  [psi phi dt_a dphidt] = update_phi( Img, psi, phi, f_psi, f_phi,...
      redist_iters)
    mu_i = trapz(trapz(Heavi( phi ) .* Img)) / trapz(trapz(Heavi( phi ) ) );
    mu_o = trapz(trapz( (1-Heavi( phi )) .* Img)) / trapz(trapz( (1-Heavi( phi )) ) );
    
    GofIandPhi = (Img - mu_i).^2 - (Img - mu_o).^2;
    Gmax_now   = max(abs(GofIandPhi(:)));
    assert( Gmax_now <= Gmax );
    g_alpha = -GofIandPhi + f_phi;
    
    lambda_now = 0*lambda;
    g_source   = delta(phi) .*g_alpha;
    dphi       = g_source;
    
    
    both_maxes = [max(abs(dphi(:)))]; % max of dphi and dpsi
    dt_a  = dt0 / max(both_maxes);
    
    phi    = phi + dt_a * dphi;
    
    dphidt= dt_a * dphi;
    
    if( redist_iters > 0 )
      dX = 1/sqrt(2);
      phi   =  reinitializeLevelSetFunction(phi,1,dX,redist_iters,3,3,false() );
      
    end
    fprintf('');
  end

  function displayLevelSets()
    img_show = repmat(img0,[1 1 3]);
    imgb = img_show(:,:,3);
    imgg = img_show(:,:,2);
    imgr = img_show(:,:,1);
    
    % zero out the non-active colors for psi1 (active red), phi2 (active green)
    imgr( abs( phi2 ) < phi_show_thresh ) = 0;
    imgg( abs( psi1 ) < phi_show_thresh ) = 0;
    %imgb( abs( phi2 ) < phi_show_thresh ) = 0;
    %imgb( abs( psi1 ) < phi_show_thresh ) = 0;
    
    imgb( abs(U)>Umax*0.1 ) = (imgb(abs(U)>Umax*0.1)/2 + imgb(abs(U)>Umax*0.1)/Umax );
    
    imgr( abs( psi1 ) < phi_show_thresh) = (imgr( abs( psi1 ) < phi_show_thresh) .* ...
      abs( psi1(abs(psi1) < phi_show_thresh ) )/phi_show_thresh  + ...
      1 * (phi_show_thresh - abs( psi1(abs(psi1) < phi_show_thresh ) ) )/phi_show_thresh );
    
    imgg( abs( phi2 ) < phi_show_thresh) = (imgg( abs( phi2 ) < phi_show_thresh) .* ...
      abs( phi2(abs(phi2) < phi_show_thresh ) )/phi_show_thresh  + ...
      1 * (phi_show_thresh - abs( phi2(abs(phi2) < phi_show_thresh ) ) )/phi_show_thresh );
    
    
    
    %imgr( abs(U)>5 ) = 0; imgg( abs(U)>5 ) = 0;
    %imgb( abs(U)>0 ) = (imgb(abs(U)>0)/2 + abs(U(abs(U)>0))/max(abs(C12(:))) );
    
    img_show(:,:,1) = imgr; img_show(:,:,2) = imgg; img_show(:,:,3) = imgb;
    img_show(img_show>1)=1; img_show(img_show<0)=0;
    sh=sfigure(1,3.5,1.75); subplot(1,2,2); imshow(img_show);
    title(['image and contours, ||U||_2=' num2str(norm(U)) ', t=' num2str_fixed_width(tt,7), ...
      ', steps = ' num2str_fixed_width(steps), ', dt = ' num2str_fixed_width(dt0) ]);
    fprintf( 'max-abs-phi = %f, t= %f, steps = %d \n',max(abs(psi1(:))),tt, steps);
    
    %imwrite(img_show,[ 'closed_loop_bridge_demo_' num2str_fixed_width(steps) '.png']);
    save_title = [ 'data_closed_loop_bridge_demo_' num2str_fixed_width(steps) '.mat'];
    save(save_title,'img','img_show','phi2','phi_star','psi1','U','eU','xi');
    
    sfigure(1); subplot(1,2,1);
    plot( t_all,Dval_all,'r+' ); hold on;
    plot( t_all,Fval_all,'g+');
    hold off;
    legend('D(\phi,\psi)','F(\psi,U)'); % '||\delta(\phi)||_{L2}',
    xlabel('time (sec)');
    title('error signals ');
    sh2=sfigure(2,1.25,2.2);
    subplot(2,1,1);
    
    plot( t_all,Norm_du_all,'b-.'); hold on;
    plot( t_all,Fval_all,'g--');  hold off;
    %legend('||U||_{L2}','F(t)');
    xlabel('time (sec)');
    title('U_t norm and F(t)');
    sfigure(2);  subplot(2,1,2);
    mesh( U ); %axis([1 256 1 256 -mean(U(:)) mean(U(:))]);
    title(sprintf('U(x,t); U_{max} := %g, test max = %g',Umax,max(abs(U(:)))));
    xlabel('x'); ylabel('y');
    
    if( steps == 10 )
      phi2_init = phi2;
      img_show_init = img_show;
      imwrite(img_show,'img_show_init.png');
    elseif( steps == 100 )
      phi2_mid = phi2;
      img_show_mid = img_show;
      imwrite(img_show,'img_show_mid.png');
      fprintf('');
    end
    
    drawnow;
    
    
  end

end

%     semilogy( t_all,delta_rel1,'r-.' ); hold on;
%     semilogy( t_all,delta_rel2,'g--');
%     semilogy( t_all,delta_abs1,'m-.' );
%     semilogy( t_all,delta_abs2,'c--');



