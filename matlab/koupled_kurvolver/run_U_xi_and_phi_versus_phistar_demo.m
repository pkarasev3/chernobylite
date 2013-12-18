function run_U_xi_and_phi_versus_phistar_demo( params )
set(0,'defaultaxesfontsize',16);
set(0,'defaulttextfontsize',18);
set(0,'defaulttextfontname','Arial');
set(0,'defaultaxesfontweight','bold');
set(0,'defaultlinelinewidth',2);
set(0,'defaultlinemarkersize',4);

dt_init         = 0.9;

global control_is_on;
global Tord;
global Xord;
global bSaveVerbose;
global Lambda1;
global ImgScale;
Tord = 1; 
Xord = 2;
bSaveVerbose = false();
if nargin==0
  control_is_on = true;
  Lambda1  = 0.1;
  ImgScale = 5.0;
else
  control_is_on = params(1);
  Lambda1       = params(2);
  ImgScale       = params(3);
end

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


function [img phi_star phi2 strtitle] = get_img_phantom_bridge()
  if exist('initial_data_phantom_bridge.mat','file')
    fprintf(' !!!***  loading pre-existing initial data.\nDelete save mat files if you dont want this. \n');
    data = load('initial_data_phantom_bridge.mat');
    img  = data.img; phi_star = data.phi_star; phi2=data.phi2;strtitle=data.strtitle;
  else
    img      = phantom(); 
    img(img==0) = 0.1;
    [star_i star_j] = find( ( (img < 1e-1) - (img == 0) ) > 1e-3  );
    keep_idx  = find( star_j < 128 );
    star_i    = star_i(keep_idx);
    star_j    = star_j(keep_idx);
    idx_left  = sub2ind( size(img), star_i, star_j);
    left_stub = img*0-1; left_stub( idx_left ) = 1;

    before_bridge        = img(142:160,118:147);
    img(142:160,118:147) =  img(150,124); % Make a 'bridge' connecting the two chunks
    [bx by] = meshgrid( linspace(-1,1,size(before_bridge,2)),linspace(-1,1,size(before_bridge,1)));
    bridge_weight = min( 0.5+0*bx, max(max( abs(bx), abs(by) ),0*bx) );
    img(142:160,118:147) = img(142:160,118:147).*(1-bridge_weight) + before_bridge.*(bridge_weight);

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
    strtitle = 'bridge';
    save('initial_data_phantom_bridge.mat','phi_star','phi2','img','strtitle');
    fprintf('saved initial data to file:\n initial_data_phantom_bridge.mat\n');
  end
  
end

function [img phi_star phi2 strtitle] = get_img_phantom_split()

img      = phantom();
img(img==0) = 0.1;
[star_i star_j] = find( ( (img < 1e-1) - (img == 0) ) > 1e-3  );
star_i_all = star_i; star_j_all = star_j;
keep_idx  = find( star_j < 128 );
star_i    = star_i(keep_idx);
star_j    = star_j(keep_idx);
idx_left  = sub2ind( size(img), star_i, star_j);
left_stub = img*0-1; left_stub( idx_left ) = 1;
idx_both  = sub2ind( size(img), star_i_all, star_j_all );
both_stubs= img*0-1; both_stubs( idx_both ) = 1;


img = img + (randn(size(img))*5e-2);
img = abs(img+0.1).^(1.5);
img(img>1)=1;

[m n] = size(img);
[xx yy] = meshgrid(1:m,1:n);

xy3     = [130,130];

RadInit = 45;%15;

d3      = RadInit - sqrt( ((xx-xy3(1))).^2 + ((yy-xy3(2))).^2 );

  phi_star  = 1e2*tanh(imfilter( both_stubs, fspecial('gaussian',[3 3],1.5),'replicate') );
  phi2      = 1e2*tanh(d3);

  phi_init = phi_star;
  phi_star = reinitializeLevelSetFunction(phi_init, 1, 1.0, 500, Xord, Tord);
  phi_init = phi2;
  phi2     = reinitializeLevelSetFunction(phi_init, 1, 1.0, 500, Xord, Tord);

strtitle = 'split';Xord, Tord
  
end

function [img phi_star phi2 strtitle] = get_img_phantom_white()
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
  strtitle = 'white';
end

function [Dval_all t_all psi1 phi2 img_show U tt xx yy] = run_core( rho_argin, dt_init )
% run demo func in-place:
% [psi1 phi2 img_show] = run_lskk_demo();
global control_is_on;
global Tord;
global Xord;
global bSaveVerbose;
global Lambda1;
global ImgScale;
dbstop if error;
addpath('~/source/chernobylite/matlab/util/');
addpath('~/source/chernobylite/matlab/display_helpers/');
addpath('~/source/chernobylite/matlab/LevelSetMethods/');
addpath('~/source/chernobylite/matlab/LSMlibPK/');

[img phi_star phi2 strtitle] = get_img_phantom_bridge(); %get_img_phantom_split(); 
%[img phi_star phi2 strtitle]  = get_img_phantom_white();
img = img * ImgScale;
control_on = control_is_on;
if( ~control_on )
  strtitle = [strtitle '_OpenLoop']; fprintf('***** x_x ');
end
disp(['run config is: '  strtitle]); pause(1);
psi1 = phi_star;




sfigure(1); clf;

epsilon   = sqrt(2)*sqrt(2)*sqrt(2); %1.1;%0.8;
Heavi     = @(z)  1 * (z >= epsilon) + (abs(z) < epsilon).*(1+z/epsilon+1/pi * sin(pi*z/epsilon))/2.0;
delta     = @(z)  (abs(z) < epsilon).*(1 + cos(pi*z/epsilon))/(epsilon*2.0);
deltaPrime= @(z)  ( pi/(2*epsilon^2)).*(-sin(pi*z/epsilon) ) ;

% the cost of overlap that we want to shrink
overlap   = @(p1,p2) trapz(trapz( (Heavi(p1*1e2).*Heavi(p2*1e2)).^2 ) );

% the cost of empty-gap to shrink. note: not symmetric!
gap_pointwise_cost  = @(p1,p2,qE)  ( (Heavi(p1+qE) - Heavi(p1)).*(1-Heavi(p2)) ).^2;
emptygap            = @(p1,p2,qE)  trapz(trapz(  gap_pointwise_cost(p1,p2,qE) ) );

x=linspace(-epsilon*1.25,epsilon*1.25,1000);
sfigure(1); subplot(2,1,2);  
plot(x,Heavi(x),'r--'); hold on; plot(x,delta(x),'b-.'); plot(x,deltaPrime(x),'g--'); hold off;



tt = 0;
img0 = img;
img  = imfilter( img * 255, fspecial('gaussian',[5 5],1),'replicate');
img  = img - min(img(:));
img  = 10*img / max(img(:)) ;
% img    = img / max(img(:));
% img(:) = (histeq(img(:))*5).^2;
% img  = imfilter( img * 10, fspecial('gaussian',[3 3],1),'replicate');
% img  = img - min(img(:));
% img  = 10*img / max(img(:));   img = img - mean(img(:));
% img0 = img; img0 = img0 - min(img0(:)); img0 = img0/max(img0(:));

img_show_mid  = img * 0;
img_show_init = img * 0;
phi2_init     = 0 * img;
phi2_mid      = 0 * img;

lambda     = 0.1*mean(abs(img(:)))^2;
kappa_eu  = 0*psi1;
kappa_xi   = 0*psi1;
delta_rel1 = [1];
delta_rel2 = [1];
delta_abs1 = [1];
delta_abs2 = [1];

t_all      = [0];
relTol     = 1e-4;
absTol     = 1e2;
phi_show_thresh = max([0.95,epsilon]);
tsum            = 0;
U               = 0 * psi1;
eps_u           = 1e-1;
steps           = 0;
MaxSteps        = 500;
psi1            = imfilter(phi2,fspecial('gaussian',[5,5],2),'replicate'); % with observer, we start out equal ... no "ground truth"
Dval            = eval_label_dist(psi1,phi2);
Fval            = eval_label_dist(psi1,U,U);
Dval_all        = [Dval];
Fval_all        = [Fval];
Norm_U_all      = [0.5];
Norm_du_all     = [0.0];
lambda_all      = [0.0];
poincare_all    = [1.0];
mean_i_all      = [0.0];
mean_o_all      = [0.0];
deltasqr_vol_all= sqrt([trapz(trapz((delta(psi1)).^2))]);

xi_integral    = 0 * phi2;

Gmax            = (max(img(:))-min(img(:)))^2; % maximum the G(\phi,I) term can ever be

dt0             = dt_init;
MaxTime         = 0.25;

rho             =  rho_argin; %(1/2); % 1/16 %(1/4);
alphaPsi        =  1e-1;
Umax            =  sqrt(Gmax/rho + 1)


[~, ~, ~, gval ]  = update_phi( img, psi1, phi2, 0*psi1, 0*psi1, 0);
f_phi=0*phi2;
phi0 = phi2;

imgForU = imfilter(img,ones(3,3)/9,'replicate');
G       = 0*phi2;
kappa_phi=0*G;

while( (steps < MaxSteps) )
  
  % Generate and accumulate user inputs
  num_inputs = 3;
  if( steps > 180 )
    num_inputs = 1;
    %lambda          =  (Gmax + 1);
  end
  if( steps > MaxSteps*0.75)
    Tord = 3; Xord = 3;
  end
  k = 1;
  U_   = U;
  stepAtInputStart = 60; % 50 for bridge
  stepAtInputStop  = 300;
  while( control_on && ( steps > stepAtInputStart ) && (k <= num_inputs) )  
    % User is the only place that reference phi_star exists !
    % Simulate their input after some time. 
    
    idx_u = find( abs( (phi_star > 0).*(  0 > phi2 ) - ...
      (phi_star < 0).*(  0 < phi2 ) ) > 0 );
    
    %idx_dont_u = find(  (abs(U)<1e-9).*(phi_star>epsilon) > 0 );
    %idx_u = setdiff( idx_u, idx_dont_u );
    
    if( (steps > stepAtInputStop) || (numel(idx_u) < k ) )
      px = 1; py = 1;
    else
      idx_u   = idx_u( randperm(numel(idx_u)) );
      [py px] = ind2sub( size( phi2 ),idx_u(k) );
    end
    U_ = updateU( U_, phi_star,phi2,px,py,imgForU,Umax);
    diffU=norm( U(:)-U_(:) );
    
    if( k<= 1);      fprintf('diffU = ');    end;
    fprintf(' %6.2g ,  ',diffU);
    k  = k+1;
    if( k==num_inputs);      fprintf('\n');    end;
    
  end
  
  % Update U
  U               = U_;
  deltaU          = U-U_;
  
  % Update error signals
  xi                   = Heavi(phi2)-Heavi(psi1);
  eU                   = Heavi(psi1)-Heavi(U);
  
 % lamNum = (trapz(trapz(delta(phi2).^2 .* G.^2 )))^(1/2);
  %lamDen = (trapz(trapz(delta(phi2).^2 .* abs(xi) )))^(1/2);
  L0     = 1; 
  L1     = 1;
 % lambda = 0*L1 * lamNum  / ( L0 + lamDen );
  alphaDgain = Lambda1 * Gmax; % sweep: .005, 0.01, 0.02, 0.04  (200,100,50,25 Dval limits)
  rvalEst    = 1.1;
  lambda = L0 + rvalEst * alphaDgain * Dval;
  lambda_all = [lambda_all, lambda];
  
  % Update phi
  f1                      = -(U.^2).*(xi);
  % dont need it! [kappa_phi(1:numel(phi2))] = kappa(phi2);
  argKappa = delta(phi2).^2 .* xi;
  [kappa_xi(:) , ~, ~, dx2, dy2 ]          = kappa(argKappa); %#ok<NASGU>
  
  pCgradTerm = trapz(trapz( sqrt(dx2 + dy2) ) );
  pCfuncTerm = trapz(trapz( abs( argKappa ) ) ) ;
  c          = find ( abs(phi2)<=epsilon );
  poincare_r = pCfuncTerm  /  pCgradTerm;
  poincare_all= [poincare_all, poincare_r];
  
  f2                   = lambda * (kappa_xi);
  f_phi                = f1 + f2;     assert( sum(isnan(f_phi(:))) == 0 );
%   if( ~control_on )
%     lambda0 = 2;
%     f_phi = 0*f_phi + lambda0 * kappa_phi;
%     lambda=L1;
%   end
  
  redist_iters          = 2;
  phi2_prev             = phi2;
  [psi1 phi2 tb1 dphi G]= update_phi( img, psi1, phi2, 0*psi1, f_phi, redist_iters );
  Dval                  = eval_label_dist(psi1,phi2);  
  
  
  % Update psi
  psi_iters = 1;
  for mp=1:psi_iters
    xi                     = Heavi(phi2)-Heavi(psi1);
    eU                     = Heavi(psi1)-Heavi(U);
   
    %kappa_psi(1:numel(phi2)) = kappa(phi2,1:numel(phi2));
   % kappa_eu(:)            = kappa(delta(psi1).^2 .* eU,1:numel(eU(:))); %#ok<NASGU>
    
  %  lambda2                = Fval * Lambda2
    g1                     = 1.0 * xi;% + lambda2 * kappa_eu;
    %alphaPsi              = 0.5 / Umax;  % D(t) to zero  F bounded
    alphaPsi               = 1.0 / Umax;  % F(t) to zero, D bounded

    g2                   = -eU .* ( alphaPsi * U).^2;
    f_psi                = g1+g2;
    psi1_prev            = psi1;
    [psi1 phi2 tb2 ~ ]   = update_psi( img, psi1, phi2, f_psi, dphi, redist_iters );
  end
 
  tb = min([tb1 tb2]);
  mf1=max(abs(f_phi(:))); mf2 = max(abs(f_psi(:)));
  fprintf('max f_phi = %f, max f_psi = %f, lambda =%f \n',mf1,mf2,lambda );
  
  %phi2 = phi2_prev * ( mf2/(mf1+mf2) ) + phi2 * (mf1/(mf2+mf1));
  %H_U = Heavi(U.^2-epsilon);
  
  % % % % Evaluate whether we're really shrinking D(\phi,\phi^*) % % % %
  
  Fval        = eval_label_dist(psi1,U,alphaPsi*U);
  Fval_all    = [Fval_all, (Fval)];
  Norm_U_all  = [Norm_U_all,  sqrt( trapz(trapz( U.^2 ) ) )];
  Norm_du_all  = [Norm_du_all, sqrt(trapz(trapz(deltaU.^2)))];
  deltasqr_vol= (trapz(trapz( (delta(phi2) ).^2 ) ))+(trapz(trapz( (delta(psi1) ).^2 ) ));
  deltasqr_vol_all = [deltasqr_vol_all, deltasqr_vol];
  Dval_all    = [Dval_all, Dval];                                      %#ok<AGROW>
  
  xi                   = Heavi(phi2)-Heavi(psi1);
  eU                   = Heavi(psi1)-Heavi(U);
%   Fhat       = trapz(trapz( delta(psi1).^2 .* (alphaPsi*U).^2 ...
%     .* eU.^2 ) );
%   FprimeLHS  = trapz(trapz( delta(psi1).^2 .* U.^2 ...
%     .* xi .* eU ) );
%   FprimeRHS  = trapz(trapz( delta(psi1).^2 .* U.^2 ...
%     .* eU.^2 .* (alphaPsi*U).^2 ) );
  fprintf('\nDval = %g ,  Fval = %g ,  rval = %g \n',...
    Dval, Fval, poincare_r );
  
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
      img_show_init psi1 phi2 lambda_all img img_show U poincare_all ...
      eU Umax alphaPsi tt  steps Dval_all Norm_du_all Norm_U_all
    setenv('Lval',[strtitle '_' num2str(Lambda1)])
    !cp -v run_whole_shebang_demo.mat  "phantom_demo_lambda1=${Lval}_`date +%d%b%Y-%H-%M`.mat"
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
    
    dpsi  = delta(psi) .* ( f_psi  ) ;
    
    
    both_maxes = max(abs([dpsi(:); dphi(:)])); % max of dphi and dpsi
    dt_a  = dt0 / both_maxes;
    psi    = psi + dt_a * dpsi;
    
    if( redist_iters > 0 )
      dX = 1.0;
      psi   =  reinitializeLevelSetFunction(psi,1,dX,redist_iters,Xord,Tord,false() );
    end
    fprintf('');
  end

  function  [psi phi dt_a dphidt G] = update_phi( Img, psi, phi, f_psi, f_phi,...
      redist_iters)
    mu_i = trapz(trapz(Heavi( phi ) .* Img)) / trapz(trapz(Heavi( phi ) ) );
    mu_o = trapz(trapz( (1-Heavi( phi )) .* Img)) / trapz(trapz( (1-Heavi( phi )) ) );
    
    GofIandPhi = (Img - mu_i).^2 - (Img - mu_o).^2;
    Gmax_now   = max(abs(GofIandPhi(:)));
    assert( Gmax_now <= Gmax );
    g_alpha = -GofIandPhi + f_phi;
    G          = -GofIandPhi;
    g_source   = delta(phi) .*g_alpha;
    dphi       = g_source;
    
    
    both_maxes = [max(abs(dphi(:)))]; % max of dphi and dpsi
    dt_a  = dt0 / max(both_maxes);
    
    phi    = phi + dt_a * dphi;
    
    dphidt= dt_a * dphi;
    
    if( redist_iters > 0 )
      dX = 1.0; 
      phi   =  reinitializeLevelSetFunction(phi,1,dX,redist_iters,Xord, Tord,false() );
      
    end
    fprintf('');
  end

  function displayLevelSets()

    phi_show_thresh = 2.0;
    img_show = draw_contour_pair_on_image(img0,psi1,phi2,phi_show_thresh,U,Umax);

    img_show = img_show(17:end-16,65:end-64,:); % truncate ROI  %img_show(9:(end-8),33:(end-32),:);
    
    sh=sfigure(1,3.5,1.75); subplot(1,2,2); imshow(img_show);
    title(['image and contours, ||U||_2=' num2str(norm(U)) ', t=' num2str_fixed_width(tt,7), ...
      ', steps = ' num2str_fixed_width(steps), ', dt = ' num2str_fixed_width(dt0) ]);
    fprintf( 'max-abs-phi = %f, t= %f, steps = %d \n',max(abs(psi1(:))),tt, steps);
    
    
    if bSaveVerbose
    imwrite(img_show,[ strtitle '_imgs_' num2str_fixed_width(steps) '.png']);
    save_title = [ strtitle '_data_' num2str_fixed_width(steps) '.mat'];
    save(save_title,'img','lambda','img_show','phi2','phi_star','psi1','U','eU','xi');
    end
    
    sfigure(1); subplot(1,2,1);
    plot( t_all,Dval_all,'r+' ); hold on;
    plot( t_all,Fval_all,'g+');
    plot( t_all,lambda_all,'m--');
    plot( t_all,poincare_all,'k-s');
    hold off;
    legend('D(\phi,\psi)','F(\psi,U)','\lambda(t)','Int. ratio func/grad'); % '||\delta(\phi)||_{L2}',
    xlabel('time (sec)');
    title('error signals '); grid on;
    sh2=sfigure(2,1.25,2.2);
    subplot(2,1,1);
    
    plot( t_all,Norm_du_all,'b-.'); hold on;
    plot( t_all,Fval_all,'g--');  hold off; grid on;
    %legend('||U||_{L2}','F(t)');
    xlabel('time (sec)');
    title('U_t norm and F(t)');
    if( (steps == MaxSteps) ||  mod(steps,5)==0 )
      sfigure(2);  subplot(2,1,2);
      mesh( U ); %axis([1 256 1 256 -mean(U(:)) mean(U(:))]);
      title(sprintf('U(x,t); U_{max} := %g, test max = %g',Umax,max(abs(U(:)))));
      xlabel('x'); ylabel('y');
    end
    
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



