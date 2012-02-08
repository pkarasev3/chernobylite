function run_U_xi_and_phi_versus_phistar_demo()
  set(0,'defaultaxesfontsize',16);  
  set(0,'defaulttextfontsize',18);
  set(0,'defaulttextfontname','Arial');
  set(0,'defaultaxesfontweight','bold');
  set(0,'defaultlinelinewidth',2);
  set(0,'defaultlinemarkersize',4);
  
  dt_init         = sqrt(2)/2; 
    
  %[Dval_alla t_alla] = run_core( sqrt(1/(2)) , dt_init);
  %sfigure(2); semilogy(t_alla,Dval_alla,'--','color',[0 0 0.8]); hold on;  
  
  %[Dval_allb t_allb] = run_core(     (1/(2)) ,dt_init);
  % sfigure(2); semilogy(t_allb,Dval_allb,'-.','color',[0 0.4 .6]); hold on;  
   
  [Dval_allc t_allc] = run_core(     1/4 ,dt_init);
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



function [Dval_all t_all psi1 phi2 img_show U tt xx yy] = run_core( rho_argin, dt_init )
% run demo func in-place:
% [psi1 phi2 img_show] = run_lskk_demo();

dbstop if error;
addpath('~/source/chernobylite/matlab/util/');
addpath('~/source/chernobylite/matlab/display_helpers/');
addpath('~/source/chernobylite/matlab/LevelSetMethods/');
addpath('~/source/chernobylite/matlab/LSMlibPK/');

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
bUseLSM = true();
if( exist('initial_data_phi_phi-star_img.mat','file' ) ) 
  load initial_data_phi_phi-star_img.mat
  assert( numel(phi2)==numel(img) ); assert( numel(phi_star)==numel(img) ); %#ok<NODEF>
else
  phi_star  = 1e2*tanh(imfilter( left_stub, fspecial('gaussian',[3 3],1.5),'replicate') );
  phi2      = 1e2*tanh(d2);
  if(bUseLSM)
     % phi, "ghost">=1,dX,iters,spatial order, time order
     phi_init = phi_star;
     phi_star = reinitializeLevelSetFunction(phi_init, 1, 1.0, 500, 3, 2);
     phi_init = phi2;
     phi2     = reinitializeLevelSetFunction(phi_init, 1, 1.0, 500, 3, 2);
  else
    phi_star  = reinit_SD(phi_star, 1, 1, 0.5, 'ENO3', 300);
    phi2      = reinit_SD(phi2,1,1,0.5,'ENO3',300);
  end
  
  
end

psi1 = phi_star;

if( ~exist('initial_data_phi_phi-star_img.mat','file' ) )
  save('initial_data_phi_phi-star_img.mat','phi_star','phi2','img');
end


sfigure(1); clf;

epsilon   = sqrt(2); %1.1;%0.8;
Heavi     = @(z)  1 * (z >= epsilon) + (abs(z) < epsilon).*(1+z/epsilon+1/pi * sin(pi*z/epsilon))/2.0;
delta     = @(z)  1 * (z == 0) + (abs(z) < epsilon).*(1 + cos(pi*z/epsilon))/(epsilon*2.0);

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
img  = img - mean(img(:));

img_show_mid  = img * 0;
img_show_init = img * 0;
phi2_init     = 0 * img;
phi2_mid      = 0 * img;

lambda     = mean(abs(img(:)))^2;
kappa_phi  = 0*psi1;
kappa_psi  = 0*psi1;
kappa_eta  = 0*psi1;
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
MaxSteps        = 500;
psi1            = imfilter(phi2,fspecial('gaussian',[5,5],2),'replicate'); % with observer, we start out equal ... no "ground truth"
Dval            = eval_label_dist(psi1,phi2);
Dval_all        = [Dval];
Norm_eU_all     = [0.5];     
Norm_U_all      = [0.5];     
deltasqr_vol_all= sqrt([trapz(trapz((delta(psi1)).^2))]);

eta_integral    = 0 * phi2;

Gmax            = (min(img(:))-mean(img(:))).^2 + (max(img(:))-mean(img(:))).^2; % maximum the G(\phi,I) term can ever be

dt0             = dt_init;
MaxTime         = 0.25;

rho             =  rho_argin; %(1/2); % 1/16 %(1/4);
Umax            =  sqrt(Gmax)/sqrt(rho);
Umax_           =  Umax; %//5.0;

[~, ~, ~, gval ]  = update_phi( img, psi1, phi2, 0*psi1, 0*psi1, 0);
f_phi=0*phi2;
phi0 = phi2;
while( (tt < MaxTime) || (steps < MaxSteps) )
  
  num_inputs = 5;
  k = 1; 
  U_ = U ; %/ Umax * Umax_;
  h_of_u_sum = 0*U;
  
  % eta_integral = eta_integral + 0.25 * ( Heavi(phi2) - Heavi(phi_star) );
     
  while( (steps > 50) && (k < num_inputs) )  % User is the only place that reference phi_star exists ! 
    phi0rightsign = ( (phi_star .* phi0) >= epsilon/2 );
    stucknearphi0bdry = 1+(phi_star.*phi0>=-epsilon).*(phi2.*phi0<=0);
    idx_u = find( abs( (phi_star > 0).*(0 > phi2 ) - ...
                     (phi_star < 0).*(0 < phi2 ) ).* ... 
   max( phi0rightsign , stucknearphi0bdry ) >0 );
    if(numel(idx_u) < k )
      break;
    end
    idx_u   = idx_u( randperm(numel(idx_u)) );
    [py px] = ind2sub( size( phi2 ),idx_u(k) );
    curr_BW = 7;
    
%     if( (steps < 100) && (isempty( intersect(py,142:166) ) || ...
%                           isempty( intersect(px,115:150) ) ) )
%         continue;                
%     end
    
    % add in to h_of_u for better metrics... (img(py,px)-img).^2 .* 
    %h_of_u = Heavi( (curr_BW - ( ( (xx - px).^2 + (yy - py).^2 ) ) ) );
    h_of_u   = exp( -( (xx - px).^2 + (yy - py).^2 )/curr_BW );
    u_in   = (phi_star(py,px) > 0).*(0 > phi2(py,px) ) - ...
            (phi_star(py,px) < 0).*(0 < phi2(py,px) );
    curr_err = abs( phi_star(py,px) - phi2(py,px) );          
%     if ( (curr_err <= 2*epsilon ) )
%               u_in = 0; % they dont click if its doing good
%     end
    
    h_of_u = h_of_u * u_in * Umax * 5e-2;
    k = k+1;
    U_  = U_ + h_of_u;
    
    dtU= 0.05;
    laplacian_of_U = 4*del2(imfilter(U_,fspecial('gaussian',[3 3],0.5),'replicate')); 
    L_ = laplacian_of_U(3:end-2,3:end-2); laplacian_of_U = 0*laplacian_of_U;
    laplacian_of_U(3:end-2,3:end-2) = L_;
    [Ux Uy]        = gradient(U_);
    normGradU      = 0*(Ux.^2+Uy.^2); normGradU(1,:) = 0; normGradU(end,:) = 0; normGradU(:,1) = 0; normGradU(:,end) = 0;
    dU= (laplacian_of_U) .* ( Heavi( (U_.^2 - Umax_.^2) ) ) - U_ .* Heavi( U_.^2 - Umax_.^2).*(1+ normGradU);
    U_ = U_ + dtU * dU;
    
    h_of_u_sum = h_of_u_sum + dtU * dU;
    curr_maxU_ = max(abs(U_(:)));  
     if( abs(u_in) > 0 )
      %fprintf('fixing nominal dynamics at x=%f, y=%f, max(U_)=%f \n ',px, py,curr_maxU_);
    end
  end
  if( 0 ) % show the delta U term 
    sfigure(3); imagesc( img/max(img(:)) + 2*h_of_u_sum/max(1e-9+abs(h_of_u_sum(:))) );  drawnow;
  end 
  
  U   = U_ ;%/ max(U_(:)) * Umax ;% / Umax_ * Umax;
 
    
  % Doesnt hold now assert( sum( U(:) < sqrt( abs(gval(:))/(rho) ) ) == 0 );
  
  prev_psi1            = psi1;
  eta                  = (Heavi(psi1)-Heavi(phi2)); 
  %Del_eta              = 4*del2(eta);
  f1                   = -(U.^2).*(eta);
  %f2                   =  rho*(U.^2).*Del_eta; % rho * (U.^2).*(eta)./( abs(eta).^(3/2)+alph1 );
  kappa_eta(:)         = kappa(eta,1:numel(eta(:)));
  f2                   = rho*Gmax*kappa_eta;
  f_phi                = f1 + f2;
    
  prev_phi2  = phi2;
  redist_iters = round(2*dt0/0.5); redist_iters = min([2,max([redist_iters,1])]);
  
  [psi1 phi2 tb1 gval ]  = update_phi( img, psi1, phi2, 0*psi1, f_phi, redist_iters );
  eta                  = (Heavi(psi1)-Heavi(phi2)); 
  eU                   = (Heavi(psi1)-Heavi(U)).*(U.^2);
  a2                   = Gmax*1e-1;
  a3                   = 1e-1;
  %kappa_eta(:)         = kappa((psi1-phi2),1:numel(eta(:)));
  f_psi                = -( eU )*a3 - eta * a2 + 0*kappa_eta;
  [psi1 phi2 tb2 ~ ]  = update_psi( img, psi1, phi2, f_psi, 0*f_phi, redist_iters );
  tb = max([tb1 tb2]);
  fprintf('max f_phi = %f, max f_psi = %f with Phi part %f \n',...
  max(abs(f_phi(:))),max(abs(f_psi(:))), max(abs(eta(:)*a2)) );  
  
  % % % % Evaluate whether we're really shrinking D(\phi,\phi^*) % % % %
  Dval        = eval_label_dist(psi1,phi2);
  Norm_eUsqrd = (trapz(trapz(  0.5*(Heavi(psi1)-Heavi(U)).^2 .* (abs(U)>1e-9) )));
  Norm_eU_all = [Norm_eU_all, (Norm_eUsqrd)];
  Norm_U_all  = [Norm_U_all,  sqrt( trapz(trapz( U.^2 ) ) )];
  deltasqr_vol= (trapz(trapz( (delta(phi2) ).^2 ) ))+(trapz(trapz( (delta(psi1) ).^2 ) ));
  deltasqr_vol_all = [deltasqr_vol_all, deltasqr_vol];
  Dval_all    = [Dval_all, Dval];                                      %#ok<AGROW>
  fprintf('\nDval = %f\n',Dval);
  
  tt         = tt + tb;
  tsum       = tsum + tb;
  t_all      = [t_all, tt]; %#ok<*AGROW>
  delta_abs1 = [delta_abs1, norm( prev_psi1(:) - psi1(:) )];
  delta_abs2 = [delta_abs2, norm( prev_phi2(:) - phi2(:) )];
  delta_rel1 = [delta_rel1, delta_abs1(end)/norm(psi1(:))];
  delta_rel2 = [delta_rel2, delta_abs2(end)/norm(phi2(:))];
  
  dt0        = dt_init * (1 - exp( -delta_abs1(end)/(sqrt(deltasqr_vol)) ) ) ; 
  dt0        = max([dt0,0.05]);
  
  %overlap_cost = overlap(psi1,phi2);
  %emptygap_cost= emptygap(psi1,phi2,5) + emptygap(phi2,psi1,5);
  % fprintf('\n overlap cost = %f. empty-gap cost = %f,',overlap_cost,emptygap_cost);
  
  % setup display image
  displayLevelSets();
  fprintf('');
  steps = steps+1;
end

result = save_all( );
fprintf('result = %f \n',result);

  function res = save_all( )
    fprintf('done! saving .... \n');
    save run_whole_shebang_demo t_all delta_abs1 delta_abs2 delta_rel1 delta_rel2 rho_argin... 
                             phi2_init phi2_mid img_show_mid img_show_init psi1 phi2 img img_show U tt xx yy  steps Dval_all 
    setenv('rhoval',num2str(rho_argin))
   % !cp -v run_whole_shebang_demo.mat  "bridge_demo_rho=${rhoval}_`date +%d%b%Y-%H-%M`.mat"
    res = 1;                               
  end
  function [Dval] = eval_label_dist( phiA, phiB )                             
    
    Dval = 0.5 * trapz(trapz( (Heavi(phiA)-Heavi(phiB)).^2 ) );
    
  end
  
function  [psi phi dt_a g_source] = update_psi( Img, psi, phi, f_psi, f_phi,... 
                                               redist_iters)
    g_source = 0*phi;   
    kappa_psi(1:numel(psi)) = kappa(psi,1:numel(psi));
    lambda_now = 0*lambda;
    dpsi  = delta(psi) .* ( f_psi   + lambda_now * kappa_psi) ;    
    both_maxes = [max(abs(dpsi(:)))]; % max of dphi and dpsi
    dt_a  = dt0 / max(both_maxes);  
    psi   = psi + dt_a * dpsi ;
    if( redist_iters > 0 )
      dX = 1/sqrt(2);
      if( bUseLSM )
        psi   =  reinitializeLevelSetFunction(psi,1,dX,redist_iters,3,3,false() );
      else
        psi =  reinit_SD(psi, 1.2, dX, dt0, 'ENO3', redist_iters);
      end
    end
    fprintf('');
end

  function  [psi phi dt_a g_source] = update_phi( Img, psi, phi, f_psi, f_phi,... 
                                               redist_iters)
    mu_i = trapz(trapz(Heavi( phi ) .* Img)) / trapz(trapz(Heavi( phi ) ) );
    mu_o = trapz(trapz( (1-Heavi( phi )) .* Img)) / trapz(trapz( (1-Heavi( phi )) ) );
    kappa_phi(1:numel(phi)) = kappa(phi,1:numel(phi));    
    GofIandPhi = (Img - mu_i).^2 - (Img - mu_o).^2;
    Gmax_now   = max(abs(GofIandPhi(:)));
    assert( Gmax_now <= Gmax );
    g_alpha = GofIandPhi + f_phi;
    
    lambda_now = 0*lambda;
    g_source= -delta(phi) .*g_alpha;
    dphi  = delta(phi) .* (-g_alpha +lambda_now * kappa_phi) ;
    fprintf('mu_i = %f, mu_o = %f, g_alpha max = %f, lam*kap max = %f,',...
      mu_i,mu_o,max(abs(g_alpha(:))),max(abs(lambda*kappa_phi(:))));
    
    both_maxes = [max(abs(dphi(:)))]; % max of dphi and dpsi
    dt_a  = dt0 / max(both_maxes);  
    phi   = phi + dt_a * dphi ;
    
    if( redist_iters > 0 )
      dX = 1/sqrt(2);
      if( bUseLSM )
        phi   =  reinitializeLevelSetFunction(phi,1,dX,redist_iters,3,3,false() );
      else
        phi =  reinit_SD(phi, 1.2, dX, dt0, 'ENO3', redist_iters);
      end
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
    %imgb( abs(U)>5 ) = (imgb(abs(U)>5)/2 + abs(U(abs(U)>5))/max(abs(U(:))) );
    
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
    sh=sfigure(1); subplot(1,2,2); imshow(img_show);
    title(['image and contours, ||U||_2=' num2str(norm(U)) ', t=' num2str_fixed_width(tt,7), ... 
             ', steps = ' num2str_fixed_width(steps), ', dt = ' num2str_fixed_width(dt0) ]);
    setFigure(sh,[10 10],3.52,1.76);
    fprintf( 'max-abs-phi = %f, t= %f, steps = %d \n',max(abs(psi1(:))),tt, steps);
    
    %imwrite(img_show,['openloop_bridge_demo_' num2str_fixed_width(steps) '.png']);
        
    
    sfigure(1); subplot(1,2,1);
    semilogy( t_all,Dval_all,'r-.' ); hold on;
    semilogy( t_all,deltasqr_vol_all,'b-' ); 
    semilogy( t_all,Norm_eU_all,'g-.'); 
    hold off;
    legend('D(\phi,\psi)','||\delta(\phi)||_{L2}','||e_U||_{L2}');
    xlabel('time (sec)');
    title('error signals and narrow-band-volume bound');
    sh2=sfigure(2);  
    setFigure(sh2,[1200 10],1.1,2.2);
    subplot(3,1,1);
    
    plot( t_all,Norm_U_all,'b--'); 
    legend('||U||_{L2}'); 
    xlabel('time (sec)');
    title('e_U and U norms');
    sfigure(2);  subplot(2,1,2);
    mesh( U ); %axis([1 256 1 256 -mean(U(:)) mean(U(:))]); 
    title('U(x,t)'); xlabel('x'); ylabel('y');
    
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



