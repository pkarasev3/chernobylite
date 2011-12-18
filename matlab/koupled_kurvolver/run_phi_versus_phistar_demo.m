function [phi1 phi2 img_show U U0 tt xx yy] = run_phi_versus_phistar_demo()
% run demo func in-place:
% [phi1 phi2 img_show] = run_lskk_demo();

dbstop if error;
addpath('~/source/chernobylite/matlab/util/');
addpath('~/source/chernobylite/matlab/display_helpers/');
addpath('~/source/chernobylite/matlab/LevelSetMethods/');

img      = phantom(); 
img(img==0) = 0.1;
[star_i star_j] = find( ( (img < 1e-1) - (img == 0) ) > 1e-3  );
keep_idx  = find( star_j < 128 );
star_i    = star_i(keep_idx);
star_j    = star_j(keep_idx);
idx_left  = sub2ind( size(img), star_i, star_j);
left_stub = img*0-1; left_stub( idx_left ) = 1;
phi_star  = 3*tanh( imfilter( left_stub, fspecial('gaussian',[3 3],1.5),'replicate') );
phi_star  = reinit_SD(phi_star, 1, 1, 0.5, 'ENO2', 100);

img(145:156,120:145) =  img(150,124); % Make a 'bridge' connecting the two chunks

img = img + (randn(size(img))*5e-2); 
img = abs(img+0.1).^(1.5);
img(img>1)=1; 

[m n] = size(img);
[xx yy] = meshgrid(1:m,1:n);

xy2     = [100,120]; xy1          = [120,92];

RadInit = 15;
d1      = RadInit - sqrt( ((xx-xy1(1))).^2 + ((yy-xy1(2))).^2 );
d2      = RadInit - sqrt( ((xx-xy2(1))).^2 + ((yy-xy2(2))).^2 );
phi1    = phi_star;
phi2    = tanh(3*d2);
phi2    = reinit_SD(phi2,1,1,0.5,'ENO2',20);

sfigure(1); clf;

epsilon   = 0.8;
zmax      = 20.0;
Heavi     = @(z)  1 * (z >= epsilon) + (abs(z) < epsilon).*(1+z/epsilon+1/pi * sin(pi*z/epsilon))/2.0;
Heavi2    = @(z)  1 * (z < zmax) .* Heavi(z);
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

lambda     = 2*mean(abs(img(:)))^2;
kappa_phi  = 0*phi1;
delta_rel1 = [1];
delta_rel2 = [1];
delta_abs1 = [1];
delta_abs2 = [1];

t_all      = [0];
relTol     = 1e-4;
absTol     = 1e2;
phi_show_thresh = 0.9;
tsum            = 0;
U               = 0 * phi1;
eps_u           = 1e-1;
steps           = 0;
MaxSteps        = 400;
Dval            = eval_label_dist(phi1,phi2);
Dval_all        = [Dval];

Gmax            = (max(img(:))-min(img(:))).^2; % maximum the G(\phi,I) term can ever be
dt_init         = 0.7;
dt0             = dt_init;


while( (tt < 3) && (steps < MaxSteps) )
  
  prev_phi1            = phi1;
  g_gain               = 1e1;
  [phi2_pred ta gval]  = update_phi( img, phi2, 0 * phi2,0);
  phi1                 = phi_star;
  U                    = max( abs(gval) .* g_gain, U );
  U                    =  2 * sqrt(Gmax) + 0*img;
  eta                  = (Heavi(phi1)-Heavi(phi2));
  f_of_U               = -(U.^2).*(eta);
  
  C21=delta(phi2*0.5).*f_of_U;
  C12=21;
  
  prev_phi2  = phi2;
  [phi2 tb]  = update_phi( img, phi2, f_of_U,2);
  
  % % % % Evaluate whether we're really shrinking D(\phi,\phi^*) % % % %
  Dval_prv = Dval;
  Dval     = eval_label_dist(phi1,phi2);
  Dval_all = [Dval_all, Dval];                                      %#ok<AGROW>
  fprintf('\nDval = %f\n',Dval);
  Dval_pred = eval_label_dist(phi1,phi2_pred);
  bBadDval = false();
  if( Dval > Dval_prv )
    fprintf('Warning, Dval did not decrease, previous value %f !? ',Dval_prv);
    bBadDval = true();
  elseif( Dval_pred < Dval )
    fprintf('Warning, Dval rate worse with f(U), predicted value %f !? ',Dval_pred);
  end
  if( bBadDval ) % If D failed to shrink, reduce the time step
    dt0 = dt0 * 0.5;
    phi_star  = reinit_SD(phi_star, 1, 1, dt0, 'ENO2', 3);
    phi1      = phi_star;
    phi2      = reinit_SD(phi2,     1, 1, dt0, 'ENO2', 3);
    fprintf(', dt0 = %f \n', dt0 );
    fprintf('');
  else           % Raise time-step back up if we're shrinking
    dt0 = dt0 + 0.25 * ( dt_init - dt0 );
  end
  
  tt         = tt + ta + tb;
  tsum       = tsum + ta + tb;
  t_all      = [t_all, tt]; %#ok<*AGROW>
  delta_abs1 = [delta_abs1, norm( prev_phi1(:) - phi1(:) )];
  delta_abs2 = [delta_abs2, norm( prev_phi2(:) - phi2(:) )];
  delta_rel1 = [delta_rel1, delta_abs1(end)/norm(phi1(:))];
  delta_rel2 = [delta_rel2, delta_abs2(end)/norm(phi2(:))];
  
  
  
  overlap_cost = overlap(phi1,phi2);
  emptygap_cost= emptygap(phi1,phi2,5) + emptygap(phi2,phi1,5);
  fprintf('\n overlap cost = %f. empty-gap cost = %f,',overlap_cost,emptygap_cost);
  
  % setup display image
  displayLevelSets();
  fprintf('');
  steps = steps+1;
end

result = save_all( );
fprintf('result = %f \n',result);

  function res = save_all( )
    fprintf('done! saving .... \n');
    save run_openloop_bridge_demo t_all delta_abs1 delta_abs2 delta_rel1 delta_rel2 ... 
                             phi2_init phi2_mid img_show_mid img_show_init phi1 phi2 img img_show U tt xx yy  steps Dval_all 
    res = 1;                               
  end
  function [Dval] = eval_label_dist( phiA, phiB )                             
    
    Dval = trapz(trapz( (Heavi(phiA)-Heavi(phiB)).^2 ) );
    
  end
                             
  function  [phi dt_a g_source] = update_phi( Img, phi, Coupling, ... 
                                               redist_iters)
    
    
    
    mu_i = trapz(trapz(Heavi( phi ) .* Img)) / trapz(trapz(Heavi( phi ) ) );
    mu_o = trapz(trapz( (1-Heavi( phi )) .* Img)) / trapz(trapz( (1-Heavi( phi )) ) );
    
    kappa_phi(1:numel(phi)) = kappa(phi,1:numel(phi));
    
    GofIandPhi = (Img - mu_i).^2 - (Img - mu_o).^2;
    assert( max(abs(GofIandPhi(:)))  <= Gmax );
    g_alpha = GofIandPhi + Coupling;
    
    g_source= -g_alpha + 0 * lambda * kappa_phi ;
    dphi  = delta(phi) .* (-g_alpha + lambda * kappa_phi) ;
    
    
    fprintf('mu_i = %f, mu_o = %f, g_alpha max = %f, lam*kap max = %f,',...
      mu_i,mu_o,max(abs(g_alpha(:))),max(abs(lambda*kappa_phi(:))));
    
    dt_a  = dt0 / max(abs(dphi(:)));  % can go above 1 but then rel-step gets jagged...
    phi   = phi + dt_a * dphi;
    if( redist_iters > 0 )
      phi =  reinit_SD(phi, 1, 1, dt0, 'ENO2', redist_iters);
    end
    
    
    
  end

  function displayLevelSets()
    img_show = repmat(img0,[1 1 3]);
    imgb = img_show(:,:,3);
    imgg = img_show(:,:,2);
    imgr = img_show(:,:,1);
    
    % zero out the non-active colors for phi1 (active red), phi2 (active green)
    imgr( abs( phi2 ) < phi_show_thresh ) = 0;
    imgg( abs( phi1 ) < phi_show_thresh ) = 0;
    imgb( abs( phi2 ) < phi_show_thresh ) = 0;
    imgb( abs( phi1 ) < phi_show_thresh ) = 0;
    
    imgr( abs( phi1 ) < phi_show_thresh) = (imgr( abs( phi1 ) < phi_show_thresh) .* ... 
      abs( phi1(abs(phi1) < phi_show_thresh ) )/phi_show_thresh  + ...
      1.5 * (phi_show_thresh - abs( phi1(abs(phi1) < phi_show_thresh ) ) )/phi_show_thresh );
    
    imgg( abs( phi2 ) < phi_show_thresh) = (imgg( abs( phi2 ) < phi_show_thresh) .* ... 
      abs( phi2(abs(phi2) < phi_show_thresh ) )/phi_show_thresh  + ...
      1.5 * (phi_show_thresh - abs( phi2(abs(phi2) < phi_show_thresh ) ) )/phi_show_thresh );
    
   
    %imgb( abs(C21)>0 ) = (imgb(abs(C21)>0)/2 + abs(C21(abs(C21)>0))/max(abs(C21(:))) );
    %imgb( abs(C12)>0 ) = (imgb(abs(C12)>0)/2 + abs(C12(abs(C12)>0))/max(abs(C12(:))) );
    
    img_show(:,:,1) = imgr; img_show(:,:,2) = imgg; img_show(:,:,3) = imgb;
    img_show(img_show>1)=1; img_show(img_show<0)=0;
    sh=sfigure(1); subplot(1,2,2); imshow(img_show);
    title(['image with coupled contours, ||U||_2=' num2str(norm(U)) ', t=' num2str(tt) ]);
    setFigure(sh,[10 10],3.5,1.9);
    fprintf( 'max-abs-phi = %f, t= %f, steps = %d \n',max(abs(phi1(:))),tt, steps);
    
    %imwrite(img_show,['openloop_bridge_demo_' num2str_fixed_width(steps) '.png']);
        
    
    sfigure(1); subplot(1,2,1);
    semilogy( t_all,Dval_all,'r-.' ); hold on;
    hold off;
    legend('D(\phi,\phi^*');
    title('relative levelset function change');
    
    
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



