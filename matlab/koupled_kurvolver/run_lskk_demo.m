function [phi1 phi2 img_show U U0 tt xx yy] = run_lskk_demo()
% run demo func in-place:
% [phi1 phi2 img_show] = run_lskk_demo();

addpath('~/source/chernobylite/matlab/util/');
addpath('~/source/chernobylite/matlab/display_helpers/');
addpath('~/source/chernobylite/matlab/LevelSetMethods/');

img   = phantom(); img = img + (randn(size(img))*1e-1); 
img = abs(img+0.1).^(1.5);
img(img>1)=1; 

[m n] = size(img);
[xx yy] = meshgrid(1:m,1:n);

xy2     = [110,112]; xy1          = [120,92];

RadInit = 25;
d1      = RadInit - sqrt( ((xx-xy1(1))).^2 + ((yy-xy1(2))).^2 );
d2      = RadInit - sqrt( ((xx-xy2(1))).^2 + ((yy-xy2(2))).^2 );
phi1    = tanh(3*d1);
phi2    = tanh(3*d2);

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
while( (min([delta_rel1(end),delta_rel2(end)]) > relTol)  || ... 
        (min([delta_abs1(end),delta_abs2(end)]) > absTol) ||  tt < 1 )
  
  % Create instantaneous state change every so often
  bTriggerInput1 = 0;
  if( tsum > 0.01 )
    bTriggerInput1 = 1;
    tsum           = 0;
    Uxy1 = [106,101]; % Input U(x,y,t)
    U0   = 10;
    Rin  = 10;
    dU   = U0*Heavi( Rin - sqrt( ((xx-Uxy1(1))).^2 + ((yy-Uxy1(2))).^2 ) );
    U    = U + dU - (eps_u * U).^3;
    phi1( 0 < (dU < 0).*(phi1>0)  ) = -1;
    phi1( 0 < (dU > 0).*(phi1<=0) ) = +1;
    phi1 =  reinit_SD(phi1, 1, 1, 0.8, 'ENO2', 10);
    
    fprintf('added input at time %f, max U = %f, norm U = %f \n', ...
      tt, max(abs(U(:))), norm(U(:)) );
    
    fprintf('');
    
  end
  
  
  CouplingSymmetric = (Heavi(phi1*1e2).*Heavi(phi2*1e2));
  C12        = CouplingSymmetric + (U.^2).*(-Heavi(U)+Heavi(phi1));
  prev_phi1  = phi1;
  [phi1 ta]  = update_phi( img, phi1, 1e2*C12);
  
  CouplingSymmetric = (Heavi(phi1*1e2).*Heavi(phi2*1e2));
  C21        = CouplingSymmetric;
  prev_phi2  = phi2;
  [phi2 tb]  = update_phi( img, phi2, 1e2*C21 );
  
  tt         = tt + ta + tb;
  tsum       = tsum + ta + tb;
  t_all      = [t_all, tt];
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
  
end

fprintf('done! saving .... \n');
save run_lskk_demo_out phi1 phi2 img_show U U0 tt xx yy

  function  [phi dt_a] = update_phi( Img, phi, Coupling )
    
    
    
    mu_i = trapz(trapz(Heavi( phi ) .* Img)) / trapz(trapz(Heavi( phi ) ) );
    mu_o = trapz(trapz( (1-Heavi( phi )) .* Img)) / trapz(trapz( (1-Heavi( phi )) ) );
    
    kappa_phi(1:numel(phi)) = kappa(phi,1:numel(phi));
    
    g_alpha= (Img - mu_i).^2 - (Img - mu_o).^2 + Coupling;
    dphi  = delta(phi) .* (-g_alpha + lambda * kappa_phi) ;
    
    fprintf('mu_i = %f, mu_o = %f, g_alpha max = %f, lam*kap max = %f,',...
      mu_i,mu_o,max(abs(g_alpha(:))),max(abs(lambda*kappa_phi(:))));
    
    dt0   = 0.8;
    dt_a  = dt0 / max(abs(dphi(:)));  % can go above 1 but then rel-step gets jagged...
    phi   = phi + dt_a * dphi;
    phi =  reinit_SD(phi, 1, 1, dt0, 'ENO2', 2);
    
    
    
  end

  function displayLevelSets()
    img_show = repmat(img0,[1 1 3]);
    imgr = img_show(:,:,1); imgr( abs( phi1 ) < phi_show_thresh) = (imgr( abs( phi1 ) < phi_show_thresh) + 1.5)/2;
    imgg = img_show(:,:,2); imgg( abs( phi2 ) < phi_show_thresh) = (imgg( abs( phi2 ) < phi_show_thresh) + 1.5)/2;
    imgb = img_show(:,:,3);
    imgb( abs(C21)>0 ) = (imgb(abs(C21)>0) + abs(C21(abs(C21)>0))/max(abs(C21(:))) )/2;
    imgb( abs(C12)>0 ) = (imgb(abs(C12)>0) + abs(C12(abs(C12)>0))/max(abs(C12(:))) )/2;
    
    img_show(:,:,1) = imgr; img_show(:,:,2) = imgg; img_show(:,:,3) = imgb;
    img_show(img_show>1)=1; img_show(img_show<0)=0;
    sh=sfigure(1); subplot(2,1,2); imagesc(img_show);  
    title(['image with coupled contours, ||U||_2=' num2str(norm(U)) ', t=' num2str(tt) ]);
    setFigure(sh,[10 10],1.3,3);
    fprintf( 'max-abs-phi = %f, t= %f \n',max(abs(phi1(:))),tt );
    
    sfigure(1); subplot(2,1,1); 
    semilogy( t_all,delta_rel1,'r-.' ); hold on; 
    semilogy( t_all,delta_rel2,'g--'); 
    semilogy( t_all,delta_abs1,'m-.' ); 
    semilogy( t_all,delta_abs2,'c--'); 
    hold off;
    legend('\Delta-rel for \phi_1','\Delta-rel for \phi_2', ...
           '\Delta-abs for \phi_1','\Delta-abs for \phi_2'); 
    title('relative levelset function change');
    
    drawnow;
  
  
  end

end





