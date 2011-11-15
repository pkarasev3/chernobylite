function [phi1 phi2 img_show] = run_lskk_demo()
% run demo func in-place:
% [phi1 phi2 img_show] = run_lskk_demo();

addpath('~/source/chernobylite/matlab/util/');
addpath('~/source/chernobylite/matlab/display_helpers/');
addpath('~/source/chernobylite/matlab/LevelSetMethods/');

img   = phantom(); img = img + 0*(randn(size(img))*1e-1 + 0.25).^2; img(img>1)=1;
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
Heavi     = @(z)  1 * (z >= epsilon) + (abs(z) < epsilon).*(1+z/epsilon+1/pi * sin(pi*z/epsilon))/2.0;
delta     = @(z)  1 * (z == 0) + (abs(z) < epsilon).*(1 + cos(pi*z/epsilon))/(epsilon*2.0);

% the cost of overlap that we want to shrink
overlap   = @(p1,p2) trapz(trapz( (Heavi(p1*1e2).*Heavi(p2*1e2)).^2 ) );

% the cost of empty-gap to shrink. note: not symmetric!
gap_pointwise_cost  = @(p1,p2,qE)  ( (Heavi(p1+qE) - Heavi(p1)).*Heavi(p2) ).^2;
emptygap            = @(p1,p2)  trapz(trapz(  gap_pointwise_cost(p1,p2,3*epsilon) ) );

x=linspace(-1,1,100);
sfigure(1); subplot(1,2,2);  plot(x,Heavi(x),'r--'); hold on; plot(x,delta(x),'b-.'); hold off;

lambda     = 0.1;
kappa_phi  = 0*phi1;


tt = 0;
img0 = img;
img  = imfilter( img * 255, fspecial('gaussian',[5 5],2),'replicate');
img  = img - mean(img(:));
delta_rel1 = [1];
delta_rel2 = [1];
relTol     = 1e-3;
phi_show_thresh = 0.9;
while( min([delta_rel1(end),delta_rel2(end)]) > relTol &&  tt < 100 )
 
  
  
  CouplingSymmetric = (Heavi(phi1).*Heavi(phi2));
  Coupling12        = Heavi(phi1).*(Heavi(phi2+3*epsilon)-Heavi(phi2)); 
  C12        = CouplingSymmetric + 1e2*Coupling12;
  prev_phi1  = phi1;
  [phi1 ta]  = update_phi( img, phi1, C12 );
  
  CouplingSymmetric = (Heavi(phi1).*Heavi(phi2));
  Coupling21        = Heavi(phi2).*(Heavi(phi1+3*epsilon)-Heavi(phi1)); 
  C21        = CouplingSymmetric + 1e2*Coupling21;
  prev_phi2  = phi2;
  [phi2 tb]  = update_phi( img, phi2, C21 );
  
  tt         = tt + ta + tb;
  delta_rel1 = [delta_rel1, norm( prev_phi1 - phi1 )/norm(phi1)];
  delta_rel2 = [delta_rel2, norm( prev_phi2 - phi2 )/norm(phi2)];
  
  overlap_cost = overlap(phi1,phi2);
  emptygap_cost= emptygap(phi1,phi2) + emptygap(phi2,phi1);
  sfigure(1); subplot(1,2,1); semilogy( delta_rel1,'r-.' ); hold on; semilogy( delta_rel2,'g--'); hold off;
  legend('\Delta-rel for \phi_1','\Delta-rel for \phi_2'); title('relative levelset function change');
  fprintf('\n overlap cost = %f. empty-gap cost = %f,',overlap_cost,emptygap_cost);
  
  img_show = repmat(img0,[1 1 3]);
  imgr = img_show(:,:,1); imgr( abs( phi1 ) < phi_show_thresh) = (imgr( abs( phi1 ) < phi_show_thresh) + 1.5)/2;
  imgg = img_show(:,:,2); imgg( abs( phi2 ) < phi_show_thresh) = (imgg( abs( phi2 ) < phi_show_thresh) + 1.5)/2;
  img_show(:,:,1) = imgr; img_show(:,:,2) = imgg; img_show(img_show>1)=1; img_show(img_show<0)=0;
  sh=sfigure(1); subplot(1,2,2); imshow(img_show);  title('image with coupled contours');
  setFigure(sh,[10 10],2.25,1.3);
  
  fprintf( 'max-abs-phi = %f, t= %f \n',max(abs(phi1(:))),tt );
  
  drawnow;
  fprintf('');
  
end

  function  [phi dt_a] = update_phi( Img, phi, Coupling )
    
    
    
    mu_i = trapz(trapz(Heavi( phi ) .* Img)) / trapz(trapz(Heavi( phi ) ) );
    mu_o = trapz(trapz( (1-Heavi( phi )) .* Img)) / trapz(trapz( (1-Heavi( phi )) ) );
    
    %kappa_phi(1:numel(phi)) = kappa(phi1,1:numel(phi));
    fprintf('mu_i = %f , mu_o = %f ,',mu_i,mu_o);
      
    g_alpha= (Img - mu_i).^2 - (Img - mu_o).^2 + Coupling;
    dphi  = delta(phi) .* (-g_alpha + lambda * kappa_phi) ;
    
    dt_a  = 1 / max(abs(dphi(:)));
    phi   = phi + dt_a * dphi;
    
    phi =  reinit_SD(phi, 1, 1, 0.8, 'ENO1', 3);
    
  end


end





