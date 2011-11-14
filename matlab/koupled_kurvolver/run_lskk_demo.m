function run_lskk_demo()

addpath('~/source/chernobylite/matlab/util/');
addpath('~/source/chernobylite/matlab/display_helpers/');
addpath('~/source/chernobylite/matlab/LevelSetMethods/');

img   = phantom(); img = img + 0*(randn(size(img))*1e-1 + 0.25).^2; img(img>1)=1;
[m n] = size(img);
[xx yy] = meshgrid(1:m,1:n);

xy2     = [110,112]; xy1          = [120,92];
phi1    = -1 + 0*img;    phi2    = -1 + 0*img;

RadInit = 15;
d1      = RadInit - sqrt( ((xx-xy1(1))).^2 + ((yy-xy1(2))).^2 );
d2      = RadInit - sqrt( ((xx-xy2(1))).^2 + ((yy-xy2(2))).^2 );
phi1    = tanh(3*d1); phi2    = tanh(3*d2);
%phi1 = -reinit_SD(d1, 1, 1, 0.5, 'ENO2', 10);
%phi2 = -reinit_SD(d2, 1, 1, 0.5, 'ENO2', 10);

img_show = repmat(img,[1 1 3]);
imgr = img_show(:,:,1); imgr( abs( phi1 ) < 0.99) = 1.0;
imgg = img_show(:,:,2); imgg( abs( phi2 ) < 0.99) = 1.0;
img_show(:,:,1) = imgr; img_show(:,:,2) = imgg; img_show(img_show>1)=1; img_show(img_show<0)=0;
sfigure(1); clf; imshow(img_show);

epsilon   = 0.9;
Heavi = @(z)  1 * (z >= epsilon) + (abs(z) < epsilon).*(1+z/epsilon+1/pi * sin(pi*z/epsilon))/2.0;
delta = @(z)  1 * (z == 0) + (abs(z) < epsilon).*(1 + cos(pi*z/epsilon))/(epsilon*2.0);

x=linspace(-1,1,100);
sfigure(2); plot(x,Heavi(x),'r--'); hold on; plot(x,delta(x),'b-.'); hold off;

lambda     = 0.1;
kappa_phi  = 0*phi1;
smoother   = fspecial('gaussian',[3 3],0.25);

dt = 0.1;
tt = 0;
img0 = img;
img  = imfilter( img * 255, fspecial('gaussian',[5 5],2),'replicate');
img  = img - mean(img(:));
while( tt < 50 )
  tt    = tt + dt;
  
  
  
  
  Coupling12 = exp( 5 * (Heavi(phi1-phi2) ) );
  phi1       = update_phi( phi1, Coupling12 );

  Coupling21 = exp( 5 * (Heavi(phi1-phi2) ) );
  phi2       = update_phi( phi2, Coupling21 );
  
  img_show = repmat(img0,[1 1 3]);
  imgr = img_show(:,:,1); imgr( abs( phi1 ) < 0.99) = (imgr( abs( phi1 ) < 0.99) + 1.5)/2;
  imgg = img_show(:,:,2); imgg( abs( phi2 ) < 0.99) = (imgg( abs( phi2 ) < 0.99) + 1.5)/2;
  img_show(:,:,1) = imgr; img_show(:,:,2) = imgg; img_show(img_show>1)=1; img_show(img_show<0)=0;
  sfigure(1); imagesc(img_show); drawnow;
  
  fprintf( 'max-abs-phi = %f, t= %f \n',max(abs(phi1(:))),tt );
  fprintf('');
  
end

  function  phi = update_phi( phi, Coupling )
    
    
    
    mu_i = trapz(trapz(Heavi( phi ) .* img)) / trapz(trapz(Heavi( phi ) ) );
    mu_o = trapz(trapz( (1-Heavi( phi )) .* img)) / trapz(trapz( (1-Heavi( phi )) ) );
    kappa_phi(1:numel(phi)) = kappa(phi1,1:numel(phi));
    fprintf('mu_i = %f , mu_o = %f ,',mu_i,mu_o);
      
    g_alpha= (img - mu_o).^2 - (img - mu_i).^2 - Coupling;
    dphi  = delta(phi) .* (g_alpha + lambda * kappa_phi)  ;
    phi   = phi + dt * dphi;
    
    %phi = -reinit_SD(phi, .1, .1, 0.5, 'ENO2', 10);
    phi =  imfilter(phi,smoother,'replicate');
  end


end





