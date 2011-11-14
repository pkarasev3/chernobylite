addpath('~/source/chernobylite/matlab/util/');
addpath('~/source/chernobylite/matlab/display_helpers/');

img   = phantom();
[m n] = size(img);
[xx yy] = meshgrid(1:m,1:n);

xy2     = [105,112]; xy1          = [120,92];
phi1    = -1 + 0*img;    phi2    = -1 + 0*img;

RadInit = 15;
d1      = RadInit - sqrt( ((xx-xy1(1))).^2 + ((yy-xy1(2))).^2 );
d2      = RadInit - sqrt( ((xx-xy2(1))).^2 + ((yy-xy2(2))).^2 );

phi1 = tanh(3*d1);
phi2 = tanh(3*d2);

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

lambda     = 0.3;
kappa_phi1 = 0*phi1;

dt = 0.1;
tt = 0;
img0 = img;
img  = img * 255;
img  = img - mean(img(:));
while( tt < 50 )
  tt    = tt + dt;
  
  mu1_i = trapz(trapz(Heavi( phi1 ) .* img)) / trapz(trapz(Heavi( phi1 ) ) );
  mu1_o = trapz(trapz( (1-Heavi( phi1 )) .* img)) / trapz(trapz( (1-Heavi( phi1 )) ) );
  kappa_phi1(1:numel(phi1)) = kappa(phi1,1:numel(phi1));
  fprintf('mu1_i = %f , mu1_o = %f ,',mu1_i,mu1_o);
  
  g_alpha= (img - mu1_o).^2 - (img - mu1_i).^2;
  
  dphi1  = delta(phi1) .* (g_alpha + lambda * kappa_phi1);
  phi1   = phi1 + dt * dphi1; 
  phi1   = imfilter(phi1,fspecial('gaussian',[3 3],0.25),'replicate');
  
  img_show = repmat(img0,[1 1 3]);
  imgr = img_show(:,:,1); imgr( abs( phi1 ) < 0.99) = (imgr( abs( phi1 ) < 0.99) + 1.0)/2;
  imgg = img_show(:,:,2); imgg( abs( phi2 ) < 0.99) = (imgg( abs( phi2 ) < 0.99) + 1.0)/2;
  img_show(:,:,1) = imgr; img_show(:,:,2) = imgg; img_show(img_show>1)=1; img_show(img_show<0)=0;
  sfigure(1); imagesc(img_show); drawnow;
  
  fprintf( 't= %f \n',tt );
  fprintf('');
  
end


img_show = repmat(img0,[1 1 3]);
imgr = img_show(:,:,1); imgr( abs( phi1 ) < 0.5) = 1.0;
imgg = img_show(:,:,2); imgg( abs( phi2 ) < 0.5) = 1.0;
img_show(:,:,1) = imgr; img_show(:,:,2) = imgg; img_show(img_show>1)=1; img_show(img_show<0)=0;
sfigure(1); imagesc(img_show); drawnow;


