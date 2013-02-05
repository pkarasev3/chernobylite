function [img_out] = getMetaHorizon( g_WC, img_in, f )
% g_WC:  maps world coords into camera coords
  
  
  dbstop if error
  bTesting = false;
  if nargin == 0
    [g_WC,img_in,K] = getTestArgs();
    bTesting = true;
  else
    K = [ f 0  (size(img_in,2)-1)/2; 
          0 f  (size(img_in,1)-1)/2;
          0 0    1 ];
  end
  G  = (g_WC^(-1));
    
  [u0,v0]=meshgrid(1:size(img_in,2),1:size(img_in,1));
  uv0 = [ u0(:)'; v0(:)'; ones(1,numel(v0))];
  uv1 = [ K^-1 * uv0; ones(1,numel(v0)) ];
 
  % if world coord lies in ground plane somewhere at (x,y,0,1), 
  %       it must have a particular value for z in cam coords
  abcd    = G(3,:);
  z0_proj = -abcd(4)./(-abcd(1)*uv1(1,:)-abcd(2)*uv1(2,:)+abcd(3));

  % below horizon cutoff, record the z-values
  img_out            = zeros( size(img_in,1),size(img_in,2) );
  img_out(z0_proj<0) = abs(z0_proj(z0_proj<0));
  
  
  
  zmax =  500.0;   % cutoff at some max z-value; why not make it be the zFarClip?
  Umax =  500.0;   % set outside range to this
  img_out( img_out > zmax ) = Umax;
  img_out = flipdim(img_out,1);   
  
  if bTesting
    sfigure(1); imagesc(img_out); title('metadata U_{horizon}');
  end
  
  fprintf('');

  
end

function [g_WC,img_in,K] = getTestArgs()

g_WC = [   -0.7071    0.7071         0    0.0002 ;
   -0.1225   -0.1225    0.9849   -9.9491 ;
    0.6964    0.6964    0.1732 -116.6234 ;
         0         0         0    1.0000 ];
       
g_roll = [ [expm( skewsym( 0.5*[1 1 1] ) ) , [0; 0 ;0]]; [0 0 0 1] ];
g_WC = g_WC * g_roll;
       
img_in = zeros(480,640,3);
K = [600 0 319.5; 0 600 239.5; 0 0 1];
end
