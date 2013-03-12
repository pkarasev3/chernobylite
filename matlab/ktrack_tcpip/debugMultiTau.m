function debugMultiTau()

dbstop if error
load shiftDataT_big;

phi_init = TKR.phi;
phi1     = double(TKR.psi); % transformed phi_init or true mask

% bSwap = true; % Swap the frames because undoing large roll makes more practical sense
% tmp   = phi1; phi1 = phi_init; phi_init=tmp; tmp = TKR.img0; TKR.img0(:)=TKR.img1(:); TKR.img1(:)=tmp(:);
% TKR.g_f2f=(TKR.g_f2f)^-1;

[epsilon,dX] = get_params();

phi_init = reinitializeLevelSetFunction( phi_init, 2, dX,20, 1, 1, true() );

tkr       = TKR;
tkr.g_f2f = (TKR.g_f2f)^-1;
tkr.f     = TKR.f;
m            = tkr.img_size(1);
n            = tkr.img_size(2);
[xx0 yy0]    = meshgrid(linspace(1,n,n),linspace(1,m,m));


% for particle filter, use multiple z0 to sample

%zmin = 40; zmax = 120; zstep = 20;
Tvals = [0:2:10];
showColor = 1;
cumimg    = 0.1*repmat(TKR.img0,[1 1 3]);

wv_ctrl = real(logm(TKR.g_ctrl))
wv_ctrl = wv_ctrl * 0.02 / max(abs(wv_ctrl(:)));
xyF = get_center(phi_init);

nDelays=-2:1:2;
cvals  = { [0,117,220],[255,0,16],[43,206,72],[224,255,102],[220,163,255] };
assert(numel(nDelays)==numel(cvals));
for nd=1:numel(nDelays)
  cvals{nd} = cvals{nd}/255;
  nDelay   = nDelays(nd);
  
  showColor    = 1;%0*(showColor>0) + 1*(showColor<=0);
  z0           = -100;
  xx           = -(xx0 - (n-1)/2) / tkr.f * z0;
  yy           =  (yy0 - (m-1)/2) / tkr.f * z0;
  
  fx  = TKR.fx ;
  fy  = TKR.fy ;
  
  g_ctrl = expm( nDelay * wv_ctrl );
  g_f2fb = TKR.g_f2f * g_ctrl^-1; % g_k = gkd * gkC ; apply inverse of control on the right! 
  g_comp =  (g_f2fb^-1); % *
  uv     = g_comp * [ xx(:)'; yy(:)'; z0 * ones(1,numel(xx)); ones(1,numel(yy)) ];
  
  xx1 =  z0 * uv(1,:)./uv(3,:) ;
  yy1 =  z0 * uv(2,:)./uv(3,:) ;
  xx1  = reshape(xx1,size(xx));
  yy1  = reshape(yy1,size(yy));
  
  phi = interp2(xx,yy,phi_init, xx1, yy1,'*linear',-100);
  phi = reinitializeLevelSetFunction( phi, 2, dX,20, 1, 1, true() );
  
  showLS( 1, TKR.img1, phi_init, cvals{nd});
  cumimg = showLS( 2, cumimg, phi, cvals{nd});
  
  pause(0.025);
  
end

sfigure(2); drawnow; pause(0.01);
%matlab2tikz('ktrack_multi_z0_tikz.tex','relativePngPath','./figs',...
%                                          'width','10.24cm','height','7.68cm');

sfigure(3);
roix=1:32:480; roiy = 1:32:640;
quiver( xx(roix,roiy),yy(roix,roiy), fx(roix,roiy), fy(roix,roiy),2 );
fprintf('');
end

function [xyF] = get_center( phi_in)
  [yc xc] = find( phi_in > 0 ); 
  if isempty(xc) || isempty(yc) 
    fprintf('warning, no phi > 0 region !?\n'); 
    xyF = [tkr.img_size(2)/2 ; tkr.img_size(1)/2];
    return;
  end % strictly speaking mean() is flawed, but OK if blob target like warped circle
  xyF= [mean(xc(:));mean(yc(:))];
  
end

function [img_out] = showLS( nFig, img_in, phi_in,cval)
% better to show as contour ?
if nargin<4
  cval=[0,1,0];
end

[epsilon,dX] = get_params();
sfigure(nFig); clf;
if ndims(img_in) == 2
  img_show = repmat( 0.1*img_in, [1 1 3] );
else
  img_show = img_in; assert(ndims(img_in)==3);
end
imgg= img_show(:,:,2);
phi = phi_in;
imgr= img_show(:,:,1); imgb=img_show(:,:,3); phi_show_thresh = 0.7;%sqrt(2);
imgr( abs( phi ) < phi_show_thresh ) = cval(1);
imgb( abs( phi ) < phi_show_thresh ) = cval(3);
imgg( abs( phi ) < phi_show_thresh) = (imgg( abs( phi ) < phi_show_thresh) .* ...
  abs( phi(abs(phi) < phi_show_thresh ) )/phi_show_thresh  + ...
  1.5 * (phi_show_thresh - abs( phi(abs(phi) < phi_show_thresh ) ) )/phi_show_thresh ) * (0.5 + cval(2));
img_show(:,:,1) = imgr; img_show(:,:,2) = imgg; img_show(:,:,3) = imgb;
img_show(img_show>1)=1; img_show(img_show<0)=0;
imshow((img_show)); %);
drawnow; pause(0.05);
img_out = img_show;
%
% figure(2);
% matlab2tikz('ktrack_ZBuffExample_tikz.tex','relativePngPath','./figs',...
%                                         'width','10.24cm','height','7.68cm');
end

function [epsilon,dX] = get_params()
epsilon   = sqrt(2);
dX        = 0.5 * 1/sqrt(2);
end

function z = Heavi(z)
[epsilon,dX] = get_params();
z = 1 * (z >= epsilon) + (abs(z) < epsilon).*(1+z/epsilon+1/pi * sin(pi*z/epsilon))/2.0;
end

function z = delta(z)
[epsilon,dX] = get_params();
z = (abs(z) < epsilon).*(1 + cos(pi*z/epsilon))/(epsilon*2.0);
end
