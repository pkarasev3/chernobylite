function [phi1 phi2 img_show U U0 tt xx yy] = run_openloop_chanvese_demo()
% run demo func in-place:
% [phi1 phi2 img_show] = run_lskk_demo();
set(0,'defaultaxesfontsize',16);
set(0,'defaulttextfontsize',18);
set(0,'defaulttextfontname','Arial');
set(0,'defaultaxesfontweight','bold');
set(0,'defaultlinelinewidth',2);
set(0,'defaultlinemarkersize',4);
addpath('../LSMlibPK/');
addpath('../util/');
addpath('../display_helpers/');
addpath('../LevelSetMethods/');

bSaveFinal = false();
bSaveAll   = bSaveFinal && ( false() );

m       = 256;
n       = 256;
[xx yy] = meshgrid(linspace(-1,1,m),linspace(-1,1,n));

img1     = ((0.3 - sqrt( xx.^2 + yy.^2 )) > 0 ) * 1.0;
img2     = 1.0*( xx < 0 );
bSmoothImg2 = false();
if( bSmoothImg2 )
  img2     = imfilter(img2,fspecial('gaussian',[m/8 n/8],m/8/2),'replicate');
end
img      = 0.5 * (img1 + img2);
img      = img - min(img(:)); img = img / max(img(:));


bCauseOverlap = true();
if bCauseOverlap
  xy1     = [.25,0.00];
  xy2     = [.10,-0.20];
else
  xy1     = [0,-0.5];  % center of curve 1
  xy2     = [0.15,0];  % center of curve 2
end

RadInit = 0.21;
d1      = RadInit - sqrt( ((xx-xy1(1))).^2 + ((yy-xy1(2))).^2 );
d2      = RadInit - sqrt( ((xx-xy2(1))).^2 + ((yy-xy2(2))).^2 );
phi1    = sqrt(m*n)*d1;
phi2    = sqrt(m*n)*d2;


epsilon   = sqrt(2);
Heavi     = @(z)  1 * (z >= epsilon) + (abs(z) < epsilon).*(1+z/epsilon+1/pi * sin(pi*z/epsilon))/2.0;
delta     = @(z)  (abs(z) < epsilon).*(1 + cos(pi*z/epsilon))/(epsilon*2.0); %fucking zombie bugs...

% the cost of overlap that we want to shrink
overlap   = @(p1,p2) trapz(trapz( (Heavi(p1*1e2).*Heavi(p2*1e2)).^2 ) );

% the cost of empty-gap to shrink. note: not symmetric!
gap_pointwise_cost  = @(p1,p2,qE)  ( (Heavi(p1+qE) - Heavi(p1)).*(1-Heavi(p2)) ).^2;
emptygap            = @(p1,p2,qE)  trapz(trapz(  gap_pointwise_cost(p1,p2,qE) ) );

x=linspace(-1,1,100);
sfigure(3); clf; plot(x,Heavi(x),'r--'); hold on; plot(x,delta(x),'b-.'); hold off;
sfigure(1); clf; sfigure(2); clf;

phi1_star = (img > 0.75)*2.0 - 1.0; phi1_star=10*phi1_star;
phi2_star = (img > 0.45).*(img < 0.55).*(xx>0)*2.0 - 1.0; phi2_star=10*phi2_star;
sfigure(3); imagesc( [ phi1_star, phi2_star] ); title('[ \phi_1^* , \phi_2^* ]'); axis equal
tt = 0;
img0 = img;
imgForU = imfilter( img, fspecial('gaussian',[3 3],1),'replicate');


lambda     = mean(abs(img(:)));
kappa_phi  = 0*phi1;
delta_rel1 = [1];
delta_rel2 = [1];
delta_abs1 = [1];
delta_abs2 = [1];

[mu_i1, mu_o1] = compute_means(img,phi1);
[mu_i2, mu_o2] = compute_means(img,phi2);
mu1_in_out = [mu_i1;mu_o1];
mu2_in_out = [mu_i2;mu_o2];
t_all      = [0];
relTol     = 1e-4;
absTol     = 1e2;
phi_show_thresh = 0.9;
tsum            = 0;
U1               = 0 * phi1;
U2               = 0 * phi2;
eps_u           = 1e-1;
steps           = 0;
MaxSteps        = 1000;
Umax            = 2.0;

while( ( (min([delta_rel1(end),delta_rel2(end)]) > relTol)  || ...
    (min([delta_abs1(end),delta_abs2(end)]) > absTol) ) &&  (steps < MaxSteps) )
  
  
  steps = steps + 1 ;
  
  
  % Update Phi1
  [U1 deltaU num_inputs] = GenerateUserInput( phi1_star, phi1, U1, imgForU, Umax ); %#ok<ASGLU>
  CouplingSymmetric =  (Heavi(phi1*1e2).*Heavi(phi2*1e2));
  C12        = CouplingSymmetric ;%+ (U.^2).*(-Heavi(U)+Heavi(phi1));
  prev_phi1  = phi1;
  [phi1 ta mui1 muo1]  = update_phi( img, phi1, 1e2*C12);
  
  % Update Phi2
  [U2 deltaU num_inputs] = GenerateUserInput( phi2_star, phi2, U2, imgForU, Umax ); %#ok<ASGLU>
  CouplingSymmetric = (Heavi(phi1*1e2).*Heavi(phi2*1e2));
  C21        = CouplingSymmetric;
  prev_phi2  = phi2;
  [phi2 tb mui2 muo2]  = update_phi( img, phi2, 1e2*C21 );
  
  sfigure(4); clf;imagesc(U1);
  
  mu1_in_out = [mu1_in_out, [mui1;muo1]];
  mu2_in_out = [mu2_in_out, [mui2;muo2]];
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

fprintf('done! .... \n');

if bSaveFinal
  save run_openloop_chanvese_demo t_all delta_abs1 delta_abs2 ...
                 delta_rel1 delta_rel2 phi1 phi2 img img_show U ...
                 U0 tt xx yy  mu1_in_out mu2_in_out steps
end

  

  function  [mu_i mu_o] = compute_means( Img,phi )
     mu_i = trapz(trapz(Heavi( phi ) .* Img)) / trapz(trapz(Heavi( phi ) ) );
     mu_o = trapz(trapz( (1-Heavi( phi )) .* Img)) / trapz(trapz( (1-Heavi( phi )) ) );
  end

  function  [phi dt_a mu_i mu_o] = update_phi( Img, phi, Coupling )
    
    
    [mu_i, mu_o] = compute_means(Img,phi);
   
    
    kappa_phi(1:numel(phi)) = kappa(phi,1:numel(phi));
    
    g_alpha= (Img - mu_i).^2 - (Img - mu_o).^2 + Coupling;
    dphi  = delta(phi) .* (-g_alpha + lambda * kappa_phi) ;
    
    fprintf('mu_i =%2.2f, mu_o =%2.2, g_alpha max =%2.2, lam*kap max =%2.2,',...
      mu_i,mu_o,max(abs(g_alpha(:))),max(abs(lambda*kappa_phi(:))));
    
    dt0   = 0.8;
    dt_a  = dt0 / max(abs(dphi(:)));  % can go above 1 but then rel-step gets jagged...
    phi   = phi + dt_a * dphi;
    phi   =  reinitializeLevelSetFunction(phi,1,1,2,3,2,false() ); %phi =  reinit_SD(phi, 1, 1, dt0, 'ENO2', 2);
    
    
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
    
   
    imgb( abs(C21)>0 ) = (imgb(abs(C21)>0)/2 + abs(C21(abs(C21)>0))/max(abs(C21(:))) );
    imgb( abs(C12)>0 ) = (imgb(abs(C12)>0)/2 + abs(C12(abs(C12)>0))/max(abs(C12(:))) );
    
    img_show(:,:,1) = imgr; img_show(:,:,2) = imgg; img_show(:,:,3) = imgb;
    img_show(img_show>1)=1; img_show(img_show<0)=0;
    sh=sfigure(1); %subplot(1,2,2); 
    imshow(img_show);
    title(['image with coupled contours, ||U_1||_2=' num2str(norm(U1(:))) ', t=' num2str(tt) ]);
    %setFigure(sh,[10 10],3.2,1.5);
    fprintf( 'max-abs-phi = %f, t= %f, steps = %d \n',max(abs(phi1(:))),tt, steps);
    
    if bSaveAll
      imwrite(img_show,['openloop_chanvese_demo_' num2str_fixed_width(steps) '.png']);
    end
    
    sfigure(2); %subplot(1,2,1);
    semilogy( t_all,delta_rel1,'r-.' ); hold on;
    semilogy( t_all,delta_rel2,'g--');
    semilogy( t_all,delta_abs1,'m-.' );
    semilogy( t_all,delta_abs2,'c--');    hold off;
    legend('\Delta-rel for \phi_1','\Delta-rel for \phi_2', ...
      '\Delta-abs for \phi_1','\Delta-abs for \phi_2');
    title('relative levelset function change');
    
    pause(0.001); drawnow; 
    
    
  end

end





