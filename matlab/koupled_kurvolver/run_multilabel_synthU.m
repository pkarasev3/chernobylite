function [phi1 phi2 img_show tt xx yy] = run_multilabel_synthU()
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

% to clear GUI prefs:  rmpref('mygraphics')

bSaveFinal = true();%false();
bSaveAll   = bSaveFinal && true%( false() );

m       = 256;
n       = 256;
[xx yy] = meshgrid(linspace(-1,1,m),linspace(-1,1,n));

img1     = ((0.6 - sqrt( xx.^2 + yy.^2 )) > 0 ) * 1.0;
img2     = 1.0*( xx < 0 ).*(1 - 0.15*(yy< 0.25));
bSmoothImg2 = false();
if( bSmoothImg2 )
  img2     = imfilter(img2,fspecial('gaussian',[m/8 n/8],m/8/2),'replicate');
end
img      = 0.5 * (img1 + img2);
img      = img .* (1+randn(size(img))*0.0);
img      = img - min(img(:)); img = img / max(img(:));


bCauseOverlap = true();
if bCauseOverlap
  xy1     = [.35,0.05];
  xy2     = [.15,-0.45];
else
  xy1     = [0,-0.5];  % center of curve 1
  xy2     = [0.15,0];  % center of curve 2
end

RadInit = 0.35;
d1      = RadInit - sqrt( ((xx-xy1(1))).^2 + ((yy-xy1(2))).^2 );
d2      = RadInit - sqrt( ((xx-xy2(1))).^2 + ((yy-xy2(2))).^2 );
phi1    = sqrt(m*n)*d1;
phi2    = sqrt(m*n)*d2;


epsilon   = sqrt(2)*sqrt(2);
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

phi2_star = (img < 0.95).*(img > 0.8)*2.0 - 1.0; phi2_star=10*phi2_star;
phi1_star = (img > 0.45).*(img < 0.55).*(xx>0)*2.0 - 1.0; phi1_star=10*phi1_star;

% Choose: synthetic or from hans johnson
[selectedButton,dlgShown]=uigetpref('mygraphics',... % Group
       'XKsliceMode',...           % Preference
       'Select Mode',...                    % Window title
       {'Select a mode:'},...
       {'synthetic','hjikbrain';'Synthetic Two Contours','Brain Example IKHJ'} );
if strcmp('hjikbrain',selectedButton)     
  [phi1,phi2,phi1_star, phi2_star, img] = load_IKHJ( );
end


sfigure(3); imagesc( [ phi1_star, phi2_star] ); title('[ \phi_1^* , \phi_2^* ]'); axis equal
tt = 0;
img0 = img;
imgForU = imfilter( 5*img, fspecial('gaussian',[3 3],1),'replicate')+1e-3;





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
phi_show_thresh = epsilon;
tsum            = 0;
U1               = 0 * phi1;
U2               = 0 * phi2;
eps_u           = 1e-1;
steps           = 0;
MaxSteps        = 600;
Umax            = 1.5;

Idx_Active      = 0; % 0 => none are active
numInputs       = 10;

prev_xi1 = phi1-phi1_star;
prev_xi2 = phi2-phi2_star;

while( ( (min([delta_rel1(end),delta_rel2(end)]) > relTol)  || ...
     (    min([delta_abs1(end),delta_abs2(end)]) > absTol) ) &&  (steps < MaxSteps) )
  
  steps = steps + 1 ;
  
  if steps > 50
    numInputs = 1+2*(steps<200);
  end
  if steps > MaxSteps/2
    numInputs = 1;
  end
  
  D1 = eval_label_dist(phi1_star,phi1) / trapz(trapz(Heavi(phi1_star))); 
  D2 = eval_label_dist(phi2_star,phi2) / trapz(trapz(Heavi(phi2_star))); 
  
  if steps>3
    if D1 > D2
      Idx_Active = 1;
    else
      Idx_Active = 2;
    end
    fprintf('D1=%3.3f,D2=%3.3f,IdxActive=%01d\n',D1,D2,Idx_Active);
  end
  
  % Update Phi1
  C12        = 0*phi2;%Heavi(1e2*phi2); 
  prev_phi1  = phi1;
  if Idx_Active == 1  
    [U1 deltaU num_inputs] = GenerateUserInput( phi1_star, phi1, prev_xi1, U1, imgForU, ...
                                                       Umax, numInputs ); %#ok<ASGLU>                                                     
    [phi1 ta mui1 muo1]  = update_active_phi( img, phi1, phi2, U1 );
  elseif Idx_Active == 2
    [phi1 ta mui1 muo1]  = update_not_active_phi( img, phi1, phi2, U1 );
  else
    [phi1 ta mui1 muo1]  = update_phi( img, phi1, C12);
  end
  prev_xi1 = (prev_phi1)-(phi1_star);  
  
  if sum( isnan(U1(:)) ) > 0 
    fprintf('4To 3a xep ?\n')
    keyboard;
  end
  
  
  % Update Phi2
  C21        = 0*phi1;%Heavi(1e2*phi1);
  prev_phi2  = phi2;
  if Idx_Active == 2
    [U2 deltaU num_inputs] = GenerateUserInput( phi2_star, phi2, prev_xi2, U2, imgForU, Umax, numInputs ); %#ok<ASGLU>
    [phi2 tb mui2 muo2] = update_active_phi( img, phi2, phi1, U2 );
  elseif Idx_Active == 1
    [phi2 tb mui2 muo2] = update_not_active_phi( img, phi2, phi1, U2 );
  else
    [phi2 tb mui2 muo2]  = update_phi( img, phi2, C21 );
  end
  prev_xi2 =(prev_phi2)-(phi2_star);

  
  sfigure(4); subplot(2,1,1); imagesc(U1); axis image; subplot(2,1,2); imagesc(U2); axis image;
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
                 delta_rel1 delta_rel2 phi1 phi2 img img_show U1 ...
                 U2 tt xx yy  mu1_in_out mu2_in_out steps
end

  function c = eval_label_dist( phiA, phiB, W )
    if( nargin < 3 )
      c = 0.5 * trapz(trapz( (Heavi(phiA)-Heavi(phiB)).^2 ) );
    else
      c = 0.5 * trapz(trapz( W.^2 .* (Heavi(phiA)-Heavi(phiB)).^2 ) );
    end
  end

  function  [mu_i mu_o] = compute_means( Img,phi )
     mu_i = trapz(trapz(Heavi( phi ) .* Img)) / trapz(trapz(Heavi( phi ) ) );
     mu_o = trapz(trapz( (1-Heavi( phi )) .* Img)) / trapz(trapz( (1-Heavi( phi )) ) );
  end

  function  [mu_i mu_o] = compute_local_means( Img,phi )
     HP   = Heavi(phi).*Heavi(abs(phi)+5.0);
     HM   =(1- Heavi(phi)).*Heavi(abs(phi)+5.0);
     mu_i = trapz(trapz(HP.* Img)) / trapz(trapz(HP) );
     mu_o = trapz(trapz( HM .* Img)) / trapz(trapz( HM ) );
  end

  function  [phi dt_a mu_i mu_o] = update_phi( Img, phi, Coupling )
    [mu_i, mu_o] = compute_means(Img,phi);
    kappa_phi    = kappa(phi);
    g_alpha      = (Img - mu_i).^2 - (Img - mu_o).^2;
    dphi         = delta(phi) .* (-g_alpha + lambda * kappa_phi + Coupling) ;
    fprintf('NoneActive, mu_i =%2.2f, mu_o =%2.2, g_alpha max =%2.2, lam*kap max =%2.2,',...
      mu_i,mu_o,max(abs(g_alpha(:))),max(abs(lambda*kappa_phi(:))));
    
    dt0   = 0.8;
    dt_a  = dt0 / max(abs(dphi(:)));  
    phi   = phi + dt_a * dphi;
    phi   =  reinitializeLevelSetFunction(phi,1,1,2,3,2,true() ); %phi =  reinit_SD(phi, 1, 1, dt0, 'ENO2', 2);
    
  end

  function  [phi dt_a mu_i mu_o] = update_active_phi( Img, phi, phiInactive, U)
    %[mu_i, mu_o] = compute_means(Img,phi);
    [mu_i, mu_o] = compute_local_means(Img,phi);
    kappa_phi(:) = kappa(phi,1:numel(phi));
    g_alpha      = (Img - mu_i).^2 - (Img - mu_o).^2 ;
    dphi         =  (-g_alpha + lambda * kappa_phi) ;
    
    gammaSqr = 0.1;
    dcpl     = ( Heavi(phiInactive) ).*(-gammaSqr + g_alpha - lambda*kappa_phi );
    dUin     = abs(U).* (Heavi(U) - Heavi(phi));
    dphi     = delta(phi) .* (dphi + dcpl + dUin);
    
    fprintf('Active: mu_i =%2.2f, mu_o =%2.2, G(I)max =%2.2, f(U)max =%2.2\n',...
      mu_i,mu_o,max(abs(g_alpha(:))),max(abs(dUin(:))));
    maxGUc = [max(abs(g_alpha(:))),max(abs(dUin(:))),max(abs(dcpl(:)))] %#ok<NOPRT>
    
    dt0   = epsilon*0.5;
    dt_a  = dt0 / max(abs(dphi(:)));  
    phi   = phi + dt_a * dphi;
    phi   =  reinitializeLevelSetFunction(phi,1,1,2,2,2,true() ); 
    
  end

function  [phi dt_a mu_i mu_o] = update_not_active_phi( Img, phi, phiActive, U)
    [mu_i, mu_o] = compute_means(Img,phi);
    kappa_phi(:) = kappa(phi,1:numel(phi));
    g_alpha      = (Img - mu_i).^2 - (Img - mu_o).^2 ;
    dphi         = 0*(-g_alpha + lambda * kappa_phi) ; 
    
    gammaSqr = 0.25;
    dcpl     =  ( Heavi(phiActive) ).*(-gammaSqr);
    dUin     =  abs(U).* (Heavi(U) - Heavi(phi));
    dphi     = delta(phi) .* (dphi + dcpl + dUin);
    
     fprintf('Inactive: mu_i =%2.2f, mu_o =%2.2, G(I)max =%2.2, f(U)max =%2.2\n',...
       mu_i,mu_o,max(abs(g_alpha(:))),max(abs(dUin(:))));
     maxGUc = [max(abs(g_alpha(:))),max(abs(dUin(:))),max(abs(dcpl(:)))] %#ok<NOPRT>

    dt0   = epsilon*0.5;
    dt_a  = dt0 / (1e-3+max(abs(dphi(:))));  
    phi   = phi + dt_a * dphi;
    phi   =  reinitializeLevelSetFunction(phi,1,1,2,2,2,true() ); 
    
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
    
    Cdisp= C21.*C12;
    imgb( abs(Cdisp)>0 ) = (imgb(abs(Cdisp)>0)/2 + abs(C21(abs(Cdisp)>0))/max(abs(Cdisp(:))) );
    imgb( abs(Cdisp)>0 ) = (imgb(abs(Cdisp)>0)/2 + abs(C12(abs(Cdisp)>0))/max(abs(Cdisp(:))) );
    
    img_show(:,:,1) = imgr; img_show(:,:,2) = imgg; img_show(:,:,3) = imgb;
    img_show(img_show>1)=1; img_show(img_show<0)=0;
    sh=sfigure(1); %subplot(1,2,2); 
    imshow(img_show);
    title(['image and contours. ||U_1||_2=' num2str(norm(U1(:))) ', t=' num2str(tt) ]);
    xlabel(['active index = ' num2str_fixed_width(Idx_Active) ]);
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

function [Xphi1, Xphi2, Xphi1_star, Xphi2_star, img_out] = load_IKHJ( )
    addpath('./IK_HJ_MultiBrainEx');
    imgVol=mha_read_volume('t1_crop_dec.mha');

    %the two structures of interest are Label=22 and Label=8
    labVol=mha_read_volume('labels_crop_dec.mha'); 
    sliceInterest=35;
    imgSlice=imgVol(:,:,sliceInterest);
    labSlice=labVol(:,:,sliceInterest);
    imgSlice=imgSlice(1:61,4:64);
    labSlice=labSlice(1:61,4:64);
    
    sz=256;
    
    % set up truth
    Xphi1_star = 2.0*(labSlice==22)-1.0; Xphi1_star=imresize(Xphi1_star,sz/61,'bilinear');
    Xphi2_star = 2.0*(labSlice==8)-1.0;  Xphi2_star=imresize(Xphi2_star,sz/61,'bilinear');
    Xphi1_star = reinitializeLevelSetFunction(1e2*Xphi1_star,1,1,20,2,2,true() ); 
    Xphi2_star = reinitializeLevelSetFunction(1e2*Xphi2_star,1,1,20,2,2,true() ); 
    
    img_out=imresize(double(imgSlice),sz/61,'bicubic');
    img_out=(img_out-min(img_out(:)))/(max(img_out(:))-min(img_out(:)));
    
    Xphi1 = Xphi1_star-max(Xphi1_star(:))+20.0;
    Xphi2 = Xphi2_star;
%     Xphi1(Xphi1>2)=2;
%     Xphi2(Xphi2>2)=2;
%     Xphi1(Xphi1<-2)=-2;
%     Xphi2(Xphi2<-2)=-2;
    
    fprintf('');
    
end







