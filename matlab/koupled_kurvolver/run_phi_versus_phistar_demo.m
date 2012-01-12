function run_phi_versus_phistar_demo()
  set(0,'defaultaxesfontsize',16);  
  set(0,'defaulttextfontsize',18);
  set(0,'defaulttextfontname','Arial');
  set(0,'defaultaxesfontweight','bold');
  set(0,'defaultlinelinewidth',2);
  set(0,'defaultlinemarkersize',4);
  
  dt_init         = 0.7; 
  
  [Dval_alla t_alla] = run_core( sqrt(1/(2)) , dt_init);
  sfigure(2); semilogy(t_alla,Dval_alla,'--','color',[0 0 0.8]); hold on;  
  
  [Dval_allb t_allb] = run_core(     (1/(2)) ,dt_init);
   sfigure(2); semilogy(t_allb,Dval_allb,'-.','color',[0 0.4 .6]); hold on;  
   
  [Dval_allc t_allc] = run_core(     (1/(4)) ,dt_init);
   sfigure(2); semilogy(t_allc,Dval_allc,'--','color',[0 0.8 .2]); hold on;  
  
   [Dval_alld t_alld] = run_core(     (1/(16)) ,dt_init);
   sfigure(2); semilogy(t_alld,Dval_alld,'-.','color',[0.6 0.2 .2]); hold on;  
   
  [Dval_alle t_alle] = run_core(     (1/(256)) ,dt_init);
   sfigure(2); semilogy(t_alle,Dval_alle,'--','color',[0.9 0.4 .2]); 
   
   legend('\rho=(1/2)^{1/2}','\rho=(1/2)','\rho=(1/4)','\rho=(1/16)','\rho=(1/256)');
   xlabel('time (sec)');
   ylabel('labeling error');
   title('Labeling Error: D(\phi,\phi^*)'); grid on;
   
   hold off;  
   
   save run_phi_versus_phistar_demo_BAK
  
end

%load_and_plot_multi_rho( )
function load_and_plot_multi_rho( )
  files = {'bridge_demo_rho=0.0039062_19Dec2011-02-51.mat',...
  'bridge_demo_rho=0.0625_19Dec2011-02-40.mat',...
  'bridge_demo_rho=0.25_19Dec2011-02-22.mat',...
  'bridge_demo_rho=0.5_19Dec2011-02-11.mat',...
  'bridge_demo_rho=0.70711_19Dec2011-02-00.mat'};
  sfigure(2); clf;  xlabel('time (sec)'); ylabel('labeling error');
  N=numel(files);
  for k = 1:N
    data = load(files{k});
    if(mod(k,2)==0)
      sym='-.';
    else
      sym='--';
    end
    if(mod(k,3)==0)
      sym='-';
    end
    semilogy(data.t_all,data.Dval_all,sym,'color',[(0.2+0.8*k/N) (0.4*(-1)^k+0.5)  1-k/N]); hold on;
    fprintf('');
  end
  legend('\rho=1/64','\rho=1/16','\rho=1/4','\rho=1/2','\rho=1/2^{1/2}'); grid on;
  axis([0 0.2 10^(-7) 10^4 ] ); hold off;
  drawnow();
  fprintf('');
end

function [Dval_all t_all phi1 phi2 img_show U tt xx yy] = run_core( rho_argin, dt_init )
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

if( exist('initial_data_phi_phi-star_img.mat','file' ) ) 
  load initial_data_phi_phi-star_img.mat
  assert( numel(phi2)==numel(img) ); assert( numel(phi_star)==numel(img) ); %#ok<NODEF>
else
  phi_star  = 3*tanh( imfilter( left_stub, fspecial('gaussian',[3 3],1.5),'replicate') );
  phi_star  = reinit_SD(phi_star, 1, 1, 0.5, 'ENO3', 300);
  phi2    = tanh(3*d2);
  phi2    = reinit_SD(phi2,1,1,0.5,'ENO3',300);
end

phi1 = phi_star;

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

lambda     = 2*mean(abs(img(:)))^2;
kappa_phi  = 0*phi1;
delta_rel1 = [1];
delta_rel2 = [1];
delta_abs1 = [1];
delta_abs2 = [1];

t_all      = [0];
relTol     = 1e-4;
absTol     = 1e2;
phi_show_thresh = 0.95;
tsum            = 0;
U               = 0 * phi1;
eps_u           = 1e-1;
steps           = 0;
MaxSteps        = 500;
Dval            = eval_label_dist(phi1,phi2);
Dval_all        = [Dval];

Gmax            = (max(img(:))-min(img(:))).^2; % maximum the G(\phi,I) term can ever be

dt0             = dt_init;
MaxTime         = 0.25;
dt_min          = 1e-3;

rho                  =  rho_argin; %(1/2); % 1/16 %(1/4);
U               =  sqrt(Gmax)/sqrt(rho) + 0*img;

while( (tt < MaxTime) && (steps < MaxSteps) )
  
  prev_phi1            = phi1;
  g_gain               = 1e1;
  [phi2_pred ta gval]  = update_phi( img, phi2, 0 * phi2,0);
  phi1                 = phi_star;
  
  assert( sum( U(:) < sqrt( abs(gval(:))/(rho) ) ) == 0 );
  
  
  eta                  = (Heavi(phi1)-Heavi(phi2));
  alph1                = 0.25;
  f_of_U               = -(U.^2).*(eta) - (U.^2).*(eta)./(eta.^2+alph1);
  
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
  if( bBadDval ) % If D failed to shrink, get the previous phi, redistance a bit, reduce the time step
    dt0 = max([dt_min, dt0 - 2e-1*dt0]);
    phi2      = prev_phi2;
    phi2      = reinit_SD(phi2,     1, 1, 0.25, 'ENO3', 2);
    fprintf(', dt0 = %f \n', dt0 );
    fprintf('');
  else           % Raise time-step back up if we're shrinking
    dt0 = min([dt0+1.3e-1*dt0, dt_init]);
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
    save run_openloop_bridge_demo t_all delta_abs1 delta_abs2 delta_rel1 delta_rel2 rho_argin... 
                             phi2_init phi2_mid img_show_mid img_show_init phi1 phi2 img img_show U tt xx yy  steps Dval_all 
    setenv('rhoval',num2str(rho_argin))
    !cp -v run_openloop_bridge_demo.mat  "bridge_demo_rho=${rhoval}_`date +%d%b%Y-%H-%M`.mat"
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
      phi =  reinit_SD(phi, 1, 1, dt0, 'ENO3', redist_iters);
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
      1 * (phi_show_thresh - abs( phi1(abs(phi1) < phi_show_thresh ) ) )/phi_show_thresh );
    
    imgg( abs( phi2 ) < phi_show_thresh) = (imgg( abs( phi2 ) < phi_show_thresh) .* ... 
      abs( phi2(abs(phi2) < phi_show_thresh ) )/phi_show_thresh  + ...
      1 * (phi_show_thresh - abs( phi2(abs(phi2) < phi_show_thresh ) ) )/phi_show_thresh );
    
   
    %imgb( abs(C21)>0 ) = (imgb(abs(C21)>0)/2 + abs(C21(abs(C21)>0))/max(abs(C21(:))) );
    %imgb( abs(C12)>0 ) = (imgb(abs(C12)>0)/2 + abs(C12(abs(C12)>0))/max(abs(C12(:))) );
    
    img_show(:,:,1) = imgr; img_show(:,:,2) = imgg; img_show(:,:,3) = imgb;
    img_show(img_show>1)=1; img_show(img_show<0)=0;
    sh=sfigure(1); subplot(1,2,2); imshow(img_show);
    title(['image and contours, ||U||_2=' num2str(norm(U)) ', t=' num2str(tt), ' steps = ' num2str_fixed_width(steps) ]);
    setFigure(sh,[10 10],3.5,1.9);
    fprintf( 'max-abs-phi = %f, t= %f, steps = %d \n',max(abs(phi1(:))),tt, steps);
    
    %imwrite(img_show,['openloop_bridge_demo_' num2str_fixed_width(steps) '.png']);
        
    
    sfigure(1); subplot(1,2,1);
    semilogy( t_all,Dval_all,'r-.' ); hold on;
    hold off;
    legend('D(\phi,\phi^*)'); xlabel('time (sec)');
    title('labeling error');
    
    
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



