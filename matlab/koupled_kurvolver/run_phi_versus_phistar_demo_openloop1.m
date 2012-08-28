function run_phi_versus_phistar_demo_openloop1()
  set(0,'defaultaxesfontsize',16);  
  set(0,'defaulttextfontsize',18);
  set(0,'defaulttextfontname','Arial');
  set(0,'defaultaxesfontweight','bold');
  set(0,'defaultlinelinewidth',2);
  set(0,'defaultlinemarkersize',4);
  
  dt_init         = 0.35; 
    
  %[Dval_alla t_alla] = run_core( sqrt(1/(2)) , dt_init);
  %sfigure(2); semilogy(t_alla,Dval_alla,'--','color',[0 0 0.8]); hold on;  
  
  %[Dval_allb t_allb] = run_core(     (1/(2)) ,dt_init);
  % sfigure(2); semilogy(t_allb,Dval_allb,'-.','color',[0 0.4 .6]); hold on;  
   
  [Dval_allc t_allc] = run_core(     (1/(4))^4 ,dt_init);
  
  % sfigure(2); semilogy(t_allc,Dval_allc,'--','color',[0 0.8 .2]); hold on;  
  
%    [Dval_alld t_alld] = run_core(     (1/(16)) ,dt_init);
%    sfigure(2); semilogy(t_alld,Dval_alld,'-.','color',[0.6 0.2 .2]); hold on;  
%    
%   [Dval_alle t_alle] = run_core(     (1/(256)) ,dt_init);
%    sfigure(2); semilogy(t_alle,Dval_alle,'--','color',[0.9 0.4 .2]); 
   
%    legend('\rho=(1/2)^{1/2}','\rho=(1/2)','\rho=(1/4)','\rho=(1/16)','\rho=(1/256)');
%    xlabel('time (sec)');
%    ylabel('labeling error');
%    title('Labeling Error: D(\phi,\phi^*)'); grid on;
%    
%    hold off;  
%    
%    save run_phi_versus_phistar_demo_BAK
  
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
addpath('~/source/chernobylite/matlab/LSMlibPK/');

img      = phantom(); 
img(img==0) = 0.1;
[star_i star_j] = find( ( (img < 1e-1) - (img == 0) ) > 1e-3  );
keep_idx  = find( star_j < 128 );
star_i    = star_i(keep_idx);
star_j    = star_j(keep_idx);
idx_left  = sub2ind( size(img), star_i, star_j);
left_stub = img*0-1; left_stub( idx_left ) = 1;

before_bridge        = img(142:160,118:147);
img(142:160,118:147) =  img(150,124); % Make a 'bridge' connecting the two chunks
[bx by] = meshgrid( linspace(-1,1,size(before_bridge,2)),linspace(-1,1,size(before_bridge,1)));
bridge_weight = min( 0.5+0*bx, max(max( abs(bx), abs(by) ),0*bx) );
img(142:160,118:147) = img(142:160,118:147).*(1-bridge_weight) + before_bridge.*(bridge_weight);

img = img + (randn(size(img))*5e-2); 
img = abs(img+0.1).^(1.5);
img(img>1)=1; 

[m n] = size(img);
[xx yy] = meshgrid(1:m,1:n);

xy2     = [100,120]; xy1          = [120,92];

RadInit = 25;%15;
d2      = RadInit - sqrt( ((xx-xy2(1))).^2 + ((yy-xy2(2))).^2 );
bUseLSM = true();
if( exist('initial_data_phi_phi-star_img.mat','file' ) ) 
  load initial_data_phi_phi-star_img.mat
  assert( numel(phi2)==numel(img) ); assert( numel(phi_star)==numel(img) ); %#ok<NODEF>
else
  phi_star  = 1e2*tanh(imfilter( left_stub, fspecial('gaussian',[3 3],1.5),'replicate') );
  phi2      = 1e2*tanh(d2);
  fprintf('initializing level sets... \n');
  if(bUseLSM)
     % phi, "ghost">=1,dX,iters,spatial order, time order
     ghost_width  = 32; dX = 1/sqrt(2); spatial_order = 5; tvdrk_order = 3; iters = 500;
     phi_init = phi_star;
     phi_star = reinitializeLevelSetFunction(phi_init, ghost_width,dX, iters, spatial_order, tvdrk_order);
     phi_init = phi2;
     phi2     = reinitializeLevelSetFunction(phi_init, ghost_width, dX, iters, spatial_order, tvdrk_order);
  else
    phi_star  = reinit_SD(phi_star, 1/sqrt(2), 1/sqrt(2), 0.35, 'ENO3', 500);
    phi2      = reinit_SD(phi2,1/sqrt(2),1/sqrt(2),0.35,'ENO3',500);
  end
  
  
end

phi1 = phi_star;

if( ~exist('initial_data_phi_phi-star_img.mat','file' ) )
  save('initial_data_phi_phi-star_img.mat','phi_star','phi2','img');
end


sfigure(1); clf;

epsilon   = sqrt(2); %1.1;%0.8;
Heavi     = @(z)  1 * (z >= epsilon) + (abs(z) < epsilon).*(1+z/epsilon+1/pi * sin(pi*z/epsilon))/2.0;
delta     = @(z)  1 * (z == 0) + (abs(z) < epsilon).*(1 + cos(pi*z/epsilon))/(epsilon*2.0);

x=linspace(-3,3,100);
sfigure(1); plot(x,Heavi(x),'r--'); hold on; plot(x,delta(x),'b-.'); hold off;



tt = 0;
% DEMO 1:  Just bridge, original image
img  = img - mean(img(:));
img  = imfilter( 255*img, fspecial('gaussian',[5 5],0.0625),'replicate');
img0 = img-min(img(:)); 
img0 = img0/max(img0(:)); 

img_show_mid  = img * 0;
img_show_init = img * 0;
phi2_init     = 0 * img;
phi2_mid      = 0 * img;

lambda     = 0.25*max(abs(img(:)))^2;
kappa_phi  = 0*phi1;
delta_rel1 = [1];
delta_rel2 = [1];
delta_abs1 = [1];
delta_abs2 = [1];

t_all      = [0];

phi_show_thresh = max([0.95,epsilon]);

steps           = 0;
MaxSteps        = 500;
Dval            = eval_label_dist(phi1,phi2);

Dval_all        = [Dval];


Gmax            = (max(img(:))-min(img(:))).^2; % maximum the G(\phi,I) term can ever be
lambda          = 0.03 * Gmax;
dt0             = dt_init;
MaxTime         = 2.0;

imgFunc_all     = [MeansCost(img,phi2,lambda,Heavi)];

while( (tt < MaxTime) && (steps < MaxSteps) )

  phi1                 = phi_star;
  
  [phi2 tb]   = update_phi( img, phi2, 0*phi2,2);
  
  imgFunc_all = [imgFunc_all , MeansCost(img,phi2,lambda,Heavi) ];
  tt          = tt + tb;
  t_all       = [t_all, tt];
  
  % setup display image
  displayLevelSets();

  
  fprintf('');
  
  steps = steps+1;
end

result = save_all( );
fprintf('result = %f \n',result);

  function res = save_all( )
    fprintf('done! saving .... \n');
    save run_openloop_bridge_demo t_all ... 
                          imgFunc_all   phi2_init phi2_mid img_show_mid img_show_init phi1 phi2 img img_show  tt xx yy  steps Dval_all 
    setenv('rhoval',num2str(rho_argin))
    !cp -v run_openloop_bridge_demo.mat  "bridge_demo_rho=${rhoval}_`date +%d%b%Y-%H-%M`.mat"
    res = 1;                               
  end
  function [Dval] = eval_label_dist( phiA, phiB )                             
    
    Dval = 0.5 * trapz(trapz( (Heavi(phiA)-Heavi(phiB)).^2 ) );
    
  end
                             
  function  [phi dt_a g_source] = update_phi( Img, phi, Coupling, ... 
                                               redist_iters)
    
    
    
    mu_i = trapz(trapz(Heavi( phi ) .* Img)) / trapz(trapz(Heavi( phi ) ) );
    mu_o = trapz(trapz( (1-Heavi( phi )) .* Img)) / trapz(trapz( (1-Heavi( phi )) ) );
    
    %kappa_phi(1:numel(phi)) = anti_kappa(phi,1:numel(phi));
    kappa_phi(1:numel(phi)) = kappa(phi,1:numel(phi));
    
    GofIandPhi = (Img - mu_i).^2 - (Img - mu_o).^2;
    assert( max(abs(GofIandPhi(:)))  <= Gmax );
    g_alpha = GofIandPhi + Coupling;
    
    g_source= -g_alpha + 0 * lambda * kappa_phi ;
    dphi  = delta(phi) .* (-g_alpha + lambda * kappa_phi) ;
    
    
    fprintf('mu_i = %f, mu_o = %f, g_alpha max = %f, lam*kap max = %f,',...
      mu_i,mu_o,max(abs(g_alpha(:))),max(abs(lambda*kappa_phi(:))));
    
    dt_a  = 1 / max(abs(dphi(:)));  % can go above 1 but then rel-step gets jagged...
    if( redist_iters > 0 )
      if( bUseLSM )
        ghost_width  = 16; dX = 1/sqrt(2); spatial_order = 5; tvdrk_order = 3;
        levelset_rhs =  computeLevelSetEvolutionEqnRHS(-phi, {-dt_a*dphi}, ghost_width, dX, spatial_order );
        phi_rd       =  reinitializeLevelSetFunction(phi+dt0*levelset_rhs,1,dX,...
                             redist_iters,spatial_order,tvdrk_order,false() );
        delta_phi_redist = norm(phi_rd-phi,'fro'); 
        fprintf('.. |dPhi_dt|=%f .. ', delta_phi_redist); 
        phi          = phi_rd;
      else
        phi   = phi + dt_a * dphi;
        phi =  reinit_SD(phi, 1, 1, dt0, 'ENO3', redist_iters);
      end
    end
    fprintf('');
    
    
  end

  function displayLevelSets()
    
    img_show = draw_contour_pair_on_image( img0, phi1, phi2, phi_show_thresh);
    img_show = img_show(9:(end-8),33:(end-32),:); % truncate ROI

    sfigure(1); imshow(img_show); axis equal;
    title(['image and contours, t=' num2str(tt), ' steps = ' num2str_fixed_width(steps) ]);

    fprintf( 'max-abs-phi = %f, t= %f, steps = %d \n',max(abs(phi1(:))),tt, steps);
    
    % write images to disk ??
    imwrite(img_show,['openloop_bridge_demo_' num2str_fixed_width(steps) '.png']);
   
   
    sfigure(2); plot(t_all,imgFunc_all-imgFunc_all(1),'m-.'); xlabel('time'); ylabel('image functional');
    
  
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



