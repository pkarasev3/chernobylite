function [phi1 phi2 img_show U tt xx yy] = run_3Dopenloop_chanvese_demo()
% run demo func in-place:
% [phi1 phi2 img_show] = run_lskk_demo();
set(0,'defaultaxesfontsize',16);  
  set(0,'defaulttextfontsize',18);
  set(0,'defaulttextfontname','Arial');
  set(0,'defaultaxesfontweight','bold');
  set(0,'defaultlinelinewidth',2);
  set(0,'defaultlinemarkersize',4);
addpath('~/source/chernobylite/matlab/util/');
addpath('~/source/chernobylite/matlab/display_helpers/');
addpath('~/source/chernobylite/matlab/LevelSetMethods/');
addpath('~/source/chernobylite/matlab/LSMlibPK/');

  load synth_data_phantom3D  % data file must exist from kslice synthetic data generator!
  disp(header'); % print the would-be mha header
  whos img3D
  
  img = double(img3D);
  phi = zeros(size(img));
  [M N K] = size(phi);
  [xx yy zz] = meshgrid( linspace(-1,1,N), linspace(-1,1,M), linspace(-1,1,K) );
  phi = 0.1 - sqrt(xx.^2+(yy+0.25).^2+5*zz.^2);
  dX  = 2/M;
  phi = reinitializeLevelSetFunction(phi,1,dX,2,3,3,false() );
  
  sfigure(1); imagesc( img(:,:,K/2)+ max(img(:))*(phi(:,:,K/2)>=0 ) );
  
  prev_phi  = phi;
  [phi ta mui1 muo1]  = update_phi( img, phi, 0*phi);
  
  mu1_in_out = [mu1_in_out, [mui1;muo1]];
  
  % setup display image
  displayLevelSets();
  fprintf('');
  
  fprintf('done! saving .... \n');


  

end


function  [mu_i mu_o] = compute_means( Img,phi )
     mu_i = trapz(trapz(trapz(Heavi( phi ) .* Img))) / trapz(trapz(trapz(Heavi( phi ) ) ));
     mu_o = trapz(trapz(trapz( (1-Heavi( phi )) .* Img))) / trapz(trapz(trapz( (1-Heavi( phi )) ) ));
  end

  function  [phi dt_a mu_i mu_o] = update_phi( Img, phi, Coupling )
        
    [mu_i, mu_o] = compute_means(Img,phi);
       
    kappa_phi(1:numel(phi)) = kappa(phi,1:numel(phi));
    
    g_alpha= (Img - mu_i).^2 - (Img - mu_o).^2 + Coupling;
    dphi  = delta(phi) .* (-g_alpha + lambda * kappa_phi) ;
    
    dt0   = 0.9;
    dt_a  = dt0 / max(abs(dphi(:)));  
    phi   = phi + dt_a * dphi;
        
    phi   =  reinitializeLevelSetFunction(phi,1,dX,2,3,3,false() );
    
    
    
  end

  function displayLevelSets()
        
  end




