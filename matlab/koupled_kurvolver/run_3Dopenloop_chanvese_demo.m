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
  img = zeros(size(img3D,1)/2, size(img3D,2)/2, size(img3D,3)/2 );
  for k =1:(size(img3D,3)/2)
    img(:,:,k) = imresize( double( 0.5*(img3D(:,:,k*2)+img3D(:,:,k*2-1))),0.5 );
  end
  img = img - min(img(:));
  img = 1000*img/max(img(:));
  whos img
  
  phi = zeros(size(img));
  [M N K] = size(phi);
  [xx yy zz] = meshgrid( 20*linspace(-1,1,N), 20*linspace(-1,1,M), 20*linspace(-1,1,K) );
  phi = 2 - sqrt(xx.^2+(yy+8).^2+5*zz.^2);
  dX  = 2/M;
  phi = reinitializeLevelSetFunction(phi,1,dX,2,3,3,false() );
  displayLevelSets(img,phi);
  Gmax = (max(img(:))-min(img(:)))^2; fprintf('Gmax=%f\n',Gmax);
  
  lambda = 0.05*Gmax;
  
  maxIters = 100;
  for iters = 1:maxIters
  
    prev_phi  = phi;
    
    for inner_iters = 1:5
      [phi ta mui1 muo1]  = update_phi( img, phi, 0*phi,lambda); %#ok<ASGLU>
      E = MeansCost(img,phi,lambda); Ediff = (mui1-muo1)^2;
      fprintf('mu_i=%g, mu_o=%g, E_cv=%g, Ediff=%g \n',mui1,muo1,E,Ediff);
    end
    
    redistIters = 2;
    phi         = reinitializeLevelSetFunction(phi,1,dX,redistIters,3,3,false() );
    
    % setup display image
    displayLevelSets(img,phi); sfigure(1); title(sprintf('E=%g, u_i=%g, u_o=%g',E,mui1,muo1));
    fprintf('%d of %d CV openloop iters ... \n',iters, maxIters);

    fprintf('');
  end


  

end


function e=getEpsilon()
  e   = sqrt(2);
end

function y=Heavi(z)
  epsilon=getEpsilon();
  y=1 * (z >= epsilon) + (abs(z) < epsilon).*(1+z/epsilon+1/pi * sin(pi*z/epsilon))/2.0;
end

function y=delta(z)
  epsilon   = getEpsilon();
  y=1 * (z == 0) + (abs(z) < epsilon).*(1 + cos(pi*z/epsilon))/(epsilon*2.0);
end




function  [mu_i mu_o] = compute_means( Img,phi )
     mu_i = trapz(trapz(trapz(Heavi( phi ) .* Img))) / trapz(trapz(trapz(Heavi( phi ) ) ));
     mu_o = trapz(trapz(trapz( (1-Heavi( phi )) .* Img))) / trapz(trapz(trapz( (1-Heavi( phi )) ) ));
end

function E = MeansCost(Img,phi,lambda)
   [mu_i, mu_o] = compute_means(Img,phi);
   
   term1 = trapz(trapz(trapz(  Heavi(phi).*(Img-mu_i).^2 ) ) );
   term2 = trapz(trapz(trapz(  (1-Heavi(phi)).*(Img-mu_o).^2 ) ) );
   [gx gy gz] = gradient( Heavi(phi) );
   term3 = lambda*trapz(trapz(trapz(  sqrt(gx.^2+gy.^2+gz.^2) ) ) );
   
   E = term1+term2+term3;
   
end
  

  function  [phi dt_a mu_i mu_o] = update_phi( Img, phi0, Coupling, lambda)
  
    [mu_i, mu_o] = compute_means(Img,phi0);
    
    kappa_phi = 0*phi0;
    kappa_phi(1:numel(phi0)) = kappa3(phi0,1:numel(phi0));
    
    g_alpha= (Img - mu_i).^2 - (Img - mu_o).^2 + Coupling;
    dphi  = delta(phi0) .* (-g_alpha + lambda * kappa_phi) ;
    
    dAdt    = trapz(trapz(trapz( abs(delta(phi0).*dphi ) ) ) );
    maxRate = 10.0;
    
    dt_a  = maxRate / dAdt;
    phi   = phi0 + dt_a * dphi;
        
    dAdt_    = trapz(trapz(trapz( Heavi(phi) ) ) ) - trapz(trapz(trapz(Heavi(phi0))));
    
    fprintf('');
    
  end

  function displayLevelSets(img,phi)
    [M N K] = size(img);
    
    slice0 = img(:,:,K/2)+ max(img(:))*(phi(:,:,K/2)>=0 );
    slice1 = img(:,:,K/2+5)+ max(img(:))*(phi(:,:,K/2+5)>=0 );
    slice2 = img(:,K/2,:)+max(img(:))*(phi(:,K/2,:) >=0 );    
    slice2 = reshape( slice2, [M N] );
    
    sfigure(1); 
    subplot(1,3,1); imagesc( slice0 );    colormap('bone');
    subplot(1,3,2); imagesc( slice1 );    colormap('bone');
    subplot(1,3,3); imagesc( slice2 );    colormap('bone');
    subplot(1,3,1); 
    fprintf(''); drawnow;
  end




