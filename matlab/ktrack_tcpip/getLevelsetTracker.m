function  tkr = getLevelsetTracker( params )  
% minimal params:  Img  HxWxC, control_is_on bool
% returns: tkr struct with update_phi method

  addpath('~/source/chernobylite/matlab/util/');
  addpath('~/source/chernobylite/matlab/display_helpers/');
  addpath('~/source/chernobylite/matlab/LSMlibPK/');
  dbstop if error;
  tkr = []; 
  tkr.CONTROL_IS_ON = params.control_is_on; 
  
  img_in       = params.Img;
  [m,n,cc]     = size(img_in); 
  [xx yy]      = meshgrid(linspace(-1,1,n),linspace(-1*m/n,1*m/n,m));
  tkr.img_size = [m,n,cc]; 
  tkr.xx       = xx;
  tkr.yy       = yy;
  
  [epsilon, dX]= get_params();
  tkr.phi      = (100)*(0.02 - xx.^2 - yy.^2);
  tkr.phi      = reinitializeLevelSetFunction(tkr.phi, 200, dX, 2, 3, 2, true() );
  tkr.U        = 0*tkr.phi;
  tkr.lambda   = 0.1 * ( max((img_in(:))) - min((img_in(:))) );
  
  % METHODS 
  tkr.update_phi = @update_phi;
  tkr.display    = @displayLevelSets; 
  tkr.get_center = @get_center;
  
  % TESTING
  img_in = max(img_in(:)) * rgb2gray(img_in);
  tkr.display(img_in);
  [dt_a mu_i mu_o g_alpha] = tkr.update_phi( img_in ); %#ok<ASGLU,NASGU>
  tkr.display(img_in); 
  
  % DONE
  fprintf('done, created tkr!\n');
  
  %---------------------------------------------%
  
  % TODO: add func to shift and initialize phi with offset
  
  function [xyF] = get_center( )
    [yc xc] = find( tkr.phi >=  max(tkr.phi(:)), 1);
    xyF= [xc;yc];
  end
  
  function  [dt_a mu_i mu_o g_alpha] = update_phi(Img)
    if ndims(Img) == 3
      Img = rgb2gray(Img);
    end
    phi           = tkr.phi;
    CONTROL_IS_ON = tkr.CONTROL_IS_ON; 
    U             = tkr.U;
    [mu_i, mu_o]  = compute_means(Img,phi);
    
    kappa_phi     = 0*phi;
    kappa_phi(1:numel(phi)) = kappa(phi,1:numel(phi));
    
    g_alpha= -(Img - mu_i).^2 + (Img - mu_o).^2;
    
    
    if CONTROL_IS_ON 
      [U_i, U_o] = compute_means(U,phi);
      uMin   =        min( 0.1, abs(U_i-U_o) );
      xi_sqr = (U-U_i).^2; 
      f_of_U = (xi_sqr > 0.1).*( min(xi_sqr - uMin, 1 ) ).*(-g_alpha);
      % f: |U-U_i| < uMin  maps to 0
      % f: |U-U_i| >=1  maps to -G 
    else
      f_of_U = 0*phi;
    end
    lambda= tkr.lambda;
    dphi  = delta(phi) .* ( g_alpha + f_of_U + lambda * kappa_phi) ;
 
    dt0   = 0.8;
    dt_a  = dt0 / max(abs(dphi(:)));  
    phi   = phi + dt_a * dphi;
    
    phi    = reinitializeLevelSetFunction(phi, 2, dX, 2, 3, 2, true() );
    tkr.phi= phi; % must force copy !?? 
    
  end

  
  function displayLevelSets(img0)
    if ndims(img0) ~= 3 
      img_show = repmat(img0,[1 1 3]);
    else
      img_show = img0;
    end
    [epsilon,dX] = get_params();
    phi_show_thresh = epsilon;
    imgb = img_show(:,:,3);
    imgg = img_show(:,:,2);
    imgr = img_show(:,:,1);
    
    tkr.phi = tkr.phi;
    imgr( abs( tkr.phi ) < phi_show_thresh ) = 0;

    imgb( abs( tkr.phi ) < phi_show_thresh ) = 0;
    imgg( abs( tkr.phi ) < phi_show_thresh) = (imgg( abs( tkr.phi ) < phi_show_thresh) .* ...
      abs( tkr.phi(abs(tkr.phi) < phi_show_thresh ) )/phi_show_thresh  + ...
      1.5 * (phi_show_thresh - abs( tkr.phi(abs(tkr.phi) < phi_show_thresh ) ) )/phi_show_thresh );
      
    img_show(:,:,1) = imgr; img_show(:,:,2) = imgg; img_show(:,:,3) = imgb;
    img_show(img_show>1)=1; img_show(img_show<0)=0;
    
    imshow(img_show);
    
  end

  function  [mu_i mu_o] = compute_means( Img,phi )
    mu_i = trapz(trapz(Heavi( phi ) .* Img)) / trapz(trapz(Heavi( phi ) ) );
    mu_o = trapz(trapz( (1-Heavi( phi )) .* Img)) / trapz(trapz( (1-Heavi( phi )) ) );
  end

  function [epsilon,dX] = get_params()
    epsilon   = sqrt(2);
    dX        = 1/sqrt(2);
  end
  
  function z = Heavi(z)
    [epsilon,dX] = get_params();
     z = 1 * (z >= epsilon) + (abs(z) < epsilon).*(1+z/epsilon+1/pi * sin(pi*z/epsilon))/2.0;
  end
  
  function z = delta(z)
    [epsilon,dX] = get_params();
    z = (abs(z) < epsilon).*(1 + cos(pi*z/epsilon))/(epsilon*2.0);
  end
  
end
      
%     fprintf('mu_i = %f, mu_o = %f, g_alpha max = %f, lam*kap max = %f,',...
%       mu_i,mu_o,max(abs(g_alpha(:))),max(abs(lambda*kappa_phi(:))));
   
