function  tkr = getLevelsetTracker( params )  
% minimal params:  Img  HxWxC, control_is_on bool
% returns: tkr struct with update_phi method

  %dbstop if error;
  addpath('~/source/chernobylite/matlab/util/');
  addpath('~/source/chernobylite/matlab/display_helpers/');
  addpath('~/source/chernobylite/matlab/LSMlibPK/');
  
  tkr = []; 
  tkr.CONTROL_IS_ON = params.control_is_on; 
  
  img_in       = params.Img;
  [m,n,cc]     = size(img_in); 
  [xx yy]      = meshgrid(linspace(-1,1,n),linspace(-1*m/n,1*m/n,m));
  tkr.img_size = [m,n,cc]; 
  tkr.xx       = xx;
  tkr.yy       = yy;
  
  [epsilon, dX]= get_params();
  tkr.phi      = (100)*(0.005 - xx.^2 - yy.^2);
  tkr.phi      = reinitializeLevelSetFunction(tkr.phi, 200, dX, 2, 3, 2, true() );
  tkr.U        = 0*tkr.phi;
  tkr.lambda   = 0.5 * ( max((img_in(:))) - min((img_in(:))) );
  tkr.xyF_prev = [tkr.img_size(2); tkr.img_size(1)];
  tkr.g_f2f    = eye(4,4);
  tkr.f        = 500.0;
  
  % METHODS 
  tkr.update_phi = @update_phi;
  tkr.display    = @displayLevelSets; 
  tkr.get_center = @get_center;
  tkr.compensate = @apply_g_compensation;
  
  % TESTING
  bTesting = false();
  if bTesting && ndims(img_in)==3
    img_in = max(img_in(:)) * rgb2gray(img_in);
    tkr.display(img_in);
    [dt_a mu_i mu_o g_alpha] = tkr.update_phi( img_in, tkr.U ); %#ok<ASGLU,NASGU>
    tkr.display(img_in); 
  end
  % DONE
  fprintf('done, created tkr!\n');
  
  %---------------------------------------------%
  

  
  function [xyF] = get_center( )
    [yc xc] = find( tkr.phi > 0 ); 
    if isempty(xc) || isempty(yc) 
      fprintf('warning, no phi > 0 region !?\n'); 
      xyF = [tkr.img_size(2)/2 ; tkr.img_size(1)/2];
      return;
    end % strictly speaking mean() is flawed, but OK if blob target like warped circle
    xyF= [mean(xc(:));mean(yc(:))];
    tkr.xyF_prev = xyF;
  end


  
  function  [dt_a mu_i mu_o g_alpha] = update_phi(Img, U)
    if ndims(Img) == 3
      Img = rgb2gray(Img);
    end
    phi           = tkr.phi;
    CONTROL_IS_ON = tkr.CONTROL_IS_ON; 
    
    if nargin < 2 
      U = 0 * phi;
    end
    
    [roi_i, roi_j]= find( phi > -dX*10 );
    left = min(roi_j(:)); right = max(roi_j(:));
    top  = min(roi_i(:)); bottom= max(roi_i(:));
    Img  = Img(top:bottom,left:right);
    phi  = phi(top:bottom,left:right);
    U    = U(top:bottom,left:right);
    
    [mu_i, mu_o]  = compute_means(Img,phi);
    
    kappa_phi     = 0*phi;
    kappa_phi(1:numel(phi)) = kappa(phi,1:numel(phi));
    
    g_alpha= -(Img - mu_i).^2 + (Img - mu_o).^2;
    
    
    if CONTROL_IS_ON  
      U = imfilter( (-1 + 2*(U>0)), fspecial('gaussian',[5 5],3),'replicate');
      [U_i, U_o] = compute_means(U,phi);
      % f: |U-U_i| < uMin  maps to 0
      % f: |U-U_i| >=1  maps to -G 
      uMin   =        0.1;
      xi_sqr = abs(U-U_i); 
      tkr.f_of_U = (xi_sqr > uMin).*( min(xi_sqr - uMin, 1.0 ) );
    else
      tkr.f_of_U = 0*phi;
    end
    f_of_U     = (tkr.f_of_U).*(-g_alpha);
    lambda     = tkr.lambda;
    
    
    dphi  = delta(phi) .* ( g_alpha + f_of_U + lambda * kappa_phi) ;
    dt0   = 0.9;
    dt_a  = dt0 / max(abs(dphi(:)));  
    phi   = phi + dt_a * dphi;
    
    tkr.phi(top:bottom,left:right) = phi; % must force copy here, not return it
    top = max(1,top-8); bottom = min(tkr.img_size(1),bottom+8);
    left= max(1,left-8); right = min(tkr.img_size(2),right+8);
    redistIters = 2; % 1 or 2, negligible effect for speed
    tkr.phi(top:bottom,left:right) = reinitializeLevelSetFunction( ... 
                   tkr.phi(top:bottom,left:right), 2, dX,redistIters, 2, 1, true() );
    
  end

  
  function displayLevelSets(img0)
    if ndims(img0) ~= 3 
      img_show = repmat(img0,[1 1 3]);
      img_show = img_show / max(img_show(:));
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
    
    % need good way to show tkr.f_of_U ...
    f_of_U = tkr.f_of_U;
    
    imagesc(img_show);
    
    
  end

  function [g_f2f] = apply_g_compensation(  )

    global TKR; 
    tkr.g_f2f = TKR.g_f2f;
    tkr.f     = TKR.f;
    
    % nullify the T, unless later we trust it ...
    g_f2f        = tkr.g_f2f;
    g_f2f(1:3,4) = 0;
    f            = tkr.f;
    if norm( g_f2f - eye(4,4),'fro') < 1e-6 
      return; % nothing to do, frame to frame is identity
    end
    
    assert( (10 < f) && (1e5 > f ) );
    
    % apply the affine warp to tkr.phi 
    u0 = 0.0; %(tkr.xyF_prev(1) - tkr.img_size(2)/2)*(1/f);
    v0 = 0.0; %(tkr.xyF_prev(2) - tkr.img_size(1)/2)*(1/f);
    z0 = 1.0;
    uv = (g_f2f^(-1)) * [ u0(:)'; v0(:)'; z0 * ones(1,numel(v0)); ones(1,numel(v0)) ];
    
    x0 =  f * uv(1,:)./uv(3,:) + tkr.img_size(2)/2;
    y0 = -f * uv(2,:)./uv(3,:) + tkr.img_size(1)/2;
    
    % apply it to phi though, really 
    dy   = y0 - tkr.img_size(1)/2;
    dx   = x0 - tkr.img_size(2)/2;
    phi1 = circshift( tkr.phi, [floor(dy), floor(dx)] );
    tkr.phi0= tkr.phi;
    tkr.phi = phi1;
    
    xy0            = [x0,y0];
    pixel_shift    = [dx,dy];  % SAVE IT!
    fprintf('pixel shift in contour compensation, dx=%3.3g, dy=%3.3g\n',dx,dy);
    
  end    
    
    
%     u0 =  (xy0(1) - 640/2) * (1/f);
%     v0 = -(xy0(2) - 480/2) * (1/f);
% 
%     uv = (g_f2f^(-1)) * [ u0(:)'; v0(:)'; ones(1,numel(v0)); ones(1,numel(v0)) ];
% 
%     x0 =  f * uv(1,:)./uv(3,:) + 640/2;
%     y0 = -f * uv(2,:)./uv(3,:) + 480/2;
% 
%     xy0_old = xy0;
%     xy0     = [x0,y0];
%     img_coord_diff = norm( xy0(:) - xy0_old(:) );



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
   
