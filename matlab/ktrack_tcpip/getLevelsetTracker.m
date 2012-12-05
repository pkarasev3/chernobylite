function  tkr = getLevelsetTracker( params )  
% minimal params:  Img  HxWxC, control_is_on bool
% returns: tkr struct with update_phi method

  %dbstop if error;
  addpath('~/source/chernobylite/matlab/util/');
  addpath('~/source/chernobylite/matlab/display_helpers/');
  addpath('~/source/chernobylite/matlab/LSMlibPK/');
  
  if nargin < 1
    params = []; params.control_is_on = false; params.Img = zeros(480,640,3);
  end
  
  tkr = []; 
  tkr.CONTROL_IS_ON = params.control_is_on; 
  
  img_in       = params.Img;
  [m,n,cc]     = size(img_in); 
  [xx yy]      = meshgrid(linspace(-1,1,n),linspace(-1*m/n,1*m/n,m));
  tkr.img_size = [m,n,cc]; 
  tkr.xx       = xx;
  tkr.yy       = yy;
  tkr.fx       = 0*xx;
  tkr.fy       = 0*yy;
  
  tkr.img0     = params.Img;
  tkr.img1     = params.Img;
  tkr.img_show = zeros(m,n,3);
  tkr.xyF      = [n/2;m/2];
  
  [epsilon, dX]= get_params();
  tkr.phi      = (100)*(0.005 - xx.^2 - yy.^2);
  tkr.phi      = reinitializeLevelSetFunction(tkr.phi, 2,dX, 2, 2, true() );
  tkr.U        = 0*tkr.phi;
  tkr.lambda   = 0.35 * ( max((img_in(:))) - min((img_in(:))) );
  tkr.xyF_prev = [tkr.img_size(2); tkr.img_size(1)];
  tkr.g_f2f    = eye(4,4);
  tkr.f        = 500.0;
  
  tkr.psi      = tkr.phi;
  tkr.phi0     = tkr.phi;
  tkr.phi1     = tkr.phi;
  
  % METHODS
  tkr.Heavi      = @Heavi;
  tkr.delta      = @delta;
  tkr.get_means  = @compute_means;
  tkr.update_phi = @update_phi;
  tkr.display    = @displayLevelSets; 
  tkr.get_center = @get_center;
  tkr.compensate = @apply_g_compensation;
  
  bTestShifter = false;
  if bTestShifter
    test_shifter(tkr.phi);
  end
  
  % TESTING 1
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


  
  function  [dt_a mu_i mu_o g_alpha] = update_phi(Img, U, t_0to1 )
    global TKR;
    global KOpts;
    if ndims(Img) == 3
      Img = rgb2gray(Img);
    end
    if nargin < 3
      t_0to1 = 0;
    end
    phi           = tkr.phi;
    CONTROL_IS_ON = tkr.CONTROL_IS_ON; 
    
    if nargin < 2 
      U = 0 * phi;
    end
    
    [roi_i, roi_j]= find( phi > -dX*20 );
    left = min(roi_j(:)); right = max(roi_j(:));
    top  = min(roi_i(:)); bottom= max(roi_i(:));
    Img  = Img(top:bottom,left:right);
    phi  = phi(top:bottom,left:right);
    fx   = TKR.fx(top:bottom,left:right); 
    fy   = TKR.fy(top:bottom,left:right);
  
    U    = U(top:bottom,left:right);
    
    [mu_i, mu_o]  = compute_means(Img,phi);
    g_alpha= -(Img - mu_i).^2 + (Img - mu_o).^2;

    [kappa_phi, dx, dy] = kappa( phi, dX );
    if KOpts.incremental_warp
      Nx =  dx ./ sqrt(dx.^2+dy.^2+1e-6) * abs( TKR.xyF(1) - [n/2] );
      Ny = -dy ./ sqrt(dx.^2+dy.^2+1e-6) * abs( TKR.xyF(2) - [m/2] );
      [m,n] = size( TKR.img0 );
%      Kxy= (norm( TKR.xyF(:) - [n/2;m/2] ) / 100.0 ) * 200.0 ;
     % fprintf('Kxy= %4.4f  ; ',Kxy );
      fgain = (1-t_0to1^2) * max(abs(g_alpha(:))) / KOpts.contour_iters;
      fdotN = -( fx .* Nx + fy .* Ny )*fgain;
    else
      fdotN = 0*phi;
    end
    
    
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
    
    if KOpts.incremental_warp
      fprintf('max g_alpha, fdotN = %4.4f, %4.4f\n',...
      max(abs(g_alpha(:))),max( (abs(fdotN(:))) ) );
    end
    dphi  = delta(phi) .* ( g_alpha + fdotN + f_of_U + lambda * kappa_phi) ;
    dt0   = 0.95;
    dt_a  = dt0 / max(abs(dphi(:)));  
    phi   = phi + dt_a * dphi;
    
    phiArea = trapz(trapz(Heavi( phi ) ) );
    minArea = 17*17;
    phi     = phi + (phiArea < minArea)*(1.0);
    
    tkr.phi(top:bottom,left:right) = phi; % must force copy here, not return it
    top = max(1,top-8); bottom = min(tkr.img_size(1),bottom+8);
    left= max(1,left-8); right = min(tkr.img_size(2),right+8);
    redistIters = 2; % 1 or 2, negligible effect for speed
    tkr.phi(top:bottom,left:right) = reinitializeLevelSetFunction( ... 
                   tkr.phi(top:bottom,left:right), 2, dX,redistIters, 2, 1, true() );
                 
  end

  
  function displayLevelSets(img0)
    global TKR; % d'oh 
    global KOpts;
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
    
    phi = tkr.phi;
  
    if KOpts.getPsiTru && isfield(TKR,'psi') % the "true" sdf for target, draw in red
      %psi     = TKR.psi; 
      psi     = TKR.phi0;
      tkr.psi = psi; % force copy ...
      imgg( abs( tkr.psi ) < phi_show_thresh ) = 0;
      imgb( abs( tkr.psi ) < phi_show_thresh ) = 0;
      imgr( abs( tkr.psi ) < phi_show_thresh) = (imgr( abs( tkr.psi ) < phi_show_thresh) .* ...
          abs( tkr.psi(abs(tkr.psi) < phi_show_thresh ) )/phi_show_thresh  + ...
          1.5 * (phi_show_thresh - abs( tkr.psi(abs(tkr.psi) < phi_show_thresh ) ) )/phi_show_thresh );
    end
    
    % draw contour in green for tracker
    imgr( abs( phi ) < phi_show_thresh ) = 0;
    imgb( abs( phi ) < phi_show_thresh ) = 0;
    imgg( abs( phi ) < phi_show_thresh) = (imgg( abs( phi ) < phi_show_thresh) .* ...
      abs( phi(abs(phi) < phi_show_thresh ) )/phi_show_thresh  + ...
      1.5 * (phi_show_thresh - abs( phi(abs(phi) < phi_show_thresh ) ) )/phi_show_thresh );
      
 
    
    img_show(:,:,1) = imgr; img_show(:,:,2) = imgg; img_show(:,:,3) = imgb;
    img_show(img_show>1)=1; img_show(img_show<0)=0;
    
    % need good way to show tkr.f_of_U ...
    % f_of_U = tkr.f_of_U;
    
    %imshow(img_show);
    
    TKR.img1(:)     = TKR.img0(:);
    TKR.img0(:)     = img0(:);
    TKR.img_show(:) = img_show(:);
    
    TKR.phi1(:)     = TKR.phi(:);
    TKR.phi0(:)     = phi(:);
    TKR.phi(:)      = phi(:);
    
  end

  function [g_f2f] = apply_g_compensation(  )

    global TKR; 
    tkr.g_f2f = TKR.g_f2f;
    tkr.f     = TKR.f;
    fprintf('| got f = %4.4f |, ',tkr.f);
    
    m            = tkr.img_size(1);
    n            = tkr.img_size(2);
    [xx yy]      = meshgrid(linspace(1,n,n),linspace(1,m,m));
    
    % !!!! Crucial source of uncertainty
    % Can't warp for large translation if this is off, low error margin
    z0           = -100.0;
    xx           = -(xx - (n-1)/2) / tkr.f * z0;
    yy           =  (yy - (m-1)/2) / tkr.f * z0;
    tkr.xx       = xx;
    tkr.yy       = yy;
    
    g_f2f        = tkr.g_f2f; % nullify the T, unless later we trust it ... % g_f2f(1:3,4) = 0;
    f            = tkr.f;
    if norm( g_f2f - eye(4,4),'fro') < 1e-6 
      fprintf('Not compensating, seems like g == identity\n');
     % return; % nothing to do, frame to frame is identity
    end
    
    assert( (10 < f) && (1e5 > f ) );
    
    % apply the affine warp to tkr.phi 
    u0 = xx;
    v0 = yy;
    uv = (g_f2f^-1) * [ u0(:)'; v0(:)'; z0 * ones(1,numel(v0)); ones(1,numel(v0)) ];
    
    u1 =  z0 * uv(1,:)./uv(3,:) ;
    v1 =  z0 * uv(2,:)./uv(3,:) ;
    
    xx1  = reshape(u1,size(xx));
    yy1  = reshape(v1,size(yy));
    phi1 = interp2(xx,yy,tkr.phi, xx1, yy1,'*linear',-100); 
    dx   = f*( xx1(m/2,n/2)-xx(m/2,n/2) )/z0;
    dy   = f*( yy1(m/2,n/2)-yy(m/2,n/2) )/z0;
    
    tkr.phi0= tkr.phi;
    tkr.phi = phi1;
    
    % Initialization for next ls iters
    TKR.phi0= tkr.phi;
    
    TKR.fx  = xx1-xx ;
    TKR.fy  = yy1-yy ;
    
    zf2f=real(logm(g_f2f));
    roll_ang=abs( zf2f(2,1) ) * 180/pi; 
    fprintf('roll_ang=%4.4f ; ',roll_ang);
    if roll_ang  > 10.0
      fprintf('\nLarge Roll!! =%4.4f \n',roll_ang);
      breakhere=1;
    end

   
    fprintf('pixel shift in contour compensation, dx=%3.3g, dy=%3.3g\n',dx,dy);
  end    




  function  [mu_i mu_o] = compute_means( Img,phi )
    mu_i = trapz(trapz(Heavi( phi ) .* Img)) / trapz(trapz(Heavi( phi ) ) );
    mu_o = trapz(trapz( (1-Heavi( phi )) .* Img)) / trapz(trapz( (1-Heavi( phi )) ) );
  end

  function [epsilon,dX] = get_params()
    epsilon   = sqrt(2);
    dX        = 1/2 * 1/sqrt(2);
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

%     z = real(logm(TKR.g_f2f)) * 180/pi;
%     if abs(z(2,1)) > 3.0
%       breakhere = true; 
%       % Save arrays and g...
%     end
      
%     fprintf('mu_i = %f, mu_o = %f, g_alpha max = %f, lam*kap max = %f,',...
%       mu_i,mu_o,max(abs(g_alpha(:))),max(abs(lambda*kappa_phi(:))));
   % apply it to phi though, really 
    %dx   = u1 - tkr.img_size(2)/2;
    %dy   = v1 - tkr.img_size(1)/2;
    %phi1 = circshift( tkr.phi, [floor(dy), floor(dx)] );
   
%u0 = 0.0; %(tkr.xyF_prev(1) - tkr.img_size(2)/2)*(1/f);
%v0 = 0.0; %(tkr.xyF_prev(2) - tkr.img_size(1)/2)*(1/f);

