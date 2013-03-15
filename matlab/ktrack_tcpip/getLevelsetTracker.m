function  tkr = getLevelsetTracker( params )  
% minimal params:  Img  HxWxC, bUseUcontrol bool
% returns: tkr struct with update_phi method
  global KOpts;
  %dbstop if error;
  addpath('~/source/chernobylite/matlab/util/');
  addpath('~/source/chernobylite/matlab/display_helpers/');
  addpath('~/source/chernobylite/matlab/LSMlibPK/');
  
  if nargin < 1
    params = []; params.bUseUcontrol = false; params.Img = zeros(480,640,3);
  end
  
  tkr = []; 
  tkr.bUseUcontrol = KOpts.bUseUControl; 
  
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
  tkr.U        = zeros(m,n);
  tkr.imgMinMax= [min(params.Img(:)); max(params.Img(:)) ];
  
  tkr.compare_w0_wX = zeros(3,2);
  tkr.curr_Nframe = 0;% source's count of image (having this is not cheating)
  tkr.prev_Nframe = 0;% source's count of image (having this is not cheating)
  tkr.true_xy  = [n/2;m/2];% for ground-truth eval (this is cheating if used)
  
  [epsilon, dX]= get_params();
  tkr.phi      = (100)*(0.005 - xx.^2 - yy.^2);
  tkr.phi      = reinitializeLevelSetFunction(tkr.phi, 2,dX, 2, 2, true() );
  tkr.fUraw    = 0*tkr.phi;
  tkr.f_of_U   = 0*tkr.phi;
  tkr.lambda   = 5 * ( max((img_in(:))) - min((img_in(:))) );
  tkr.xyF_prev = [tkr.img_size(2)/2; tkr.img_size(1)/2];
  tkr.g_f2f    = eye(4,4);
  tkr.g_ctrl   = eye(4,4);
  tkr.g_k_d    = eye(4,4);
  tkr.f        = 500.0;
  tkr.Area     = sum( tkr.phi(:) >= 0 );

  
  tkr.psi      = tkr.phi;
  tkr.phi0     = tkr.phi;
  tkr.phi1     = tkr.phi;
  
  % METHODS
  tkr.Heavi      = @Heavi;
  tkr.delta      = @delta;
  tkr.get_means  = @compute_means;
  tkr.get_area   = @compute_area;
  tkr.update_phi = @update_phi;
  tkr.display    = @displayLevelSets; 
  tkr.get_center = @get_center;
  tkr.compensate = @apply_g_compensation;
  
  if KOpts.bCalibrateDelay % forget the tracking, send fixed signal
    tkr.get_center = @get_center_DelayCalib;
  end
  
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
  
  function xyF = get_center_DelayCalib()
    global TKR;    assert( KOpts.bCalibrateDelay );
    frame_idx  = TKR.curr_Nframe;    f_idx_check=tkr.curr_Nframe;
    xyF        = [tkr.img_size(2)/2 ; tkr.img_size(1)/2];
    if mod(frame_idx,20) == 0
      xyF(1) = xyF(1)-120;
    elseif mod(frame_idx,30)==0
      xyF(1) = xyF(1)+120;
    end
  end
  
  function [xyF] = get_center( )
    [yc xc] = find( tkr.phi >= 0 ); 
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
    
    if ndims(Img) == 3
      Img = rgb2gray(Img);
    end
    if nargin < 3
      t_0to1 = 0;
    end
    phi          = tkr.phi;
    bUseUcontrol = tkr.bUseUcontrol; 
    
    if nargin < 2 
      U = 0 * phi;
    end
    %TKR.U=U;
    
    % % Extract ROI (mostly for speed)
    [roi_i, roi_j]= find( phi > -dX*10 );
    
    left  = min(roi_j(:));
    right = max(roi_j(:));
    top   = min(roi_i(:)); 
    bottom= max(roi_i(:));
    Img   = Img(top:bottom,left:right);
    phi   = phi(top:bottom,left:right);
    U     = U(top:bottom,left:right);
    
    % % Image-Dependent part
    [mu_i, mu_o]  = compute_means(Img,phi);
    g_alpha= -(Img - mu_i).^2 + (Img - mu_o).^2;

    % % Regularization part 
    [kappa_phi, phix, phiy,phix2,phiy2] = kappa( phi );
    
    iMinMax   = TKR.imgMinMax; 
    Gmax      = max(  (iMinMax(2)-Img).^2 , (iMinMax(1)-Img).^2 );
      
    isItReallyGmax = ( sum( abs(g_alpha(:)) > Gmax(:) ) < 1 ); 
    assert( isItReallyGmax );
      
    if bUseUcontrol  
      tkr.f_of_U = get_fU( phi, U, Img, g_alpha, phix./sqrt(phix2+phiy2+1e-12), ...
                                  -phiy./sqrt(phix2+phiy2+1e-12) ); 
    else
      tkr.f_of_U = 0*phi;
    end
    
    
    lambda = tkr.lambda;
    dphi  = delta(phi) .* ( g_alpha + tkr.f_of_U + lambda * kappa_phi) ;
    dt0   = 0.75;
    dt_a  = dt0 / max(abs(dphi(:)));  
    phi   = phi + dt_a * dphi;
    
    tkr.phi(top:bottom,left:right) = phi; % must force copy here, not return it
    top = max(1,top-8); bottom = min(tkr.img_size(1),bottom+8);
    left= max(1,left-8); right = min(tkr.img_size(2),right+8);
    redistIters = 3; % 1 or 2, negligible effect for speed
    tkr.phi(top:bottom,left:right) = reinitializeLevelSetFunction( ... 
                   tkr.phi(top:bottom,left:right), 2, dX,redistIters, 1, 1, true() );
%     TKR.xroi = roi_j;
%     TKR.yroi = roi_i;
  end

  
  function displayLevelSets(img0)
    global TKR; % d'oh 
    
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
      phi0     = TKR.phi0;
      imgg( abs( phi0 ) < phi_show_thresh ) = 0;
      imgb( abs( phi0 ) < phi_show_thresh ) = 0;
      imgr( abs( phi0 ) < phi_show_thresh) = (imgr( abs( phi0 ) < phi_show_thresh) .* ...
          abs( phi0(abs(phi0) < phi_show_thresh ) )/phi_show_thresh  + ...
          1.5 * (phi_show_thresh - abs( phi0(abs(phi0) < phi_show_thresh ) ) )/phi_show_thresh );
    end
    
    % draw contour in green for tracker
    imgr( abs( phi ) < phi_show_thresh ) = 0;
    imgb( abs( phi ) < phi_show_thresh ) = 0;
    imgg( abs( phi ) < phi_show_thresh) = (imgg( abs( phi ) < phi_show_thresh) .* ...
      abs( phi(abs(phi) < phi_show_thresh ) )/phi_show_thresh  + ...
      1.5 * (phi_show_thresh - abs( phi(abs(phi) < phi_show_thresh ) ) )/phi_show_thresh );
      
    img_show(:,:,1) = imgr; img_show(:,:,2) = imgg; img_show(:,:,3) = imgb;
    img_show(img_show>1)=1; img_show(img_show<0)=0;
    
    TKR.img1(:)     = TKR.img0(:);
    TKR.img0(:)     = img0(:);
    TKR.img_show(:) = img_show(:);
    
    TKR.phi(:)      = TKR.phi1(:);
    TKR.phi1(:)     = TKR.phi0(:);
    TKR.phi0(:)     = phi(:);
    
    
    g_f2f = TKR.g_f2f;
    zf2f=real(logm(g_f2f));
    roll_ang=( zf2f(2,1) ) * 180/pi; 
    pitch_ang=( zf2f(3,2) ) * 180/pi; 
    yaw_ang =( zf2f(3,1) ) * 180/pi; 
    fprintf('roll_ang=%4.4f, yaw_ang=%4.4f, pitch_ang=%4.4f\n',roll_ang,yaw_ang,pitch_ang);
    if abs(roll_ang)  > 2.0
      fprintf('\nLarge Roll!! =%4.4f ',roll_ang);
      centroidError = norm( TKR.xyF - [320;240], 2);
      if abs(roll_ang) > 2.0 && centroidError > 200.0
        fprintf(', and Large |xyF|!! =%4.4f \n',centroidError);
        fprintf(' *** \n');
      end
      fprintf('\n');
    end
    if abs(yaw_ang)  > 1.0
      fprintf('\nLarge Yaw!! =%4.4f \n',yaw_ang);
      breakhere=1;
    end
    
    A_current = TKR.get_area( phi );
    TKR.Area  = A_current;
    fprintf(', A=%4.2f   , done display() ...',A_current);
  end

  function [g_f2f] = apply_g_compensation(  )

    global TKR; 
   
    tkr.g_f2f = TKR.g_f2f;
    tkr.f     = TKR.f;
    fprintf('| got f = %4.4f |, ',tkr.f);
    
    % argh, why !? 
    m            = tkr.img_size(1);
    n            = tkr.img_size(2);  
    [xx yy]      = meshgrid(linspace(1,n,n),linspace(1,m,m));
    tkr.xx = xx; 
    tkr.yy = yy;
    
    
    
    %xx = tkr.xx;
    %yy = tkr.yy;
    
    % Can't warp for large translation if this is off, low error margin
    % Interesting concept: use particle filter, use multiple z0 to sample
    z0           = -100.0; % !!!! Crucial source of uncertainty
    xx           = -(xx - (n)/2) / tkr.f * z0;
    yy           =  (yy - (m)/2) / tkr.f * z0;
    tkr.xx       = xx;
    tkr.yy       = yy;
    
    [g_ctrl, TKR] = solve_gkC( TKR, KOpts );
    %fuck          = TKR.compare_w0_wX'
    g_f2f         = tkr.g_f2f;
    
    if norm( g_f2f - eye(4,4),'fro') < 1e-6 
      fprintf('Not compensating, seems like g == identity\n'); %nothing to do, frame to frame is identity
    end
    f            = tkr.f;
    assert( (10 < f) && (1e5 > f ) );
    
    % apply the affine warp to tkr.phi 
    u0 = xx;
    v0 = yy;
    
    TKR.g_ctrl = g_ctrl;
    
    %expm(  real(logm( TKR.g_f2f ) - logm( TKR.g_ctrl ) ) );
    g_f2fb     = TKR.g_f2f*(TKR.g_ctrl)^(-1);
    
    TKR.g_k_d  = g_f2fb;
    
    g_comp =  (g_f2fb^-1); % *
    uv     = g_comp * [ u0(:)'; v0(:)'; z0 * ones(1,numel(v0)); ones(1,numel(v0)) ];
    
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
    TKR.xx  = xx;
    TKR.yy  = yy;
    
    fprintf('pixel shift in contour compensation, dx=%3.3g, dy=%3.3g\n',dx,dy);
    
    

  end    



  function A = compute_area(phi)
    A = trapz(trapz(Heavi( phi ) ) );
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

