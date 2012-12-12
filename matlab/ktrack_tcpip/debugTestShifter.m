function debugTestShifter()

load shiftDataX_big;

phi = TKR.phi0;
phi1= TKR.phi1;
    
fx  = TKR.fx ;
fy  = TKR.fy ;
xx  = TKR.xx; yy = TKR.yy;
sfigure(2); 
roix=1:32:480; roiy = 1:32:640; 
quiver( xx(roix,roiy),yy(roix,roiy), fx(roix,roiy), fy(roix,roiy),2 );

[epsilon,dX] = get_params();

phi = reinitializeLevelSetFunction( phi, 2, dX,20, 1, 1, true() );
phi1= reinitializeLevelSetFunction( phi1, 2, dX,20, 1, 1, true() );

sfigure(1); 
img_show = repmat( 0.1*TKR.img0, [1 1 3] );
imgg= img_show(:,:,2);
imgr= img_show(:,:,1); imgb=img_show(:,:,3); phi_show_thresh = epsilon/sqrt(2);
 imgr( abs( phi ) < phi_show_thresh ) = 0;
    imgb( abs( phi ) < phi_show_thresh ) = 0;
    imgg( abs( phi ) < phi_show_thresh) = (imgg( abs( phi ) < phi_show_thresh) .* ...
      abs( phi(abs(phi) < phi_show_thresh ) )/phi_show_thresh  + ...
      1.5 * (phi_show_thresh - abs( phi(abs(phi) < phi_show_thresh ) ) )/phi_show_thresh );
    img_show(:,:,1) = imgr; img_show(:,:,2) = imgg; img_show(:,:,3) = imgb;
    img_show(img_show>1)=1; img_show(img_show<0)=0;
    imshow(flipdim(img_show,1));
drawnow; pause(0.05);
matlab2tikz('ktrack_afterRoll_img2_tikz.tex','relativePngPath','./figs','width','6.4cm','height','4.8cm');


    sfigure(3);
img_show2 = repmat( 0.1*TKR.img1, [1 1 3] );
imgg= img_show2(:,:,2);
imgr= img_show2(:,:,1); imgb=img_show2(:,:,3); phi_show_thresh = epsilon/sqrt(2);
 imgr( abs( phi1 ) < phi_show_thresh ) = 0;
    imgb( abs( phi1 ) < phi_show_thresh ) = 0;
    imgg( abs( phi1 ) < phi_show_thresh) = (imgg( abs( phi1 ) < phi_show_thresh) .* ...
      abs( phi1(abs(phi1) < phi_show_thresh ) )/phi_show_thresh  + ...
      1.5 * (phi_show_thresh - abs( phi1(abs(phi1) < phi_show_thresh ) ) )/phi_show_thresh );
    img_show2(:,:,1) = imgr; img_show2(:,:,2) = imgg; img_show2(:,:,3) = imgb;
    img_show2(img_show2>1)=1; img_show2(img_show2<0)=0;
    imshow(flipdim(img_show2,1));    
drawnow; pause(0.05);
matlab2tikz('ktrack_beforeRoll_img1_tikz.tex','relativePngPath','./figs','width','6.4cm','height','4.8cm');
    

%!cp -v ./ktrack_afterRoll*   ~/source/visioncontrol/thesis-pk/figs/
%!cp -v ./ktrack_beforeRoll*   ~/source/visioncontrol/thesis-pk/figs/    

itrMax = 100;
  for itr = 1:itrMax
    
    
    % Traverse starting at the transformed phi, 
    %                        towards the "original"
    
    [~, dx, dy , ~, ~] = kappa( phi, dX );

    % doh, broke 
    %[dx,dy] = computeUpwindDerivatives2D(phi, fx, fy, 2, dX, 3);
    
    Nx =  dx ./ sqrt(dx.^2+dy.^2+1e-6);
    Ny = -dy ./ sqrt(dx.^2+dy.^2+1e-6);

    % should use upwinding here?
    fdotN = -( fx .* Nx + fy .* Ny ) ;

    dphi = fdotN .* delta(phi);
    dt   = 0.9/max(abs(dphi(:)));
    phi  = phi+dt*dphi;
    
    % Note: updated reinit func to return upwind derivs, use them! 
    [phi, phi_x_plus, phi_y_plus, phi_z_plus, ...
              phi_x_minus, phi_y_minus, phi_z_minus] ...
                      = reinitializeLevelSetFunction( phi, 2, dX,5, 2, 2, true() );
    
    sfigure(1); 
    imagesc( 1.0*phi.*(abs(phi>0)) + 1.0*phi1.*(abs(phi1>0)));
    title( sprintf( '%05d iters of %05d', itr, itrMax ) );
    drawnow; 
    
    if( max( phi(:)) < 5 );
      phi = phi+1.0; fprintf('...');
    end
    
  end

end


  function [epsilon,dX] = get_params()
    epsilon   = sqrt(2);
    dX        = 0.5 * 1/sqrt(2);
  end

  function z = Heavi(z)
    [epsilon,dX] = get_params();
    z = 1 * (z >= epsilon) + (abs(z) < epsilon).*(1+z/epsilon+1/pi * sin(pi*z/epsilon))/2.0;
  end

  function z = delta(z)
    [epsilon,dX] = get_params();
    z = (abs(z) < epsilon).*(1 + cos(pi*z/epsilon))/(epsilon*2.0);
  end
