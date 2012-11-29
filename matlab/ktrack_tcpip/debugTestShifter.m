function debugTestShifter()

load shiftDataX_big;

phi = tkr.phi;
phi1= tkr.phi0;

fx  = xx1-xx ;
fy  = yy1-yy ;
sfigure(2); 
roix=1:16:480; roiy = 1:16:640; 
quiver( xx(roix,roiy),yy(roix,roiy), fx(roix,roiy), fy(roix,roiy),3 );

[~,dX] = get_params();
itrMax = 1000;
  for itr = 1:itrMax
    
    
    % Traverse starting at the transformed phi, 
    %                        towards the "original"
    
    [~, dx, dy , ~, ~] = kappa( phi, dX );

    Nx =  dx ./ sqrt(dx.^2+dy.^2+1e-6);
    Ny = -dy ./ sqrt(dx.^2+dy.^2+1e-6);

    fdotN = -( fx .* Nx + fy .* Ny ) ;

    dphi = fdotN .* delta(phi);
    dt   = 0.9/max(abs(dphi(:)));
    phi  = phi+dt*dphi;
    phi = reinitializeLevelSetFunction( phi, 2, dX,5, 2, 2, true() );
    
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
