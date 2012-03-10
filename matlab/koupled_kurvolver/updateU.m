function [U h_out uk] = updateU( U, phi_star,phi,px,py,img, Umax)
% generate updated U given click-point (px,py) in img
% U_t = h(u_k) + div( D(U,\mathbf{x}) \cdot gradU ), nonlinear diffusion

    Uraw      = U;
    h_out     = 0*Uraw;
    uk        = (phi_star(py,px) > 0)*(0 > phi(py,px) ) - (phi_star(py,px) < 0)*(0 < phi(py,px) );     
    epsilon   = sqrt(2); %1.1;%0.8;
    
    Heavi     = @(z)  1 * (z >= epsilon) + (abs(z) < epsilon).*(1+z/epsilon+1/pi * sin(pi*z/epsilon))/2.0;
    delta     = @(z)  (abs(z) <= epsilon).*(1 + cos(pi*z/epsilon))/(epsilon*2.0); %#ok<NASGU>

    [X Y] = meshgrid(1:size(phi,2),1:size(phi,1));
    dX    = 1/epsilon;
 
    dt = dX/Umax;        
    
    if( abs(uk)>0  )
     
      alpha  = 10.0;
      
      h_of_u = (Umax/4) * uk*( (0.5*alpha) >= sqrt( ( (X - px).^2 + (Y - py).^2 ) ) ) ...
                                      .*exp( -alpha*(img(py,px)-img).^2 );
      h_out  = h_of_u;
      Uraw = Uraw + h_of_u;
      sub_iters = 5;
      smth      = fspecial('gaussian',[5 5],0.5);
      for iters = 1:sub_iters
        [Ux Uy] = gradient(imfilter(U,smth,'replicate'), dX);
        diffusionTermX        = ( (U/Umax).^2 .* Heavi( (U/Umax).^2 - 1 ).*Ux );% smth,'replicate');
        diffusionTermY        = ( (U/Umax).^2 .* Heavi( (U/Umax).^2 - 1 ).*Uy );% smth,'replicate');
        [diffXX, diffXY]      = gradient(diffusionTermX,dX); %#ok<NASGU>
        [diffYX, diffYY]      = gradient(diffusionTermY,dX); %#ok<NASGU>
        %assert( norm(diffYX(:) - diffXY(:) ) < Umax*max([norm(diffYX(:)),norm(diffYX(:))]) );
       
        dU    = dt*(diffXX + diffYY) + h_of_u;
        U     = U + dU;
        
        h_of_u = 0*h_of_u; % only done once
      end
    end
end
