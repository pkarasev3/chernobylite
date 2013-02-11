function fU = get_fU( phi, U, Img, phiNx, phiNy)
  global TKR;
  
%   [Ux Uy] = gradient(U);
%   magGradU= 1e-9 + sqrt( Ux.^2 + Uy.^2 );
%   Ux= Ux./magGradU; 
%   Uy=-Uy./magGradU;
%   nphi_dot_nU = Ux.*phiNx + Uy.*phiNy;
%   
%   fU = 0.1 + 0.9*abs( nphi_dot_nU ) ;
  fU = fU_area_regulator( phi, Img, phiNx, phiNy );


  fprintf('');

end

function f = fU_area_regulator( phi, Img, g_alpha, phiNx, phiNy )
% normal is inward (phi>0 interior)
% area increase => -N dot N
% area decrease =>  N dot N
  global TKR;
  global KOpts;
  A_current = TKR.get_area( phi );
  A_min     = pi * 15^2;
  A_max     = pi * 17^2;
  
  iMinMax   = TKR.imgMinMax; 
  Gmax      = max(  (iMinMax(2)-Img).^2, (iMinMax(1)-Img).^2 );
  
  [yc xc] = find( phi >= -2 );
  yc      = [yc(:); round(size(phi,1)/2)];
  xc      = [xc(:); round(size(phi,2)/2)];
  
  xyC     = [mean(xc(:));mean(yc(:))];
  %xyC     = xyC + 0.1*(TKR.xyF - TKR.xyF_prev);
  
  distToC    = 0*phi;
  distToC(sub2ind(size(phi),yc,xc)) = ... 
            ( abs(yc(:)-xyC(2)) + abs(xc(:)-xyC(1)) );
  
  f = 0*phi; % sqrt(A_current/pi) !? 
  if A_current < A_min 
    f         =  Gmax;
  elseif A_current > A_max
    f         = -Gmax;
  end
  rmax = 40.0;  
  f(distToC>rmax) = f(distToC>rmax) - ...
                       Gmax(distToC>rmax).* ...
                        (distToC(distToC>rmax)-rmax)*0.5;
  rmin = 3.0;  
  f(distToC<rmin) = f(distToC<rmin) + Gmax(distToC<rmin);

end




function f = fU_old_version( phi, U )
    global TKR;    
    U = imfilter( (-1 + 2*(U>0)), fspecial('gaussian',[5 5],3),'replicate');
    [U_i, U_o] = TKR.get_means(U,phi);
    % f: |U-U_i| < uMin  maps to 0
    % f: |U-U_i| >=1  maps to -G 
    uMin   =        0.1;
    xi_sqr = abs(U-U_i); 
    (xi_sqr > uMin).*( min(xi_sqr - uMin, 1.0 ) );
end

% % incremental warp (removed currently) 
% 
% fx   = TKR.fx(top:bottom,left:right); 
% fy   = TKR.fy(top:bottom,left:right);
% 
% if KOpts.incremental_warp
%       Nx =  phix ./ sqrt(phix.^2+phiy.^2+1e-6);
%       Ny = -phiy ./ sqrt(phix.^2+phiy.^2+1e-6);
%       fgain = 50.0 / KOpts.contour_iters;
%       fdotN = -( fx .* Nx + fy .* Ny )*fgain ;
%       %fdotN =  fgain / (1+max(abs(fdotN(:))));
%     else
%       fdotN = 0*phi;
%     end
%f_of_U     = (tkr.f_of_U).*(-g_alpha);
