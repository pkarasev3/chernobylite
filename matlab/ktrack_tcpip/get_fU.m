function fU = get_fU( phi, U,phiNx,phiNy)
  global TKR;
  
  [Ux Uy] = gradient(U);
  magGradU= 1e-9 + sqrt( Ux.^2 + Uy.^2 );
  Ux= Ux./magGradU; 
  Uy=-Uy./magGradU;
  nphi_dot_nU = Ux.*phiNx + Uy.*phiNy;
  
  fU = 0.1 + 0.9*abs( nphi_dot_nU ) ;
  
  fprintf('');

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

