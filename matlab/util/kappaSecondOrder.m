function [K ngrad dx dy] = kappaSecondOrder(phi, dX, bCheckPoincareIneq)
% input : MxN matrix phi, spatial step dX
% 
% output: kappa, |grad phi|_2, phi_x, phi_y
% throws error if poincare inequality is violated (indicating the data is extremely noisy
% and taking derivatives didn't make sense to begin with)
% peter karasev feb 2012

  if( nargin < 2 )
    dX = 1.0;
  end
  if( nargin < 3 )
    bCheckPoincareIneq = false();
  end
   
  uw_1   = upwindTerms( phi );
  dx     = uw_1.IxCenterDiff / dX;
  dy     = uw_1.IyCenterDiff / dX;
  ngrad  = ( (uw_1.IxCenterDiff).^2 + (uw_1.IyCenterDiff).^2 ).^(1/2) ;
  
  uw_2x  = upwindTerms( uw_1.IxCenterDiff );
  uw_2y  = upwindTerms( uw_1.IyCenterDiff );
  
  dxx    = uw_2x.IxCenterDiff /dX^2;
  dyy    = uw_2y.IyCenterDiff /dX^2;
  dxy    = (uw_2x.IyCenterDiff + uw_2y.IxCenterDiff)*0.5/dX^2;
  
  assert( norm( dxy(:)  - uw_2x.IyCenterDiff(:)/dX^2 ) <= 1e-3*norm(dxy(:)) );
  
  % derivatives
  
  dx2 = dx.^2;    dy2 = dy.^2;
  
  K = (dx2.*dyy + dy2.*dxx - 2*dx.*dy.*dxy)./(dx2 + dy2 + eps);
  
  % Note: Poincare in L1, for phi with compact support, says:
  % sum(sum( abs(phi) ) )  \leq  (d/2) sum(sum( sqrt(dx2 + dy2) ) )
  if( bCheckPoincareIneq )
    poincare_RHS = trapz(trapz( sqrt(dx2 + dy2)*dX^2 ) ) * (ndims(phi)/2);
    poincare_LHS = trapz(trapz( abs(phi) *dX^2 ) );
    assert( poincare_LHS <= poincare_RHS );
  end
  fprintf('');
  
end
