function [ngrad dx dy] = gradSecondOrder(phi, dX, bCheckPoincareIneq)
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
  
  % derivatives
  dx     = uw_1.IxCenterDiff / (2*dX);
  dy     = uw_1.IyCenterDiff / (2*dX);
  ngrad  = ( (uw_1.IxCenterDiff).^2 + (uw_1.IyCenterDiff).^2 ).^(1/2) ;
  
  % Note: Poincare in L1, for phi with compact support, says:
  % sum(sum( abs(phi) ) )  \leq  (d/2) sum(sum( sqrt(dx2 + dy2) ) )
  if( bCheckPoincareIneq )
    poincare_RHS = trapz(trapz( ngrad*dX^2 ) ) * 2/(ndims(phi));
    poincare_LHS = trapz(trapz( abs(phi) * dX^2 ) );
    assert( poincare_LHS <= poincare_RHS );
  end
  fprintf('');
  
end
