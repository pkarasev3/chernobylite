function K = kappa(phi, p, dX)
  
  if( nargin < 3 )
    dX = 1.0;
  end
  [rr cc] = ind2sub(size(phi), p);

  % *** Performance Note: *** there is much repetition of computation here! Profiler result 
  % indicates we waste massive time here doing safe_sub2ind(), when indices p are often constant 
  % between calls to this function. 
  
  % shift operations
  shiftD = @(M) M(safe_sub2ind(size(phi), rr-1, cc));
  shiftU = @(M) M(safe_sub2ind(size(phi), rr+1, cc));
  shiftR = @(M) M(safe_sub2ind(size(phi), rr,   cc-1));
  shiftL = @(M) M(safe_sub2ind(size(phi), rr,   cc+1));
  shiftUL = @(M) M(safe_sub2ind(size(phi), rr-1, cc-1));
  shiftUR = @(M) M(safe_sub2ind(size(phi), rr-1, cc+1));
  shiftDL = @(M) M(safe_sub2ind(size(phi), rr+1, cc-1));
  shiftDR = @(M) M(safe_sub2ind(size(phi), rr+1, cc+1));

  % derivative operations
  Dx  = @(M) (shiftL(M) - shiftR(M))/(2*dX);
  Dy  = @(M) (shiftU(M) - shiftD(M))/(2*dX);
  Dxx = @(M) (shiftL(M) - 2*M(p) + shiftR(M))/(dX^2);
  Dyy = @(M) (shiftU(M) - 2*M(p) + shiftD(M))/(dX^2);
  Dxy = @(M) (shiftUL(M) + shiftDR(M) - shiftUR(M) - shiftDL(M))/(4*dX^2);
  
  % derivatives
  dx  = Dx(phi);  dy  = Dy(phi);
  dx2 = dx.^2;    dy2 = dy.^2;
  dxx = Dxx(phi); dyy = Dyy(phi); dxy = Dxy(phi);

  alpha = eps;
  K = (dx2.*dyy + dy2.*dxx - 2*dx.*dy.*dxy)./(dx2 + dy2 + alpha);
  
  % Note: Poincare in L1, for phi with compact support, says:
  % sum(sum( sqrt(dx2 + dy2) ) )  \leq  (d/2) sum(sum( abs(phi) ) )
  
end

function ind = safe_sub2ind(sz, rr, cc)
  rr(rr < 1) = 1;
  rr(rr > sz(1)) = sz(1);
  cc(cc < 1) = 1;
  cc(cc > sz(2)) = sz(2);
  ind = sub2ind(sz, rr, cc);
end
