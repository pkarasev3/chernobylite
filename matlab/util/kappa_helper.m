  [rr cc] = ind2sub(size(phi), p); 
  
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
