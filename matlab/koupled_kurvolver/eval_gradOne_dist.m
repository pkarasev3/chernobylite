 function [Fval] = eval_gradOne_dist( phi )                             
   
  p       = 1:numel(phi(:));
  [rr cc] = ind2sub(size(phi), p);

  % shift operations
  shiftD = @(M) M(safe_sub2ind(size(phi), rr-1, cc));
  shiftU = @(M) M(safe_sub2ind(size(phi), rr+1, cc));
  shiftR = @(M) M(safe_sub2ind(size(phi), rr,   cc-1));
  shiftL = @(M) M(safe_sub2ind(size(phi), rr,   cc+1));
  
  % derivative operations
  Dx  = @(M) (shiftL(M) - shiftR(M))/2;
  Dy  = @(M) (shiftU(M) - shiftD(M))/2;
  
  % derivatives
  phix  = reshape(Dx(phi),size(phi));
  phiy  = reshape(Dy(phi),size(phi));
  
  % TODO: evaluate boundary conditions properly! don't force constant value ...
    
  Fval = trapz(trapz( (phix.^2+phiy.^2-1).^2 ) );
    
 end
  

function ind = safe_sub2ind(sz, rr, cc)
  rr(rr < 1) = 1;
  rr(rr > sz(1)) = sz(1);
  cc(cc < 1) = 1;
  cc(cc > sz(2)) = sz(2);
  ind = sub2ind(sz, rr, cc);
end
