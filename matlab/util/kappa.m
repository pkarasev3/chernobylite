function [K dx dy dx2 dy2 p dX] = kappa(phi, p, dX)
  if( nargin < 2 )
    dX = 1.0;
  end
  if( nargin == 2 && numel(p) == 1 )
    dX=p;
  elseif nargin == 2
    dX = 1.0;
  end
  use_all_idx = false;
  if( nargin < 2 || numel(p) == 1 )
      p = 1:numel(phi);
      use_all_idx=true;
  end
  

  % *** Performance Note: *** there is much repetition of computation here! Profiler result 
  % indicates we waste massive time here doing safe_sub2ind(), when indices p are often constant 
  % between calls to this function. 
  if ~exist('Dxx','var')
    kappa_helper; 
    %fprintf('Generated KappaHelper; should only happen once!\n');
    assert( 0<exist('Dy','var') );
  end
  
  % derivatives
  dx  = Dx(phi);  dy  = Dy(phi);
  dx2 = dx.^2;    dy2 = dy.^2;
  dxx = Dxx(phi); dyy = Dyy(phi); dxy = Dxy(phi);

  alpha = eps;
  K = (dx2.*dyy + dy2.*dxx - 2*dx.*dy.*dxy)./(dx2 + dy2 + alpha);
  
  % Note: Poincare in L1, for phi with compact support, says:
  % sum(sum( sqrt(dx2 + dy2) ) )  \leq  (d/2) sum(sum( abs(phi) ) )
  
  if use_all_idx
    K = reshape(K,size(phi));
    dx=reshape(dx,size(phi));
    dy=reshape(dy,size(phi));
    dx2=reshape(dx2,size(phi));
    dy2=reshape(dy2,size(phi));
  end
  
end


function ind = safe_sub2ind(sz, rr, cc) %#ok<*DEFNU>
  rr(rr < 1) = 1;
  rr(rr > sz(1)) = sz(1);
  cc(cc < 1) = 1;
  cc(cc > sz(2)) = sz(2);
  ind = sub2ind(sz, rr, cc);
end
