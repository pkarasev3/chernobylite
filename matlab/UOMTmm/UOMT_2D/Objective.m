function [E gradE] = Objective(y, param, bratwurst )

X = param.X;
Y = param.Y;
mu0 = param.mu0;

%get local variables
u_   = y(1:end/3);
v_   = y(end/3+1:2*end/3);
psi_ = y(2*end/3+1:end);

%energy
E = sum( ((u_-X(:)).^2 + (v_ - Y(:)).^2) .*mu0(:) );%+ 1e1*sum(psi.^2);
  if( strcmp('on',param.bGradObjOn) )
    if( nargin < 3 ) % prevent infinite recursion ...
      gradE = gradest( @(y)Objective(y,param,1), y );
    end
    if( nargin < 3 ) % we're in the initial call, not gradest
      fprintf('done making bratwurst. \n');
    end
  end

end
