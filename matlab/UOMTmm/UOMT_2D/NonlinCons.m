function [c,ceq] = NonlinCons(y, param)

m = param.m;
n = param.n;
Dx = param.Dx; % differential operator (matrix)
Dy = param.Dy; % differential operator (matrix)
X = param.X;   % grid coordinates 
Y = param.Y;   % grid coordinates 
mu0 = param.mu0; % fixed image 0
mu1 = param.mu1; % fixed image 1

%get local variables
u_   = y(1:end/3);
v_   = y(end/3+1:2*end/3);
psi_ = y(2*end/3+1:end);
u = reshape(u_, m, n);
v = reshape(v_, m, n);
psi = reshape(psi_, m, n);

%determinant of jacobian w.r.t space
%========================
dxu_ = Dx*u_;
dyu_ = Dy*u_;
dxv_ = Dx*v_;
dyv_ = Dy*v_;
detDu_ = dxu_ .* dyv_ - dxv_ .* dyu_;

mu1_u = interp2( X', Y', mu1, u', v', 'linear', min(mu1(:))); 
mu1_u_ = mu1_u(:);

%mass preservation constraint +
%boundary is mapped to boundary
%=============================

% Have to move these to the augmented lagrangian!
% Use a smoothed step function penalizer!
% Then Trust-Region-Reflective can be used with hessian pattern!
ceq = [detDu_.*mu1_u(:)  - mu0(:) - psi_];

%derivative > 0
c   = -detDu_; 


% Graveyard:
%  (u(1,:)-X(1,:))';
%        (u(end,:)-X(end,:))';
%        v(:,1)-Y(:,1);
%        v(:,end)-Y(:,end)


 
