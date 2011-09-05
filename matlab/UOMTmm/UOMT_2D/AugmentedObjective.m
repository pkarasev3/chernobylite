function E = AugmentedObjective(y,Leq,Lineq,param)
% The augmented lagrangiant L = f + L1 + L2
% Leq: lambda lagranage mults for equality constraints
% Lineeq: lambda lagranage mults for inequality constraints

X = param.X;
Y = param.Y;
mu0 = param.mu0;

%get local variables
u_   = y(1:end/3);
v_   = y(end/3+1:2*end/3);
psi_ = y(2*end/3+1:end);

%energy
E1 = sum( ((u_-X(:)).^2 + (v_ - Y(:)).^2) .*mu0(:) );%+ 1e1*sum(psi.^2);

[c,ceq] = NonlinCons(y, param);

E2 = sum( Leq .* c  + Lineq .* Ceq );

E = E1 + E2;

end
