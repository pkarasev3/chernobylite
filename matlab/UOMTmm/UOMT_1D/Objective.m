function E = Objective(y, x, mu0)

%get local variables
u   = y(1:end/2);
psi = y(end/2+1:end);

%energy
E = sum( abs(u-x) .*mu0);%+ 1e1*sum(psi.^2);