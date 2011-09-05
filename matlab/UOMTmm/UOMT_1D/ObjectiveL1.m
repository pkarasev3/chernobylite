function E = ObjectiveL1(y, x, mu0)

%get local variables
u   = y(1:end/2);
psi = y(end/2+1:end);

tau = 1e-2;

%energy
E = sum( (u-x).^2 .*mu0) + tau * sum( abs( psi ) );
