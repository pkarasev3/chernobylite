function [c,ceq] = NonlinCons(y, D, x, mu0, mu1)

%get local variables
u   = y(1:end/2);
psi = y(end/2+1:end);

%du/dx
du = D*u;

%mu1(u(x))
u_floor = floor(u);
u_ceil  = ceil(u);
u_weight = u - u_floor;
mu1_u = (1-u_weight).*mu1(u_floor)+ u_weight.* mu1(u_ceil);

%mass preservation constraint
ceq = [du.*mu1_u  - mu0 - psi; 
       u(1)-x(1); 
       u(end)-x(end) ];

%derivative > 0
TolCon = 0; % Not a hack, its the 'TolCon' setting in Options
c   = -du+ TolCon; 





 
