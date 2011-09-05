clearvars
close all;

%parameters
n  = 120;
discretize = 'forward';


%space var and space step
x  = (1:n)';
dx = 1;
x0 = 1;
x1 = n;
xc = (n-1)/2;
xd = (n-1)*2/3;

%define mu0
%spline with (mu(a)=mu'(a)=mu(b)=mu'(b) = 0, 1 in the middle)
M0 = [ x0^4 x0^3 x0^2 x0 1;...
      x1^4 x1^3 x1^2 x1 1;...
      4*x0^3 3*x0^2 2*x0 1 0;...
      4*x1^3 3*x1^2 2*x1 1 0;...
      xc^4 xc^3 xc^2 xc 1];
r0 = [1;1;0;0;1.5];
coeffs0 =  M0\r0;
mu0 = coeffs0(1)*x.^4+ coeffs0(2)*x.^3+ coeffs0(3)*x.^2+ coeffs0(4)*x+coeffs0(5);
%define mu1
%spline with (mu(a)=mu'(a)=mu(b)=mu'(b) = 0, 1 at 1/2 and 2/3 of interval)
% M1 = [ x0^5 x0^4 x0^3 x0^2 x0 1;...
%       x1^5  x1^4 x1^3 x1^2 x1 1;...
%       5*x0^4 4*x0^3 3*x0^2 2*x0 1 0;...
%       5*x1^4 4*x1^3 3*x1^2 2*x1 1 0;...
%       xc^5 xc^4 xc^3 xc^2 xc 1;
%       xd^5 xd^4 xd^3 xd^2 xd 1 ];
% r1 = [2;2;0;0;2;1.6];
% coeffs1 =  M1\r1;
% mu1 = coeffs1(1)*x.^5+coeffs1(2)*x.^4+ coeffs1(3)*x.^3+...
%       coeffs1(4)*x.^2+ coeffs1(5)*x+coeffs1(6);
psi_m = ones(n,1);
% psi_m(n/3:n/3+5) = 2;
% psi_m(n/3+10:n/3+15) = 0.5;
psi_a = zeros(n,1);
psi_a(9*n/10:9*n/10+10) = 1/100*[0 1 4 8 10 9 8 4 2 1 0];
psi_a(n/10:n/10+10) = 1/100*[0 1 4 8 10 9 8 4 2 1 0];
psi_a(n/10+10:n/10+20) = - 1/100*[0 1 4 8 10 9 8 4 2 1 0];

mu1_ = mu0.*psi_m+psi_a;
mu1  = mu1_ + smooth(randn(size(mu1_))*1e-3 * norm(mu1_),2);


%plot
sfigure(1)
plot(x,mu0)
hold on
plot(x,mu1,'r')
legend('mu0','mu1')
hold off

%compute space derivative
switch discretize
    case 'center'
        D = 1/(2*dx)*spdiags([-ones(n,1),ones(n,1)], [-1 1], n, n); %center
        D(1,1:3) = 1/(2*dx)*[-3 4 -1];
        D(end,end-2:end) = 1/(2*dx)*[1 -4 3];
    case 'forward'
        D = 1/dx*spdiags([-ones(n,1),ones(n,1)], [0 1], n, n);  %forward
        D(1,1:3) = 1/(2*dx)*[-3 4 -1];
        D(end,end-2:end) = 1/(2*dx)*[1 -4 3];
    case 'backward'
        D = 1/dx*spdiags([-ones(n,1),ones(n,1)], [-1 0], n, n);  %backward
        D(1,1:3) = 1/(2*dx)*[-3 4 -1];
        D(end,end-2:end) = 1/(2*dx)*[1 -4 3];
end

%mass difference
deltaM = sum(mu1) - sum(mu0);
%linear constraints on u and psi
if deltaM == 0
    error('masses are not unbalanced!')
elseif deltaM > 0
    LB = [ones(n,1); zeros(n,1)];
    UB = [n*ones(n,1); Inf(n,1)]; 
else
    LB = [ones(n,1); -Inf(n,1)];
    UB = [n*ones(n,1); zeros(n,1)];
end

%initial guess
y0 = zeros(2*n,1);
%y0 = [zeros(n,1); zeros(n/3,1); 0.7*ones(2*n/3,1)];

%minimize:  solve L2 version fast, initialize for L1
%==========
options=optimset('Algorithm','sqp','display','iter-detailed',... 
                          'TolFun',1e-4,'TolX',1e-4,'TolCon',1e-4, 'MaxFunEvals',24000);
[yL2] = fmincon(@(y)ObjectiveL2(y, x, mu0), y0, [], [], [], [], LB , UB, ...
            @(y) NonlinCons(y, D, x, mu0, mu1), options);
y0  = yL2;          
[y fval flag out lambda] = fmincon(@(y)ObjectiveL1(y, x, mu0), y0, [], [], [], [], LB , UB, ...
            @(y) NonlinCons(y, D, x, mu0, mu1), options);

disp('constraint violation: ' );          
disp( lambda  )

%get result
u   = y(1:end/2);
psi = y(end/2+1:end);

disp('RESULT:')
disp('-------')
%check constraints
du = D*u;
u_floor = floor(u);
u_ceil  = ceil(u);
u_weight = u - u_floor;
mu1_u = (1-u_weight).*mu1(u_floor)+ u_weight.* mu1(u_ceil);
ceq = du.*mu1_u  - mu0 - psi;
% du>0?
if any(du <=0)
    disp('BAD : du/dx <= 0!')
else
    disp('GOOD: du/dx > 0!')
end
% u between 1 and n?
if any(or(u<1, u>n)) || u(1)~= x(1) || u(end)~= x(end)
    disp('BAD : u out of bounds!')
else
    disp('GOOD: u is in bounds!')
end
% psi right sign?
if deltaM>0
    if any(psi <0)
        tol = 1e-12;
        if( min( psi ) < -tol )
          disp('BAD : psi < eps!')
        else
          disp( ['Acceptable:  ' num2str(-tol) ' < psi < 0 !' ]);
        end
    else
        disp('GOOD: psi > 0!')
    end
else
    if any(psi >0)
        disp('BAD : psi > 0!')
    else
        disp('GOOD: psi < 0!')
    end
end

%max equality constraint
disp(['MAX ERROR for equality cons: ' num2str(max(abs(ceq)))])
%functional value
disp(['Functional Value           : ' num2str(sum( (u-x).^2 .*mu0))])
%mass difference
disp(['Mass difference            : ' num2str(deltaM)]);
%total mass creation
disp(['Mass created by psi        : ' num2str(sum(psi))])

%plots
mu1_mapped =  du.*mu1_u;
sfigure(2);
subplot(3,1,1); plot(x,mu0,'b--');
hold on
subplot(3,1,1); plot(x,mu1,'r');
subplot(3,1,1); plot(x,mu1_mapped, 'g-.');
legend('mu0','mu1','mu1 mapped');
hold off

sfigure(2);
subplot(3,1,2); plot(x,psi,'b--'); hold on;
subplot(3,1,2); plot(x,yL2(end/2+1:end),'r-.'); hold off;
legend('\psi_{L1}: source term','\psi_{L2} : source term');

sfigure(2);
subplot(3,1,3); plot(x,(u-x));
legend('conservative displacement');

sfigure(3);
subplot(2,1,1); plot(x,u);
legend('mapping u');

sfigure(3); subplot(2,1,2);
plot(x,du);
legend('du/dx');







