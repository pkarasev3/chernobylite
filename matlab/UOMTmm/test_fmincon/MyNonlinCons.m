function [c,ceq] = MyNonlinCons(x, maxA)

c(1) = 2*(x(1)*x(2)+x(1)*x(3)+x(2)*x(3))-maxA;
c(2) = -norm(x)+1e-5;
ceq = [];