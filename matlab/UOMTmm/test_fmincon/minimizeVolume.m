clearvars

maxA = 6;
options=optimset('Algorithm','interior-point');

Aeq = [1 1 1];
beq = 6;

A = [0 0 1; 0 1 0];
b = [1;1];

LB = zeros(3,1);
UB = Inf(3,1);

x0 = [1 1 1];
x = fmincon(@(x)MyObjective(x), x0, A, b, Aeq, beq, LB , UB, ...
            @(x) MyNonlinCons(x,maxA), options);