function [u v psi HESSIAN] = runUOMT2D(mu0, mu1, interval)
HESSIAN  = [];
param    = []; 
small_sz = 128; % small size to use if we use lazy arg inputs
bGetGradAndHessian = 1; % whether we get grad and hess from fmincon
bUseBratwurst = false;
if(  bUseBratwurst )
  addpath('~/source/tann-lab/trunk/matlab_utils/NumericalDerivativesLib/')
  param.theMethod  = 'active-set'; % 'sqp';
  param.bGradObjOn = 'on';
else
  param.theMethod  = 'sqp';
  param.bGradObjOn = 'off';
end
uomt_info = []; % structure that gets info appended into it for save

%Default args
if nargin == 0
    mu0 = 'diffuseCirc';
    mu1 = small_sz;
    interval = [1 10];
elseif nargin == 1
    mu1 = small_sz;
    interval = [1 10];
elseif nargin == 2
    interval = [1 10];
end

minVal = interval(1);
maxVal = interval(2);
    
%LOAD IMAGES and EQUALIZE
%------------------------
if isnumeric(mu1) && ischar(mu0)
    %synthetic examples
    type = mu0;
    n = mu1;
    [mu0 mu1] = syntheticImagesUnbalanced(type, n);
else
    % read images, if necessary
    if ischar(mu0)
        mu0=imread(mu0);
    end
    if ischar(mu1)
        mu1=imread(mu1);
    end
    
    % convert to gray scale, if necessary
    if( size(mu0,3) > 1 )
        mu0 = rgb2gray(mu0);
    end
    if( size(mu1,3) > 1 )
        mu1 = rgb2gray(mu1);
    end
    mu0 = double(mu0);
    mu1 = double(mu1); 
end
fprintf('The size of stuff is: \n ');
disp( num2str( size( mu0 ) ) );

%rescale to have mass between minVal and maxVal
maxMu = max([max(mu0(:)), max(mu1(:))]);
minMu = min([min(mu0(:)), min(mu1(:))]);
mu0 = minVal+(mu0-minMu)*(maxVal-minVal)/(maxMu-minMu);
mu1 = minVal+(mu1-minMu)*(maxVal-minVal)/(maxMu-minMu);

%parameters
[m, n] = size(mu0);
x = 1:m;
y = 1:n;
[X Y] = meshgrid(x,y);
X = X'; Y = Y';

%compute space derivative (center)
Dxp = 1/2*spdiags([-ones(m,1),ones(m,1)], [-1 1], m, m);
Dxp(1,1:3) = 1/2*[-3 4 -1];
Dxp(end,end-2:end) = 1/2*[1 -4 3];
Dx = kron(speye(n), Dxp);

Dy = 1/2*spdiags([-ones(m*n,1),ones(m*n,1)], [-m m], m*n, m*n);
Dy(1:m, 1:3*m) = 1/2*[ -3*speye(m), 4*speye(m), -1*speye(m)];
Dy(end-m+1:end, end-3*m+1:end) =1/2* [1*speye(m), -4*speye(m), 3*speye(m)];

%mass difference
deltaM = sum(mu1(:)) - sum(mu0(:));
%linear constraints on u and psi
if deltaM == 0
    error('masses are not unbalanced!')
elseif deltaM > 0
    LB = [ones(m*n,1); ones(m*n,1); zeros(m*n,1)];
    UB = [m*ones(m*n,1); n*ones(m*n,1); Inf(m*n,1)]; 
else
    LB = [ones(m*n,1); ones(m*n,1); -Inf(m*n,1)];
    UB = [m*ones(m*n,1); n*ones(m*n,1); zeros(m*n,1)];
end

%initial guess
y0 = [X(:); Y(:); zeros(m*n,1)];

%setup param struct
param.m = m;
param.n = n;
param.mu0 = mu0;
param.mu1 = mu1;
param.Dx  = Dx;
param.Dy  = Dy;
param.X = X;
param.Y = Y;

%minimize:
%==========
options=optimset('Algorithm',param.theMethod,'GradObj',param.bGradObjOn,...
                 'display','iter-detailed','TolFun',1e-3,'TolX',1e-3, ...
                 'TolCon',1e-4, 'MaxFunEvals',Inf,'UseParallel','always');

%H = hessian(@(y)Objective(y, param),y0+randn(numel(y0),1)*1e-1 );
if( bGetGradAndHessian )
[y,fval,exitflag,output,lambda,grad,hessOut]  = fmincon(@(y)Objective(y, param), y0, [], [], [], [], LB , UB, ...
            @(y) NonlinCons(y, param), options);   %#ok<ASGLU>
  uomt_info.fval = fval; clear fval;
  uomt_info.grad = grad; clear grad;
  uomt_info.hessian = hessOut; clear hessOut;
else
  [y]  = fmincon(@(y)Objective(y, param), y0, [], [], [], [], LB , UB, ...
              @(y) NonlinCons(y, param), options);  
end
%return result
u_   = y(1:end/3);
v_   = y(end/3+1:2*end/3);
psi_ = y(2*end/3+1:end);
u = reshape(u_, m, n);
v = reshape(v_, m, n);
psi = reshape(psi_, m, n);

%save result
datenow = datestr(now,31);
datenow = strrep(datenow,':','-');
datenow = strrep(datenow,' ','_');
filename = ['result_' datenow]; 
save(filename);
save('last.mat');

analyzeResult(filename);

end








