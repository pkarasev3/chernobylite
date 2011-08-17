function [wx wy absw dt] =...
    computeOMTOptFlow_pklasso( I1, I2,alphaNorm,spgl1_opts)
%by PK
% Compute OMT optical flow with L1 regularizer for sparsity.
% must have custom spgl1-1.7 in your path.
% "custom": PK made some modifications for speed/ease_of_use/less_display
%
% strongly suggested before use: go to the spgl1-1.7 dir and do "spgsetup"
% to compile the fast mex/c++ version of code... GCC 4.4.x ubuntu 10.04
% works, not sure about windows compiler ...
%
%INPUT:
% - I1, I2: images (array/pathToFile, color/gray, uint8/double)
% - alphaNorm: (optional) double, alpha = alphaNorm * norm(I2-I1),
%               default 0.01
%
%OUTPUT
% - wx: double array with dim of images (flow field in x direction)
% - wy: double array with dim of images (flow field in y direction)
% - absw: double array with dim of images (magnitude of flow field)

%default alphaNorm
if nargin == 2
    alphaNorm = 1e-2;
end

% read images, if necessary
if ischar(I1)
    I1=imread(I1);
end
if ischar(I2)
    I2=imread(I2);
end

% convert to gray scale, if necessary
if( size(I1,3) > 1 )
    I1 = rgb2gray(I1);
end
if( size(I2,3) > 1 )
    I2 = rgb2gray(I2);
end
I1 = double(I1);
I2 = double(I2);

%equalize mass in I1 and I2
sumI1 = sum(I1(:));
sumI2 = sum(I2(:));
I2 = I2/sumI2*sumI1;
%rescale to have mass between 1 and 2
maxI = max([max(I1(:)), max(I2(:))]);
I1 = I1 / maxI + 1;
I2 = I2 / maxI + 1;

%define alpha in OMT functional via norm of I2-I1
alpha = alphaNorm*norm(I2(:)-I1(:));

%spatial and temporal step size
n1 = size(I1,1); n2 = size(I1,2);  nn = n1*n2;
h1 = 1; h2 = 1;
dt = 1; %dt scales output magnitude linearly

% generate differential operators
%==================================
% I_t
dI = 1/dt*(I2(:) - I1(:));

% \div(I u)
e1 = ones(n1,1);  ddx = spdiags(1/(2*h1)*[-e1,e1],[-1,1],n1,n1);
ddx(1,:) = 0; ddx(end,:) = 0; %BC
Dx = kron(speye(n2),ddx);

e2 = ones(n2,1);  ddy = spdiags(1/(2*h2)*[-e2,e2],[-1,1],n2,n2);
ddy(1,:) = 0; ddy(end,:) = 0; %BC
Dy = kron(ddy,speye(n1));

%Imid and DivI
Imid = kron(speye(2),spdiags(I2(:)+I1(:),0,nn,nn));
DivI = 0.5 * [Dx Dy] * Imid;

%setup A and b for lasso
A  = [DivI ; alpha * [Dx  Dy] ]; % gradient smoothing
b  = [-dI ; zeros( numel(dI), 1 ) ];

% initial fit: least-squars with  alpha ||x||_2
Nvars      = n1 * n2 * 2;
x0         = (A'*A + alpha * speye( Nvars ) ) \ (A' * b);
lasso_size =  norm( A*x0 - b ) * (1+alpha);

if( nargin < 4 )
    spgl1_opts = spgSetParms('iterations',500,'verbosity',0,'iter_skip',50);
end
options = spgl1_opts;

if( lasso_size > norm(b) ) % what the ??
    fprintf(' ... something obscure happend , lasso_size < ||b|| should be impossible \n' );
    lasso_size = norm(b);
end
if( norm(b) < 1e-12 )
    x = zeros( Nvars,1);
    fprintf('b is all zeros! \n');
else
    [x,r,g,info]  = spg_bpdn(A, b, lasso_size, options );
    % display_info_spgl1( info,x );
end

wx = x(1:nn);
wy = x(nn+1:end);
wx = reshape(wx, n1, n2);
wy = reshape(wy, n1, n2);
absw = sqrt(reshape(wx.^2 + wy.^2,n1,n2));
  
end

function display_info_spgl1( info,x ) %#ok<DEFNU>

sfigure(1);

subplot(3,1,1); plot( info.xNorm1, 'b.' ); 
title( ['SPGL1 Output, |x|_0 \leq ' num2str(numel(x)) ] );
legend( 'xNorm1');
subplot(3,1,2); plot( info.rNorm2, 'r.' ); 
legend( 'rNorm2');
subplot(3,1,3); plot( info.lambda, 'g.' );
legend( 'lambda' );
pause(.001);
end

  
  
