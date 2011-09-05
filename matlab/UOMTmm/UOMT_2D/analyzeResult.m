function [] = analyzeResult(filenameIn)

%initialize psi
psi=0;

if strcmp(filenameIn, 'last')
    a = load('last.mat');
    load(a.filename);
else
    load(filenameIn);
end


minVal = min(min(mu0(:)), min(mu1(:)));
maxVal = max(max(mu0(:)), max(mu1(:)));
 
disp('RESULT:')
disp('-------')
%check constraints
dxu_ = Dx*u_;
dyu_ = Dy*u_;
dxv_ = Dx*v_;
dyv_ = Dy*v_;
detDu_ = dxu_ .* dyv_ - dxv_ .* dyu_;

mu1_u = interp2( X', Y', mu1, u', v', 'linear',min(mu1(:)));  
mu1_u(isnan(mu1_u)) = mu1(isnan(mu1_u));

ceq = detDu_.*mu1_u(:)  - mu0(:) - psi_;
% detDu>0?
if any(detDu_ <=0)
    disp('Warning : some du/dx <= 0!');
    disp([ 'Relative Magnitude: '  ... 
               num2str(  abs(sum( detDu_(detDu_<0) ))/max(detDu_) ) ] );
else
    disp('GOOD: du/dx > 0!')
end
% u,v in bounds and boundary mapped to boundary?
if max( u_ ) > n+1e-3 || min(u_) < 1-1e-3 || ... 
     max( v_ ) > m+1e-3 || min(v_) < 1-1e-3 
    disp('BAD : u out of bounds!')
else
    disp('GOOD: u is in bounds!')
end
if max(abs((u(1,:) - X(1,:) ) ) > 1e-3 ) || ... 
    max(abs((u(end,:) - X(end,:) ) ) > 1e-3 ) ||...
      max(abs((v(:,1) - Y(:,1) ) ) > 1e-3 ) || ... 
       max(abs((v(:,end) - Y(:,end) ) ) > 1e-3 ) %#ok<*NODEF>
    disp('BAD : boundary not mapped to boundary!')
else
    disp('GOOD: u, v map boundary to boundary!')
end
% psi right sign?
if deltaM>0
    if any(psi_ <0)
        tol = 1e-12;
        if( min( psi_ ) < -tol )
          disp('BAD : psi < eps!')
        else
          disp( ['Acceptable:  ' num2str(-tol) ' < psi < 0 !' ]);
        end
    else
        disp('GOOD: psi > 0!')
    end
else
    if any(psi_ >0)
        tol = 1e-12;
        if( min( psi_ ) > tol )
          disp('BAD : psi < eps!')
        else
          disp( ['Acceptable:  ' '0 < psi < ' num2str(-tol)  ' !' ]);
        end
    else
        disp('GOOD: psi < 0!')
    end
end

%max equality constraint
disp(['MAX ERROR for equality cons:   ' num2str(max(abs(ceq)))])
%sum equality cons error
disp(['CUMULATIVE ERROR for equ cons: ' num2str(sum(abs(ceq)))])
%functional value
disp(['Functional Value           : ' num2str(sum( ((u_-X(:)).^2 + (v_ - Y(:)).^2) .*mu0(:) ))])
%mass difference
disp(['Mass difference            : ' num2str(deltaM)]);
%total mass creation
disp(['Mass created by psi        : ' num2str(sum(psi_))])

%plots
mu1_mapped =  reshape(detDu_.*mu1_u(:),m,n);
sfigure(1);
subplot(3,2,1);
imshow(uint8(mu0), [minVal maxVal]);
title('mu0')
subplot(3,2,3); 
imshow(uint8(mu1), [minVal maxVal]);
title('mu1')
subplot(3,2,5);
imshow(uint8(mu1_mapped), [minVal maxVal]);
title('mu1 mapped')
subplot(3,2,6);
minPsi = min(abs(psi_));
maxPsi = max(abs(psi_));
imshow(uint8(abs(psi)), [minPsi maxPsi]);
title('\psi_{L2}: source term');
subplot(3,2,2);
imshow(uint8(mu1_mapped-psi), [minVal maxVal]);
title('mu1 mapped - \psi_{L2} (should equal mu0)');

sfigure(2);
subplot(3,3,1);
imagesc(u-X);
title('distance in x')
subplot(3,3,5);
imagesc(v-Y);
title('distance in y')
subplot(3,3,9);
imagesc(sqrt((u-X).^2+(v-Y).^2));
title('norm of displacement')

sfigure(3);
imagesc(reshape(detDu_,m,n));
title('det(Du)')


% Graveyard
% u_f = floor(u_);
% u_c  = ceil(u_);
% u_w = u_ - u_f;
% v_f = floor(v_);
% v_c  = ceil(v_);
% v_w = v_ - v_f;
% 
% iff = sub2ind([m n], u_f, v_f);
% ifc = sub2ind([m n], u_f, v_c);
% icf = sub2ind([m n], u_c, v_f);
% icc = sub2ind([m n], u_c, v_c);
% mu1_u = (1-u_w).*(1-v_w) .* mu1(iff)+...
%         u_w.*(1-v_w)     .* mu1(icf)+...
%         (1-u_w).*v_w     .* mu1(ifc)+...
%         u_w.*v_w         .* mu1(icc);
