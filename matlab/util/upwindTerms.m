function [upwind_terms] = upwindTerms( I, u, v )
% returns a struct with upwind derivative terms in it 
% also gives the 2nd order central difference gradient for I
% I should be N x N
% u,v are optional. leave them out to just get central difference gradient

I = double(I);
if( ndims(I) > 2 )
  I = rgb2gray(I);
end

if( ~exist('u','var') || ~exist('v','var') )
  u = zeros(size(I)); v = zeros(size(I));
end

dummy = zeros( size(I) );
upwind_terms = struct('Ix_plus',dummy,'Ix_minus',dummy, ...
  'Iy_plus', dummy, 'Iy_minus', dummy, 'aplus', dummy, 'aminus', dummy, ...
  'bplus', dummy, 'bminus', dummy );


upwind_terms.aplus   = max( u, dummy );
upwind_terms.aminus  = min( u, dummy );
upwind_terms.bplus   = max( v, dummy );
upwind_terms.bminus  = min( v, dummy );

% % shift operators for *constant boundary* condition
  shiftU = @(M) M([1 1:end-1],:);
  shiftR = @(M) M(:,[2:end end]);
  shiftL = @(M) M(:,[1 1:end-1]);
  shiftD = @(M) M([2:end end],:);
% %

% % second order upwind scheme
upwind_terms.Ix_plus  = -shiftR(shiftR(I)) + 4*shiftR(I) - 3 * I;
upwind_terms.Ix_minus = 3 * I - 4 * shiftL(I) + shiftL( shiftL( I ) );
upwind_terms.Iy_plus  = -shiftU(shiftU(I)) + 4*shiftU(I) - 3 * I;
upwind_terms.Iy_minus = 3 * I - 4 * shiftD(I) + shiftD( shiftD( I ) );
% % 

upwind_terms.I  = I;
upwind_terms.u  = u;
upwind_terms.v  = v;

upwind_terms.IxCenterDiff = -( 1/12*(shiftR(shiftR(I))-shiftL(shiftL(I)))+2/3*(shiftL(I)-shiftR(I))) ;
upwind_terms.IyCenterDiff =  ( 1/12*(shiftU(shiftU(I))-shiftD(shiftD(I)))+2/3*(shiftD(I)-shiftU(I))) ;

end

