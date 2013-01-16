function [D b group tau_mid H xv yv] = setup_matrices_matchv( xhat, yhat, tvals, Kgroups)
% setup measurement matrix D, observation vector b
% and assign groups based on time-values tvals spaced by
% increment tspace (scalar). Also return the midpoints of group time intervals.

npts = numel(xhat); assert( npts == numel(yhat) && npts == numel(tvals) );
if( nargin < 4 )
  Kgroups = round(npts/2);
else
  assert( numel(Kgroups) == 1 );
end

integrator_mtrx = gallery('triw',npts,1,npts)';
H               = integrator_mtrx;

tnext           = [tvals(2:end)' , 2*tvals(end)-tvals(end-1)];
H               = H .* repmat( tnext(:)'- tvals(:)', npts, 1 );
H0              = H;
H2              = H*H*H;
H               = H*H;
H2              = H2(2:end,:);
H               = H(2:end,:); 

D    =  [ H , 0*H ;
        0*H , H  ; 
         H0 , H0*0 ;
         H0*0 , H0 ];

% expecting the data to start at (x,y) as zero, that's the part we drop
assert( abs(yhat(1)) < 1e-1 * mean(abs(yhat)) && abs(xhat(1)) < 1e-1 * mean(abs(xhat)) );

b    =  [ xhat(2:end);  yhat(2:end); H0(2:end,:)\smooth(xhat(2:end),10) ; H0(2:end,:)\smooth(yhat(2:end),10) ];

[num_in_group_k, tau_mid] = hist( tvals, linspace(tvals(1),tvals(end),Kgroups) );
group = zeros( npts, 1 );
kprev = 0;
for k = 1:Kgroups
  group( (kprev+1):(kprev+num_in_group_k(k) ) ) = k;
  kprev = kprev + num_in_group_k(k);
end
assert( min(group(:)) > 0 );

group = [repmat( group, 2, 1 )];

fprintf('');

end
