function [pts Nlines] = gen2D( Npoints, Nlines,sigma)

if( nargin ~= 3 )
  Npoints = 1000;
  Nlines  = 2;
  sigma   = 0.02;
end

% random number of lines
Nlines   = 3 + poissrnd( Nlines );
Nlines   = min([ 8, Nlines]);


% random number of points
Npoints  = round( Npoints + randn(1,1) * Npoints * 0.1 );

% point partition (todo- be random as well)
pts_k = round( linspace( 1, Npoints, Nlines ) );

pts = zeros( 2, Npoints);

for k = 1:(Nlines-1)
  
  pt_a = [0;0]; pt_b = [0;0];
  while( norm( pt_a - pt_b ) < 1/sqrt(2) )
    pt_a = randn(2,1) + [k;0] - [Nlines/2;0];
    pt_b = pt_a + randn(2,1) + [k;0] - [Nlines/2;0];
  end
  
  tt   = linspace(0,1,pts_k(k+1)-pts_k(k)+1);
  pts(:,pts_k(k):pts_k(k+1)) = repmat( (pt_b ),1,numel(tt)) .* repmat( tt, 2, 1) + ...
                               repmat( (pt_a),1,numel(tt)) .* repmat( 1-tt, 2, 1);
                             
  pts(:,pts_k(k):pts_k(k+1)) = pts(:,pts_k(k):pts_k(k+1)) + randn(size(pts(:,pts_k(k):pts_k(k+1))))*sigma;
  
end

Nlines = Nlines-1; % correction
fprintf('');


end
