function pts = gen2D( Npoints, Nlines)

if( nargin ~= 2 )
  Npoints = 1000;
  Nlines  = 5;
end

% random number of lines
Nlines   = poissrnd( Nlines );

% random number of points
Npoints  = Npoints + randn(1,1) * Npoints * 0.1;

% point partition (todo- be random as well)
pts_k = round( linspace( 1, Npoints, Nlines ) );

pts = zeros( 2, Npoints);

for k = 1:Nlines
  
  
  
end




end
