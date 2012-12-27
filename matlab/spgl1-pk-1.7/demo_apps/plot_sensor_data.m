 % drive 1 
foo1=csvread('pk_to_starbux_loop_raw_sensor_data.txt',1,0);  %#ok<*NOPTS> % "WGS84"


foo = foo1;

lat=foo(:,1) ;
lon=foo(:,2) ;
alt=foo(:,3) ;

tt = (foo(1:end,end) - foo(1,end) ) * 1e-9;

plot( lon-lon(1), lat-lat(1),'r.'); xlabel('longitude'); ylabel('latitude'); axis equal;

earth_flatness = 1/298.257223563
earth_equa_rad = 6378137.0000 
earth_pole_rad = 6356752.3142 
a      = earth_equa_rad;
b      = earth_pole_rad;
earth_angular_eccentricity = acos(b/a)

alpha  = earth_angular_eccentricity;

phi    = lat * pi/180;
lambda = lon * pi/180;
h      = alt;

N =  a ./ sqrt( 1 - (sin(phi)*sin(alpha)).^2 );
x = (N+h).*cos(phi).*cos(lambda);
y = (N+h).*cos(phi).*sin(lambda);
z = (cos(alpha)^2*N+h).*sin(phi);

% TODO: weight the distance accumulation by accelerometer norm values

x = smooth(x,256);
y = smooth(y,256);
z = smooth(z,256);
dx = [0 ; x(2:end)-x(1:end-1)];
dy = [0 ; y(2:end)-y(1:end-1)];
dz = [0 ; z(2:end)-z(1:end-1)];
path_length = sum( sqrt( dx.^2+dy.^2+dz.^2 ) );

title(['total distance = ' num2str(path_length)  'meters' ] );

