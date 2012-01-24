function [xhat yhat tvals foo ax ay az] = read_sensor_data( filename )

foo=csvread(filename,1,0);  %#ok<*NOPTS> % "WGS84"

lat=foo(:,1) ;
lon=foo(:,2) ;

tvals = (foo(1:end,end) - foo(1,end) ) * 1e-9;

xhat = smooth(lon,20);
yhat = smooth(lat,20);
tvals= smooth(tvals,20);
xhat = (xhat-xhat(1));
yhat = (yhat-yhat(1));
tvals= tvals-tvals(1);

ax   = foo(:,5);
ay   = foo(:,6);
az   = foo(:,7);

netA = sqrt( ax.^2+ay.^2 +  az.^2 );
figure(1); plot( tvals, netA ); 
xlabel('Time (seconds)','FontSize',18); 
ylabel('|Acceleration| (meters/seconds^2)','FontSize',18);

fprintf('');

end
