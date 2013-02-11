% A script to generate an input file for the simulator_client program.
% writes a simple CSV array to text, with appropriate extension .xyzuvw

objFileName = 'SomeObjPath';
savePath    = '/u4/home/pkarasev3/source/ktrack/testdata/' %#ok<*NOPTS>
saveToFile  = [savePath, objFileName, '.xyzuvw']                           %#ok<NOPTS>


tt = linspace(0,10,1000); tt = tt(:);

Scenario = 3;

if Scenario == 1
  x  = 10 + 0*tt;
  y  = 0*tt;
  z  = 20 + 0*tt;
  cx = 50+0*tt; cy = 50+0*tt; cz = 20+0*tt;
  ax = x + 5*cos(tt); 
  ay = y; 
  az = z + 5*sin(tt);
elseif Scenario == 2
  % Moving around, dips below horizonish
  x  = 10 * cos(0.9 * pi*tt      ) .* tanh(tt);
  y  = 10 * sin(0.9 * pi*tt * 0.7) .* tanh(tt);
  z  = 20 + 7.0 * sin(2*pi*tt * 1.2) .* tanh(tt*3); 
elseif Scenario == 3
  % Moving around, stays above horizonish
  x  = 10 * cos(0.9 * pi*tt      ) .* tanh(tt);
  y  = 10 * sin(0.9 * pi*tt * 0.7) .* tanh(tt);
  z  = 20 + 4.0 * sin(2*pi*tt * 0.6) .* tanh(tt*3);
  u  = 0*x;
  v  = pi * sin(2*pi*tt);
  w  = 0*x;
end

sfigure(1); plot3( x,y,z,'r.' ); hold on;
plot3(x(end),y(end),z(end),'mx','LineWidth',4); axis equal; hold off;

dlmwrite(saveToFile, [x y z 0*x 0*y 0*z], 'delimiter',',','precision','%4.6f');
setenv('saveToFileK',savePath);
!echo "successful save to csv file?" && ls -ltrh ${saveToFileK} | grep .xyzuvw

% Useful snippet:
%    cp -v  /u4`pwd`/*.xyzuvw ./
