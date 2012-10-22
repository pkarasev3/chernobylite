% A script to generate an input file for the simulator_client program.
% writes a simple CSV array to text, with appropriate extension .xyzuvw

objFileName = 'MiG-35.obj';
savePath    = '/u4/home/pkarasev3/source/ktrack/testdata/'
saveToFile  = [savePath, objFileName, '.xyzuvw']                           %#ok<NOPTS>


tt = linspace(0,10,5000); tt = tt(:);
x  = 10 * cos( 2*pi*tt ) .* tanh(tt);
y  = 10 * sin( 2*pi*tt * 0.7) .* tanh(tt);
z  = 20 + 6.5 * sin(2*pi*tt * 3.5) .* tanh(tt); 

sfigure(1); plot3( x,y,z,'r.' ); hold on;
plot3(x(end),y(end),z(end),'mx','LineWidth',4); axis equal; hold off;

dlmwrite(saveToFile, [x y z 0*x 0*y 0*z], 'delimiter',',','precision','%4.6f');
setenv('saveToFileK',savePath);
!echo "successful save to csv file?" && ls -ltrh ${saveToFileK} | grep .xyzuvw
