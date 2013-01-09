% test function for show_speedometer()

current_jcp   = javaclasspath();
needToAddJars = false;
if isempty(current_jcp)
  needToAddJars = true;
elseif isempty( strfind( current_jcp{1}, 'jfreechart-1.0.14' ) )
  needToAddJars = true;
end
if needToAddJars
  javaaddpath([pwd '/jfreechart-1.0.14/lib/jcommon-1.0.17.jar'])
  javaaddpath([pwd '/jfreechart-1.0.14/lib/jfreechart-1.0.14.jar'])
end
  

vin  = 50;
vref = 55;

plot_handleA= make_speedometer(0,50); % make a default figure with no input figure handle 
pause(1);

fh          = sfigure(2);
plot_handle = make_speedometer(vin,vref,fh);
pause(1); 

for vref = 0:0.25:160
  
  % Set arrow 1
  set_vref_arrow( plot_handle, vref );
  
  
  % Set arrow 2
  vin    = vref + 20*sin(pi/4 + 2*pi*vref/160);
  set_vin_arrow( plot_handle, vin );
  
  pause(0.01);
  

end

