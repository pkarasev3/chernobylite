function test()
set(0,'defaultaxesfontsize',16);
set(0,'defaulttextfontsize',18);
set(0,'defaulttextfontname','Arial');
set(0,'defaultaxesfontweight','bold');
set(0,'defaultlinelinewidth',2);
set(0,'defaultlinemarkersize',12);
dbstop if error;

PauseLen = 0.01;
x0 = [15;0];
tspan = linspace( 0, 20, 1000); %[0,20];
dstar = 1;
opts  = odeset('OutputFcn',@display_callback,'MaxStep',0.1);

sh    = sfigure(1);
sh    = drawCarAndWall( x0(1), dstar,sh);
[T,X]=ode45(@myode, tspan, x0, opts);

sfigure(2); plot(T,X(:,1)); 

breakhere=1;


  function xdot = myode( t, x )
    xdot = [0;0];
    xdot(1) = x(2);
    xdot(2) = -1.0*x(2) + 1*(dstar-x(1)); 
  end


  function status=display_callback(t,x,flag)
    switch flag
      case {'init',''}
        d  = x(1,end);
        sh = drawCarAndWall( d, dstar, sh); 
        drawnow;
        pause( PauseLen );
        status= (d<0);
      case {'done'}
        status=0; % return 1 to terminate it
    end
    
  end


end
