cvx_path = genpath('~/source/cvx/');
addpath(cvx_path);
addpath('../../display_helpers/');
addpath('../../util/');

tt = linspace(0,20,500);
dt = tt(2)-tt(1);

xy1 = zeros(2,length(tt));
uv1 = zeros(2,length(tt));
xy2 = zeros(2,length(tt));
uv2 = zeros(2,length(tt));

x01  = [-10;0];
x02  = [10;0];
xF1  = [10;10];
xF2  = [-10;10];
mass1= 1/2;
mass2= 1/2;
xy1(:,1) = x01;
xy2(:,1) = x02;

dist_12  = zeros(1,length(tt));
sparse_ak= zeros(2,length(tt));
damp2    = 0.5;
for k = 2:length(tt)
  
  tk  = tt(k);
  xF1 = xF2;
  xF2 = 10*[cos(2*pi*tk/tt(end));sin(2*pi*tk/tt(end))];
  dist_12  = sqrt( (xy1(1,:) - xy2(1,:)).^2 + (xy1(2,:) - xy2(2,:)).^2 );
  
  ak1  = ((xF1 - xy1(:,k-1))*5 - (0.5*uv1(:,k-1)).^3 + 1e-1*randn(2,1))/mass1;
  vk1  = uv1(:,k-1) + dt*ak1;
  xy1(:,k) = xy1(:,k-1) + vk1*dt + 0.5*ak1*dt^2;
  uv1(:,k) = vk1;
  
  turnk = [0;0]; M = 2;
  if( (dist_12(:,k-1)<2) && (k>4) )
    if( dist_12(:,k-M) > 2 )
      bSmart = false();
      if( bSmart )
        turnk = -uv2(:,k-1); 
        damp2 = damp2 * 1.25;
      else
        turnk = -2*( [uv1(1,k-1); uv1(2,k-1) ] - [uv1(1,k-2); uv1(2,k-2) ] );
      end
      turnk = (1/(dt))*norm(uv1(:,k-1)) * turnk / (norm(turnk)+1e-6);
    end
  end
  sparse_ak(:,k) = turnk;
  ak2  = (xF2 - xy2(:,k-1) - (damp2*uv2(:,k-1)).^3 + 1e-1*randn(2,1))/mass2;
  ak2  = ak2 + turnk;
  vk2  = uv2(:,k-1) + dt*ak2;
  xy2(:,k) = xy2(:,k-1) + vk2*dt + 0.5*ak2*dt^2;
  uv2(:,k) = vk2;
  
  
  
  sfigure(1); subplot(1,2,1);
  plot(xy1(1,1:k),xy1(2,1:k),'r-');  hold on;
  plot(xy1(1,k),xy1(2,k),'rx'); plot(xy2(1,k),xy2(2,k),'bx');
  plot(xy2(1,1:k),xy2(2,1:k),'b-'); 
  hold off; subplot(1,2,2); plot(tt(1:k),sparse_ak(1,1:k),'m*'); hold on;
  plot(tt(1:k),sparse_ak(2,1:k),'co'); hold off;
  drawnow;
  fprintf('');
  
end

save reactive_data  xy1 xy2 uv1 uv2 sparse_ak tt dt
!ls -ltrh reactive_data*


