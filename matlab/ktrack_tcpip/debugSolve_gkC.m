function debugSolve_gkC( )

dbstop if error

% W1C1 is harder to interpret ...
%data=load('gctrlResults_W1C1.mat'); TKR=data.TKR; KOpts=data.KOpts; r=data.results;

data=load('gctrlResults_W0C1.mat'); TKR=data.TKR; KOpts=data.KOpts; r=data.results;
fprintf('\n');
sfigure(1); clf; hold on;
r.nFrame_in = [0;r.nFrame_in];
wA = []; wB = []; wF = [];
for m = 1:numel(r.ang_diff)
  
  g_f2f  = r.g_f2f(:,:,m);
  g_ctrl = r.g_ctrl(:,:,m);
  tx     = r.estm_xy(m,1)-320; %tx = -tx;
  ty     = r.estm_xy(m,2)-240; %ty = -ty;
  Nt     = r.nFrame_in(m+1)-r.nFrame_in(m);
  
  w_f2f        = skewsym( g_f2f ); 
  Kt           = 3 / sqrt(TKR.img_size(1) * TKR.img_size(2) );
  w_ctrlA      = [  Kt*ty;  % 
                    Kt*tx;  
                         0 ]*pi/180;
  w_ctrlA      = Nt * w_ctrlA;                       
  w_ctrlB      = skewsym( g_ctrl );
  
  % assert( sum( 0 > w_ctrlA(1:2) .* w_ctrlB(1:2) ) == 0 )
  
  plot( w_ctrlA(1), w_ctrlA(2), 'rx' )
  plot( w_ctrlB(1), w_ctrlB(2), 'bo' )
  plot( w_f2f(1), w_f2f(2), 'ks' )
  
  wA = [wA, w_ctrlA(:)]; wB = [wB, w_ctrlB(:)]; wF = [wF, w_f2f(:)];
  
  F = w_f2f(:)'
  A = w_ctrlA(:)'
  B = w_ctrlB(:)'
  
  fprintf(''); pause(0.001);
  
end
hold off;
sfigure(1);
plot( wA(1,:),'r--x'); hold on; plot( wF(1,:), 'k-s'); plot( wB(1,:), 'b--o'); 
sfigure(2); clf;
plot( wA(2,:),'r--x'); hold on; plot( wF(2,:), 'k-s'); plot( wB(2,:), 'b--o'); 

disp( logm(TKR.g_k_d)+logm(TKR.g_ctrl) ) ; 
disp( logm( TKR.g_f2f ) )
m            = tkr.img_size(1);
n            = tkr.img_size(2);
[xx yy]      = meshgrid(linspace(1,n,n),linspace(1,m,m));

z0           = -100.0;
xx           = -(xx - (n-1)/2) / tkr.f * z0;
yy           =  (yy - (m-1)/2) / tkr.f * z0;

tx           =  TKR.xyF(1) - (n)/2;
ty           =  TKR.xyF(2) - (m)/2;

g_f2f        = tkr.g_f2f;
z_f2f        = real(logm(g_f2f));
w_f2f        = [z_f2f(3,2); -z_f2f(3,1); z_f2f(2,1)]';
Kt           = 3 / sqrt(m*n);
w_ctrl       = [  Kt*ty;  % 
                  Kt*tx;  
                         0 ]'*pi/180;

% Bound the control appropriately
yaw_and_pitch = w_ctrl(1:2)*180/pi;
yaw_and_pitch(yaw_and_pitch<-0.25)=-0.25;
yaw_and_pitch(yaw_and_pitch> 0.25)= 0.25;
w_ctrl(1:2) = yaw_and_pitch*pi/180;

tauDelay     = TKR.curr_Nframe - TKR.prev_Nframe;
g_ctrl       = expm([ tauDelay * [ skewsym(w_ctrl), [0;0;0] ]; [0 0 0 0] ])
w_f2f_hat    = w_f2f- tauDelay * w_ctrl;
fprintf( 'Ndelay=%02d, wx=%4.4f, wy=%4.4f, wz=%4.4f \n',tauDelay,...
  w_f2f_hat(1),...
  w_f2f_hat(2),...
  w_f2f_hat(3) );

g_k_d = expm( real( logm( g_f2f ) - logm( g_ctrl ) ) ); fprintf('');

% aha, issue: g_k_d had -.0078, g_ctrl had 0.0044; more effort needs to be in g_ctrl! 
tmp1=( logm( g_k_d ) +logm(g_ctrl) ) 
tmp2=( logm( g_f2f ) )
tmp3= logm( g_ctrl )

% Better solution:  
%                     minimize   -norm( w_ctrl ) + norm(w_k_d), 
%                    subject to  w_f2f - w_ctrl = w_k_d  ,  |w_ctrl| < w_max

fprintf('');

  wMax = [0.5;0.5];
  w_f2fin=w_f2f(1:2);
  [w1,w2] = meshgrid(-1:0.05:1,-1:0.05:1);
  wE      = 0*w1;
  for m = 1:size(w1,1)
    for n = 1:size(w1,2)
      x = [w1(m,n); w2(m,n)];
      wE(m,n)= wcost(x);
    end
  end
    
  fprintf('');  

  function E = wcost( x )
    
    w_ctrl = x(:);
    w_k_d  = w_f2fin(:) - tauDelay*w_ctrl;
    E      = -sqrt( sum( w_ctrl.^2 ) ) + sqrt( sum(w_k_d.^2) );
    
  end

  function [cin,ceq] = wcon(x)
    w_ctrl = x(:);
    a      = abs(w_ctrl)-wMax;
    b      = -[ty;tx].*w_ctrl(:); 
    
    cin    = [a; b];
    ceq    = [];
  end

end
