function run_gradTowardOne_test()
  set(0,'defaultaxesfontsize',16);  
  set(0,'defaulttextfontsize',18);
  set(0,'defaulttextfontname','Arial');
  set(0,'defaultaxesfontweight','bold');
  set(0,'defaultlinelinewidth',2);
  set(0,'defaultlinemarkersize',4);
  
  dt_init         = 0.7; 
  
  test_gen_phi( );

  
 
end

 

function test_gen_phi(phi_star,phi,img)

if( nargin < 3 )
  disp('no args received, attempting to load data file...');
  data=load('data_run_U_accumulate_test.mat');
  phi_star = data.phi_star;
  phi      = data.phi;
  img      = data.img;
end
epsilon   = sqrt(2); %1.1;%0.8;
Heavi     = @(z)  1 * (z >= epsilon) + (abs(z) < epsilon).*(1+z/epsilon+1/pi * sin(pi*z/epsilon))/2.0;
delta     = @(z)  1 * (z == 0) + (abs(z) < epsilon).*(1 + cos(pi*z/epsilon))/(epsilon*2.0);

[X Y] = meshgrid(1:size(phi,2),1:size(phi,1));
sfigure(1); clf; sfigure(2); clf; 
dt0 = 1e-1;
phi0 = phi;

  k = 1; kmax = 1000;
  t = 0;
  Fval_all = zeros(1,kmax);
  Fval     = eval_gradOne_dist(phi);
  Fval_all(1) = Fval; 
  while( k < kmax )
    
      k    = k+1;
      rhs  = gradTowardsOne(phi);
      dt   = dt0 / (max(abs(rhs(:))) + 1e-15);
      t    = t + dt;
      Fval0=Fval;
      phi_prv = phi;
      phi  = phi + dt * rhs;
      Fval = eval_gradOne_dist(phi);
      if( Fval > Fval0 )
        phi = phi_prv;
        t   = t - dt;
        dt0 = dt0 * 0.8; fprintf('dt0=%f\n',dt0);
      else
        dt0 = dt0 * (1 + 0.05*(dt0<(1e-1-0.05)));
      end
      sfigure(2); imagesc(phi); title(['iter = ' num2str(k) ' of ' num2str(kmax) ... 
                                           ', t= ' num2str(t) ', Fval= ' num2str(Fval) ]); 
                                               
      Fval_all(k) = Fval; sfigure(1); plot(Fval_all(1:k),'r-.');  title('\int_{\Omega} (|\nabla \phi|^2-1|)^2  ');
      drawnow;
          
    fprintf('');
    
  end
  
  disp('done with test of U accumulate... saving');
  save data_run_gradTowardOne_test   phi phi_star img 
  !ls -ltrh ./data_run_gradTowardOne* 
end
