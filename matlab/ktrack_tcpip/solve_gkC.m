function [g_ctrl, TKR ] = solve_gkC( TKR, KOpts )
    
    m            = TKR.img_size(1);
    n            = TKR.img_size(2);  
    tx           = TKR.xyF(1) - (n)/2; 
    ty           = TKR.xyF(2) - (m)/2;

    g_f2f        = TKR.g_f2f; 
    z_f2f        = real(logm(g_f2f));
    w_f2f        = [z_f2f(3,2); -z_f2f(3,1); z_f2f(2,1)]';
    Kt           = 1 / sqrt(m*n);
    w_ctrl       = [Kt*ty; Kt*tx; 0 ]'*pi/180;
   
    if ~KOpts.gkC_smart
      w_ctrl = 0*w_ctrl; 
    end
    
    % Bound the control appropriately
    yaw_and_pitch = w_ctrl(1:2)*180/pi;
    
    Wmax = 0.5;
    yaw_and_pitch(yaw_and_pitch<-Wmax)=-Wmax;
    yaw_and_pitch(yaw_and_pitch> Wmax)= Wmax;
    w_ctrl(1:2)  = yaw_and_pitch*pi/180;
    tauDelay= TKR.curr_Nframe - TKR.prev_Nframe;    
    w_ctrl_0     = w_ctrl * tauDelay;
    g_ctrl       = expm([ tauDelay * [ skewsym(w_ctrl), [0;0;0] ]; [0 0 0 0] ]);
    w_f2f_hat    = w_f2f- tauDelay * w_ctrl;

    bSolveWithFminCon = KOpts.gkC_fmincon;
    if bSolveWithFminCon
      
      wMax    =  tauDelay * [  Wmax;  Wmax ]*pi/180;
      wMin    =  tauDelay * [ -Wmax; -Wmax ]*pi/180;
      tmpFu   = [ty,tx];
      wS      = 3;
      for LoL = 1:2
        if tmpFu(LoL) > 0
          wMax(LoL) =  min( [Wmax, wS*w_ctrl_0(LoL)] );
          wMin(LoL) = 0;
        else
          refv = max( [-Wmax, wS*w_ctrl_0(LoL)] );
          % should be equivalent:
          wMin(LoL) = -min( [Wmax, wS*abs(w_ctrl_0(LoL)) ] ); 
          assert( abs(refv - wMin(LoL) ) < 1e-9 );
          wMax(LoL) = 0;
        end
      end
      %wMax
      %wMin
      assert( sum( wMax < wMin ) == 0 );
      w_f2fin =  w_f2f(1:2);
      xydirin =  [tx;ty] ./ sqrt(1e-15+tx.^2+ty.^2);

      opts    = optimset('Display','off','MaxFunEvals',10^3,.... %'final','iter-detailed'
                         'Algorithm','active-set','TolFun',1e-8,'TolX',1e-8);
      x0      = w_ctrl(1:2);
      cin0    = wcon(x0);
      [x,fval,exitflag,output,lambda]= ...
                    fmincon(@wcost,x0,[],[],[],[],wMin,wMax,@wcon,opts); %disp(output);                  
      cin = wcon(x);
      cin0_cinX = [cin0(:), cin(:)]; % compare constraint violation of initial guess and found solution
      w_ctrlX = [x(:)',0];
    else
      w_ctrlX = tauDelay*w_ctrl;
    end
    TKR.compare_w0_wX = [w_ctrl(:)*tauDelay, w_ctrlX(:)];
    w_f2f_hat = w_f2f - w_ctrlX;
    g_ctrl    = expm([ [ skewsym(w_ctrlX), [0;0;0] ]; [0 0 0 0] ]);
    fprintf( 'Ndelay=%02d, wx=%4.4f, wy=%4.4f, wz=%4.4f \n',tauDelay,...
                                                            w_f2f_hat(1),...
                                                            w_f2f_hat(2),...
                                                            w_f2f_hat(3) );

 function E = wcost( x )
      wC     = x(:);
      w_k_d  = w_f2fin(:) - wC; % tauDelay removed... already applied it at init
      E      = ( sum( w_k_d.^2 ) ) ;
    end

    function [cin,ceq] = wcon(x)
      wC = x(:);
      if ( min( [ sum(abs(xydirin)), sum(abs(wC))] ) < 1e-9 )
        cin  = [0];
      else
        xydir_w= [wC(2);wC(1)]/sqrt(sum(wC(1:2).^2));  % want <w1,w2> damn close to 1 (colinear)
        cin   = [-((xydir_w'*xydirin)-0.995) ];
      end
      ceq    = [];      %ceq(x) == 0
    end
  
end
