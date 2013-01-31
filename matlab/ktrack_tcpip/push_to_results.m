function results = push_to_results( results )
  global KOpts;
  global TKR;
  dbstop if error
  
  if isempty(results)
    fprintf('first round, returning empty resluts!\n');
    results = struct('true_xy',[TKR.true_xy(:)'],...  % received target centroid hint
                     'estm_xy',[TKR.true_xy(:)'],...  % tracker result
                     'D_ls_err',[0],... % D, error between levelsets phi,psi
                     'nFrame_in',[TKR.curr_Nframe],... % source iters, best case equal to tracker iter
                     'ang_diff',[0],...
                     'g_f2f',zeros(4,4,KOpts.max_input_frames),...
                     'g_ctrl',zeros(4,4,KOpts.max_input_frames),...
                     'w0wX',zeros(3,2,KOpts.max_input_frames));
    results.g_f2f(:,:,1)=eye(4,4);% [results.g_f2f, eye(4,4)];
    results.g_ctrl(:,:,1)=eye(4,4);%    = [results.g_ctrl, eye(4,4)];
    return;
  end

  true_xy     = TKR.true_xy;
  xyF         = TKR.xyF;
  Nframe      = TKR.curr_Nframe;
  
  % restrict these to an roi...
  phi     = TKR.phi(  TKR.yroi,TKR.xroi);
  psi     = TKR.phi0( TKR.yroi,TKR.xroi);
  
  img_show= TKR.img_show;
  iter    = size( results.true_xy,1 );

  D_err   = 0.5*trapz(trapz( (TKR.Heavi(phi) - TKR.Heavi(psi)).^2 ) );
  
  z_f2f    = real(logm(TKR.g_f2f));
  ang_diff = norm( [z_f2f(2,1), z_f2f(3,1), z_f2f(3,2)]);
  
  results.ang_diff  = [results.ang_diff; ang_diff];
  results.true_xy   = [results.true_xy; true_xy(:)'];
  results.estm_xy   = [results.estm_xy; xyF(:)'];
  results.D_ls_err  = [results.D_ls_err; D_err];
  results.nFrame_in = [results.nFrame_in; Nframe];
  results.g_f2f(:,:,iter+1)     = TKR.g_f2f;%[results.g_f2f, TKR.g_f2f];
  results.g_ctrl(:,:,iter+1)    = TKR.g_ctrl;%[results.g_ctrl, TKR.g_ctrl];
  results.w0wX(:,:,iter+1)      = TKR.compare_w0_wX;
  
  if 0 < KOpts.showImages
    if mod(iter,KOpts.showImages)==0
      sfigure(1); imshow( TKR.img_show ); 
    end
    if mod(iter,KOpts.saveImages)==0
      filename=sprintf('ktracker_%04d.png',iter); 
      imwrite( uint8(255*TKR.img_show),filename ); fprintf([' wrote to: ', filename]);
    end
    %sfigure(2); imagesc( [phi, psi] ); 
    % ... save 'em, but will it affect speed ?! 
  end
  
end
      
