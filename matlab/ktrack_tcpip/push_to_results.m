function results = push_to_results( results )
  global KOpts;
  global TKR;
  dbstop if error
  
  if isempty(results)
    fprintf('first round, returning empty resluts!\n');
    results = struct('true_xy',[],...  % received target centroid hint
                     'estm_xy',[],...  % tracker result
                     'D_ls_err',[],... % D, error between levelsets phi,psi
                     'nFrame_in',[]);  % source iters, best case equal to tracker iter
    return;
  end

  true_xy     = TKR.true_xy;
  xyF         = TKR.xyF;
  Nframe      = TKR.true_Nframe;
  
  % restrict these to an roi...
  phi     = TKR.phi(TKR.yroi,TKR.xroi);
  psi     = TKR.psi(TKR.yroi,TKR.xroi);
  
  img_show= TKR.img_show;
  iter    = size( results.true_xy,1 );

  D_err   = 0.5*trapz(trapz( (TKR.Heavi(phi) - TKR.Heavi(psi)).^2 ) );
  
  results.true_xy   = [results.true_xy; true_xy(:)'];
  results.estm_xy   = [results.estm_xy; xyF(:)'];
  results.D_ls_err  = [results.D_ls_err; D_err];
  results.nFrame_in = [results.nFrame_in; Nframe];
  
  if KOpts.saveImages
    sfigure(1); imshow( TKR.img_show ); 
    %sfigure(2); imagesc( [phi, psi] ); 
    % ... save 'em, but will it affect speed ?! 
  end
  
end
      
