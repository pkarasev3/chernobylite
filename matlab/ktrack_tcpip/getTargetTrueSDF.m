function psi = getTargetTrueSDF( Zbuffer, target_xy )
% given suggested x,y image point and depth values, 
%  extract via thresholding the target and return a "true" signed dist func

  global TKR;
  
  roi_sz = 36;
  xrange = target_xy(1)-(4/3)*roi_sz:target_xy(1)+(4/3)*roi_sz;
  yrange = target_xy(2)-roi_sz:target_xy(2)+roi_sz;
  xrange = round(xrange); yrange = round(yrange);
  xrange(xrange<1) = [];
  yrange(yrange<1) = [];
  xrange(xrange>size(Zbuffer,2)) = [];
  yrange(yrange>size(Zbuffer,1)) = [];
  assert( ~isempty(xrange) && ~isempty(yrange) );
  
  subZbuff = Zbuffer( yrange, xrange );
  minZ     = min( subZbuff(:) );
  
  zVarThresh = 10;
  psi        = -100  +     0*Zbuffer;
  subPsi     = -100  + 200.0* ( abs( subZbuff - minZ ) < zVarThresh );
  
  redistIters = 20; 
  dX = 1/sqrt(2); 
  subPsi = reinitializeLevelSetFunction( subPsi, 2, dX,redistIters, 2, 1, true() );
  
  psi(yrange,xrange) = subPsi;
  
  bDebugPsi = false;
  if bDebugPsi
    ShowPsi_Debug( psi );
  end
  
  TKR.psi = psi;
  TKR.xroi=xrange;
  TKR.yroi=yrange;
  
end

function ShowPsi_Debug( psi ) 
  sfigure(2); imagesc(psi .* (psi > -1) ); title('psi');
  drawnow;

end
