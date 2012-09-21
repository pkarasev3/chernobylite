function [xy0 g_prv g_f2f] = getCompensation( g_WC, g_prv, xy0, f )

  if isempty(g_prv) 
    g_prv = g_WC;
  end
  
  g_f2f = g_WC * (g_prv)^(-1);
  g_f2f(1:3,4) = 0; % zero out the translation ... hm...
  g_prv = g_WC; 
  
  u0 =  (xy0(1) - 640/2) * (1/f);
  v0 = -(xy0(2) - 480/2) * (1/f);

  uv = (g_f2f^(-1)) * [ u0(:)'; v0(:)'; ones(1,numel(v0)); ones(1,numel(v0)) ];
  
  x0 =  f * uv(1,:)./uv(3,:) + 640/2;
  y0 = -f * uv(2,:)./uv(3,:) + 480/2;
  
  xy0_old = xy0;
  xy0     = [x0,y0];
  img_coord_diff = norm( xy0(:) - xy0_old(:) );
  return;
  
end
