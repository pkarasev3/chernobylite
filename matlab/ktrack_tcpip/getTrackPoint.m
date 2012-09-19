function xyF = getTrackPoint( img, xy0, flag)

  xyF = xy0;
  if strcmp(flag,'local_max_bright')
    xyF = local_max_bright(img,xy0);
  end
  
end


function xyF = local_max_bright( img, xy0 )
  [H W c] = size(img);
  y0 = round(xy0(2)); x0 = round(xy0(1));
  sz = 10;
  xrange = x0-sz:x0+sz; xrange(xrange<1) = []; xrange(xrange>W) = [];
  yrange = y0-sz:y0+sz; yrange(yrange<1) = []; yrange(yrange>H) = [];
  subimg = rgb2gray( img( yrange, xrange,:) );
  subimg = imfilter(double(subimg),ones(3,3)/9,'replicate');
  [ym xm] = find( subimg == max(subimg(:)) );
  
  assert( 1-sz-1+x0 == x0-sz ); assert( 2*sz+1 -sz -1 + x0 == x0+sz );
  xyF     = [mean(xm(:))-sz-1+x0; mean(ym(:))-sz-1+y0];
  
  return;
end

