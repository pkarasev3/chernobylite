function [xyF,pandora] = getTrackPoint( img, xy0, U, flag)
  
  bRunTest = false;
  pandora  = []; 
  if nargin == 0
    bRunTest = true;
    [img,xy0,flag] = createTestData();
    U = zeros( size(img,1),size(img,2) );
  end
  
  xyF = xy0;
  if strcmp(flag,'local_max_bright')
    xyF = local_max_bright(img,xy0);
  elseif strcmp(flag,'mean_align_levelset')
    xyF = levelset_means( img, xy0, U );
  end
  
  if bRunTest
    runTest(img);
  end
  
end

function [xyF] = levelset_means( img, xy0, U_in )
  
  global TKR;   
  global KOpts;
  bUseUControl =   false(); %true(); %
  
  img = rgb2gray(double(img) * 1.0/255.0);
  img = 10.0 * (img - min(img(:)))/(max(img(:))-min(img(:))+1e-9);
  if ~isfield(TKR,'compensate')   % (a) creating tracker! 
    params = struct('bUseUControl',bUseUControl,'Img',img);
    tkr = getLevelsetTracker( params );
    TKR = tkr; fprintf('\nCreated TKR; Should only happen once! \n');
  else              % (b) already exists! 
    tkr = TKR;
  end
  
  if KOpts.compensation 
    tkr.compensate( );
  end
  
  % Force a copy to happen, hackish
  tkr.psi(:) = TKR.psi(:);
  
  itrs = KOpts.contour_iters;
  for m = 1:itrs
    tkr.update_phi(img,U_in,m/itrs);
  end
  
  xyF     = tkr.get_center();
  TKR.xyF = xyF;
  
  tkr.display(img); % Careful, this responsible stuff moving tkr => TKR 
  
  pause(0.001);
end

function xyF = local_max_bright( img, xy0 )
  [H W c] = size(img);
  y0 = round(xy0(2)); x0 = round(xy0(1));
  sz = 20;
  x0 = max( x0, 2 ); x0 = min(x0,W-1);
  y0 = max( y0, 2 ); y0 = min(y0,W-1);
  xrange = x0-sz:x0+sz; xrange(xrange<1) = []; xrange(xrange>W) = [];
  yrange = y0-sz:y0+sz; yrange(yrange<1) = []; yrange(yrange>H) = [];
  subimg = rgb2gray( img( yrange, xrange,:) );
  subimg = imfilter(double(subimg),ones(3,3)/9,'replicate');
  


  [ym xm] = find( subimg == max(subimg(:)) );
  xyF     = [mean(xm(:))-sz-1+x0; mean(ym(:))-sz-1+y0];
  assert( 1-sz-1+x0 == x0-sz ); assert( 2*sz+1 -sz -1 + x0 == x0+sz );
  
end

function [img,xy0,flag] = createTestData()
  global TKR;
  TKR = [];
  W   = 2*320;
  H   = 2*240;
  img = 5.0 * (checkerboard(10,2,2) > 0.5);
  img = imresize( img, W/size(img,2) );
  img = imfilter(img,fspecial('gaussian',[16 16],5),'replicate');
  img = (img - min( img(:) ))/(max(img(:))-min(img(:)));
  img = repmat( img(1:H,1:W), [1 1 3]);
  xy0 = [size(img,2)/2 ; size(img,1)/2];
  flag= 'mean_align_levelset';

end

function runTest(img)
  global TKR;
  tkr = TKR;
  
  %  run a loop, updating tkr struct in-place 
  
  img0 = max(img(:))*rgb2gray(img); kMax = 300;
  tic;
  sfigure(1); clf;
  for k = 1:kMax
    img = img0 .* (1 + 0.05*randn(tkr.img_size(1),tkr.img_size(2)) );
    img = circshift(img,[round(k/4),round(k/3)]);
    tkr.U = 0*tkr.phi; 
    tkr.update_phi(img);                           
    xyF = tkr.get_center( );
    tkr.display(img); title( [ num2str(k)  '  of  ' num2str(kMax) , sprintf(',  xyF= %3.3f,%3.3f',xyF(1),xyF(2) ) ] );
    hold on; plot( xyF(1), xyF(2), 'rs','LineWidth',3 ); plot( xyF(1), xyF(2), 'mx','LineWidth',1 ); hold off;
    pause(0.05); drawnow();
  end
  timeTotal=toc;
  fprintf('iterations = %03d, Hz = %4.4f \n', k, 1/(timeTotal/kMax) );
  TKR = tkr;
end

