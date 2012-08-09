function img_show = draw_contour_pair_on_image( img0, phi1, phi2, phi_show_thresh )

  img_show = repmat(img0,[1 1 3]);
  imgb = img_show(:,:,3);
  imgg = img_show(:,:,2);
  imgr = img_show(:,:,1);

  % zero out the non-active colors for phi1 (active red), phi2 (active green)
    imgr( abs( phi2 ) < phi_show_thresh ) = 0;
    imgb( abs( phi2 ) < phi_show_thresh ) = 0;
    imgg( abs( phi2 ) < phi_show_thresh ) = 0; % bizarre... why need this ? 

  % make phi1 "patterned"
  sz = 4;
  [M,N] = size(imgr);
  
  [xx yy] = meshgrid( linspace(-1,1,N), linspace(-1,1,M) );
  [yc,xc] = find( phi1 == max(phi1(:)), 1, 'first' );
  xx      = xx - xx(yc,xc);
  yy      = yy - yy(yc,xc);
  theta   = atan2( yy, xx );
  winfunc = sin( 12 * pi * theta );
  phi1(winfunc > 0.2) = sign(phi1( winfunc > 0.2 ) ) * 10.0;
  
  phi2( abs(phi1)<phi_show_thresh ) = sign( phi2(abs(phi1)<phi_show_thresh ) ).*10.0;
  
    imgg( abs( phi1 ) < phi_show_thresh ) = 0;
    imgb( abs( phi1 ) < phi_show_thresh ) = 0;
    
  imgr( abs( phi1 ) < phi_show_thresh) = (imgr( abs( phi1 ) < phi_show_thresh) .* ... 
    abs( phi1(abs(phi1) < phi_show_thresh ) )/phi_show_thresh  + ...
    1 * (phi_show_thresh - abs( phi1(abs(phi1) < phi_show_thresh ) ) )/phi_show_thresh );

  
  imgg( abs( phi2 ) < phi_show_thresh) = (imgg( abs( phi2 ) < phi_show_thresh) .* ... 
    abs( phi2(abs(phi2) < phi_show_thresh ) )/phi_show_thresh  + ...
    1 * (phi_show_thresh - abs( phi2(abs(phi2) < phi_show_thresh ) ) )/phi_show_thresh );

  
  img_show(:,:,1) = imgr; img_show(:,:,2) = imgg; img_show(:,:,3) = imgb;
  img_show(img_show>1)=1; img_show(img_show<0)=0;
  
end



  %J=checkerboard(sz,M/(sz*2),N/(sz*2)); J( J > 0.5 ) = 1;
  %phi1( J > 0.5 ) = sign(phi1(J>0.5)).*10.0;
  %phi1 = imfilter( phi1, fspecial('gaussian',[3 3], 0.1 ), 'replicate' );
  
