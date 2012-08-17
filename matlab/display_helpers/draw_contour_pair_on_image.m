function img_show = draw_contour_pair_on_image( img0, phi1, phi2, phi_show_thresh, U, Umax)
  
  bDrawScalarField_U = false();
  if exist('U','var') % if there's a U input, it is a scalar field and we need a Umax
    bDrawScalarField_U = true(); 
    assert(numel(U) == numel(phi1) ); 
    assert(true==exist('Umax','var'));
  end

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
  [ycm,xcm] = find( phi1 >= 0.0 );
  yc = round((max(ycm(:))+min(ycm(:)))/2.0); %[yc,xc] = find( phi1 == max(phi1(:)), 1, 'first' );
  xc = round((max(xcm(:))+min(xcm(:)))/2.0);
  xx      = xx - xx(yc,xc);
  yy      = yy - yy(yc,xc);
  theta   = atan2( yy, xx );
  winfunc = sin( 15 * pi * theta );
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

  if bDrawScalarField_U  % Try to color in a checker-like pattern for scalar field
    Ushow = imfilter(U,fspecial('gaussian',[5 5],2),'replicate');
    Ushow1=0*Ushow; Ushow1(1:2:end,1:2:end) = Ushow(1:2:end,1:2:end);
    Ushow2=0*Ushow; Ushow2(2:2:end,2:2:end) = Ushow(2:2:end,2:2:end);
    [idxUy1 idxUx1] = find( (abs(Ushow1)>Umax/4) > 0 );
    [idxUy2 idxUx2] = find( (abs(Ushow2)>Umax/4) > 0 );
    if ~( min([numel(idxUy1),numel(idxUy2)]) <10 )
      idxU1  = sub2ind(size(U),idxUy1,idxUx1);
      idxU2  = sub2ind(size(U),idxUy2,idxUx2);
      imgb(idxU1) = 0.5 + 2*abs(Ushow1(idxU1))/Umax;
      imgr(idxU1) = 0.5 * ( (Ushow1(idxU1)<0)); %imgr(idxU1);
      imgg(idxU1) = 0.5 * (1 - (Ushow1(idxU1)<0)); %imgg(idxU1);
      %imgb(idxU2) = abs(Ushow(idxU2))/Umax;
      %imgr(idxU2) = abs(Ushow(idxU2))/Umax;
      %imgg(idxU2) = abs(Ushow(idxU2))/Umax;
    end
  end
  
  img_show(:,:,1) = imgr; img_show(:,:,2) = imgg; img_show(:,:,3) = imgb;
  img_show(img_show>1)=1; img_show(img_show<0)=0;
  
end



  %J=checkerboard(sz,M/(sz*2),N/(sz*2)); J( J > 0.5 ) = 1;
  %phi1( J > 0.5 ) = sign(phi1(J>0.5)).*10.0;
  %phi1 = imfilter( phi1, fspecial('gaussian',[3 3], 0.1 ), 'replicate' );
  
