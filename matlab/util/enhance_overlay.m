function img = enhance_overlay( img, chan_idx )
  % input:  image  MxNx3 img (uint8) that is grey scale plus one channel overlay
  % chan_idx:  1,2,or 3 (rgb) which channel the overlay is
  % output: image with the overlay expanded / made more visible 
dbstop if error
  img = double(img); img = img / max(img(:));
  ir = img(:,:,1); ig = img(:,:,2); ib = img(:,:,3); 
  overlaychan = img(:,:,chan_idx);
  chan_idxB   = mod( (chan_idx+1),3 )+1;
  overlaychanB= img(:,:,chan_idxB );
  idx_overlay = find( 3*overlaychan(:) >= 2*(ir(:)+ig(:)+ib(:) ) ); % where is overlay
  
  [ii jj] = ind2sub(size(ir),idx_overlay);
  ii(ii<2)=2; jj(jj<2)=2; 
  ii(ii>size(img,1)-1)=size(img,1)-1; 
  jj(jj>size(img,2)-1)=size(img,2)-1;
  npts=numel(ii);
  for i = 1:npts 
    overlaychanB(ii(i)+1,jj(i)) = overlaychanB(ii(i)+1,jj(i))/2 + 0.5;
    overlaychanB(ii(i)-1,jj(i)) = overlaychanB(ii(i)+1,jj(i))/2 + 0.5;
    overlaychanB(ii(i),jj(i)+1) = overlaychanB(ii(i),jj(i)+1)/2 + 0.5;
    overlaychanB(ii(i),jj(i)-1) = overlaychanB(ii(i),jj(i)-1)/2 + 0.5;
    overlaychan(ii(i)-1,jj(i)-1) = 2*overlaychan(ii(i)-1,jj(i)-1) ;
    overlaychan(ii(i)+1,jj(i)+1) = 2*overlaychan(ii(i)+1,jj(i)+1) ;
    overlaychanB(ii(i)-1,jj(i)-1) = 1-overlaychanB(ii(i)-1,jj(i)-1) ;
    overlaychanB(ii(i)+1,jj(i)+1) = 1-overlaychanB(ii(i)+1,jj(i)+1) ;
  end
  
  img(:,:,chan_idx) = overlaychan;
  img(:,:,chan_idxB) = overlaychanB;
  img(img>1) = 1;
  img = uint8( 255*img );
  fprintf('');

end
