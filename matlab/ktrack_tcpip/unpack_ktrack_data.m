function [img, g_WC, f] = unpack_ktrack_data( data_raw, headerLen)
  dbstop if error
    meta_data  = typecast(data_raw(1:headerLen),'double');
    g_WC       = reshape(meta_data(7:7+15),[4,4])';
    
    bHasNoise  = (rand(1,1) > 0.95);
    if( bHasNoise )
      noise      = [ expm( skewsym(randn(3,1)*0.05) ) , 0.05*randn(3,1); [0 0 0 1] ];
      disp('before noise g_WC = '); disp(g_WC);
      g_WC       = g_WC * noise;
    end
    
    f          = meta_data(23); assert( (1e2 < f) && (f < 1e4) ); % ensure sane f
    disp('g_WC = '); disp(g_WC);
    img_raw    = typecast(data_raw(headerLen+1:end),'uint8');
    B=reshape( img_raw(1:3:end),[640,480])';
    G=reshape( img_raw(2:3:end),[640,480])';
    R=reshape( img_raw(3:3:end),[640,480])';
    img        = uint8(zeros(480,640,3)); 
    img(:,:,1)=R; 
    img(:,:,2)=G; 
    img(:,:,3)=B;
end
