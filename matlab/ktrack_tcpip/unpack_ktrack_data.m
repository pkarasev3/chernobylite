function [img, g_WC, f, true_xy, true_Nframe, Zbuffer ] = ...
                       unpack_ktrack_data( data_raw, headerLen)
  dbstop if error
    meta_data  = typecast(data_raw(1:headerLen),'double');
    g_WC       = reshape(meta_data(7:7+15),[4,4])';
    
    bHasNoise  = false() && (rand(1,1) > 0.95);
    if( bHasNoise )
      noise      = [ expm( skewsym(randn(3,1)*0.05) ) , 0.05*randn(3,1); [0 0 0 1] ];
      disp('before noise g_WC = '); disp(g_WC);
      g_WC       = g_WC * noise;
    end
    
    
    f          = meta_data(23); assert( (1e2 < f) && (f < 1e4) ); % ensure sane f
    disp('g_WC = '); disp(g_WC);
    
    true_xyz1_W = [meta_data(3:5)';1];
    true_xyz1_C = g_WC * true_xyz1_W;
    true_xy     = [-1;1] .* (f * true_xyz1_C(1:2)/true_xyz1_C(3)) + [320; 240];
    true_Nframe = meta_data(6); % number of frame in source stream, which runs faster than algorithm
    
    Irgb_end    = headerLen + 640*480*3;
    Zbuf_end    = headerLen + 640*480*3 + 640*480*4 ;
    
    img_raw    = typecast(data_raw(headerLen+1:Irgb_end),'uint8');
    B=reshape( img_raw(1:3:end),[640,480])';
    G=reshape( img_raw(2:3:end),[640,480])';
    R=reshape( img_raw(3:3:end),[640,480])';
    img        = uint8(zeros(480,640,3)); 
    img(:,:,1)=R; 
    img(:,:,2)=G; 
    img(:,:,3)=B;
    
    zbv_raw    = typecast(data_raw(Irgb_end+1:Zbuf_end),'single');
    Zbuffer    = reshape( zbv_raw(:), [640 480] )';
    
    % careful, must agree with source
      zNear = 10.0; 
      zFar = 1000.0; 
      
    % get the actual depth, zbuffer is stored nonlinearly
    z_n     = 2.0 * Zbuffer - 1.0;
    Zbuffer = 2.0 * zNear * zFar ./ (zFar + zNear - z_n * (zFar - zNear));
    
    bDebugZbuffer = false;
    if bDebugZbuffer
      ShowZbuff_Debug( Zbuffer );
    end
    
end

function ShowZbuff_Debug( zb ) 
  sfigure(2); imagesc(zb); title('zbuffer input.');
  zmin = min(zb(:)); 
  zmax = max(zb(:));
  fprintf(' min = %4.4f, max = %4.4f\n', zmin, zmax);
  drawnow;

end
