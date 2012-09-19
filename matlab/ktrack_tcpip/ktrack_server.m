%
% run this server first, then run something like ktrack client:
%         ./simulator_client -m ../../testdata/F22mod.3ds
import java.net.ServerSocket
import java.io.*
javaaddpath([pwd]);
addpath('~/source/chernobylite/matlab/display_helpers/');
addpath('~/source/chernobylite/matlab/util/');

b_compensateMotion = true();

opts = struct('output_port',5001,'number_of_retries',1000,...
                    'compensation', b_compensateMotion);
disp(opts);

number_of_retries = opts.number_of_retries; % set to -1 for infinite
output_port       = opts.output_port;
retry             = 0;

if ~exist('server_socket','var')
  server_socket  = [];
  io_socket      = [];
else
  server_socket.close();
  io_socket.close();
end

while true
  
  retry = retry + 1;
  
  try
    if ((number_of_retries > 0) && (retry > number_of_retries))
      fprintf(1, 'Too many retries\n');
      break;
    end
    
    fprintf(1, ['Try %d waiting for client to connect to this ' ...
      'host on port : %d\n'], retry, output_port);
    
    % wait for 5 seconds for client to connect server socket
    server_socket = ServerSocket(output_port);
    server_socket.setReceiveBufferSize( 640 * 480 * 3 + 1024 );
    server_socket.setSoTimeout(5000);
    io_socket = server_socket.accept();
    fprintf(1, 'Client connected\n');
    
    input_stream    = io_socket.getInputStream();
    d_input_stream  = DataInputStream(input_stream);
    
    frameIdx = 0;
    g_prv    = [];
    xy0      = [350,240]; % initial track-point
    gotFrame = 1;
    headerLen      = 192; % sizeof( meta_data ) in cpp
    expected_bytes = headerLen + 640 * 480 * 3; % most we can write from matlab seems to be: 250240
    
    while gotFrame
      tA = tic(); tB = toc( tA );
      RecvBytes      = 0; % how many we've gotten so far
      data_raw       = [];
      bytes_available = input_stream.available;
      fprintf(1,'Available bytes from client = %d\n', bytes_available);
      while (RecvBytes < expected_bytes) && (tB < 5.0) %(bytes_available == 0)
        bytes_available = input_stream.available;
        if bytes_available > 0
          fprintf(1,'bytes = %d .. ', bytes_available);
          data_reader = DataReader(d_input_stream);
          recv        = data_reader.readBuffer(bytes_available);
          data_raw    = [data_raw, recv(:)'];
          RecvBytes   = RecvBytes + bytes_available;
        end
        pause(0.001); tB = toc( tA );
      end
      fprintf('\n');
      if RecvBytes == expected_bytes
        gotFrame = true;
        frameIdx = frameIdx+1;
      else
        gotFrame = false;
        continue;
      end
      
      % Should be done receiving from client now, read it.
      meta_data  = typecast(data_raw(1:headerLen),'double');
      g_WC       = reshape(meta_data(7:7+15),[4,4])';
      f          = meta_data(23); assert( (1e2 < f) && (f < 1e4) ); % ensure sane f
      disp('g_WC = '); disp(g_WC);
      img_raw    = typecast(data_raw(headerLen+1:end),'uint8');
      B=reshape( img_raw(1:3:end),[640,480])';
      G=reshape( img_raw(2:3:end),[640,480])';
      R=reshape( img_raw(3:3:end),[640,480])';
      img        = uint8(zeros(480,640,3)); img(:,:,1)=R; img(:,:,2)=G; img(:,:,3)=B;
      
      xy0prev=xy0;
      if opts.compensation
        [xy0 g_prv g_f2f] = getCompensation( g_WC, g_prv, xy0, f );
      end
      
      % Run the tracker
      xyF = getTrackPoint( img, xy0, 'local_max_bright' );
      xy0_comp   = xy0;
      xy0        = xyF;
      
      % Generate the return bytes
      message = typecast( uint16([xyF(1),xyF(2)]),'uint8');
      sfigure(1); imshow(img); hold on;
      
      plot( [xy0_comp(1), xyF(1)], [xy0_comp(2), xyF(2)],...
        'c-o', 'MarkerSize',10,'LineWidth',2);
      plot( [xy0prev(1), xy0_comp(1)], [xy0prev(2), xy0_comp(2)],...
        '-mo','MarkerSize',4, ...
        'MarkerFaceColor','m','LineWidth',2);
      
      hold off; title([ sprintf('Received image %05d, x=%3.2f,y=%3.2f', ...
        frameIdx,xyF(1),xyF(2)),'Server reply:' message] );
      
      % Return data to client via stream
      output_stream   = io_socket.getOutputStream();
      d_output_stream = DataOutputStream(output_stream);
      
      % output the data over the DataOutputStream
      fprintf(1, 'Writing %d bytes\n', length(message))
      d_output_stream.write( uint8(message) );
      d_output_stream.flush;
      pause(0.01);
    end
    
    % clean up
    server_socket.close;
    io_socket.close;
    break;
    
  catch
    %s = lasterror
    if ~isempty(server_socket)
      server_socket.close
    end
    
    if ~isempty(io_socket)
      io_socket.close
    end
    
    % pause before retrying
    pause(0.01);
  end
end
%               hold on;
%img = silent_rectangle([xyF(2),xyF(1)],20,img,[255 50 100]);
%img = silent_rectangle([xyF(2),xyF(1)],18,img,[ 0 0 0]);
%               plot( xyF(1), xyF(2), 'c+', 'MarkerSize',12 ); hold off;

%plot( [xy0prev(1)], [xy0prev(2)],'r+','MarkerSize',8,'LineWidth',2);
%plot( xyF(1), xyF(2), 'gs', 'MarkerSize',14,'LineWidth',3);
%       plot( [xy0_comp(1)], [xy0_comp(2)],'s','MarkerEdgeColor','y',...
%                            'MarkerSize',8,'LineWidth',2,'MarkerFaceColor',[0 0 0]);
%
