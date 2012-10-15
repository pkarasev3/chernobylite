%
% run this server first, then run something like ktrack client:
%         ./simulator_client -m ../../testdata/F22mod.3ds
import java.net.ServerSocket
import java.io.*
javaaddpath([pwd]);
addpath('~/source/chernobylite/matlab/display_helpers/');
addpath('~/source/chernobylite/matlab/util/');

b_compensateMotion = true();
b_computeHorizon   = true();
opts = struct('output_port',5001,'number_of_retries',1000,...
                    'compensation', b_compensateMotion,...
                    'horizon', b_computeHorizon);
disp(opts);

number_of_retries = opts.number_of_retries; % set to -1 for infinite
output_port       = opts.output_port;
retry             = 0;

if ~exist('server_socket','var')
  server_socket  = [];
  io_socket      = [];
  sh             = [];
else
  server_socket.close(); %#ok<SUSENS>
  io_socket.close(); %#ok<SUSENS>
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
    xy0      = [320,240]; % initial track-point
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
          fprintf(1,'bytes = %d, total= %d .. ', bytes_available, RecvBytes);
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
      
      fprintf( 'expected bytes = %d, got %d\n', expected_bytes, numel(data_raw(:)) );
      % Should be done receiving from client now, read it.
      [img, g_WC, f] = unpack_ktrack_data( data_raw, headerLen);
      fprintf('unpacked OK!\n');
      
      xy0prev=xy0;
      if opts.compensation
        [xy0 g_prv g_f2f] = getCompensation( g_WC, g_prv, xy0, f );
        fprintf('compensated OK!\n');
      end
      
      if opts.horizon
        [horizonU]    = getMetaHorizon( g_WC, img, f );
        horizon_show = horizonU .* (255+rgb2gray( double(img) ));
        sfigure(2); 
        imagesc(horizon_show); title('horizon metadata');
        fprintf('got horizon OK!\n');
      end
      
      
      % Run the tracker
      xyF = getTrackPoint( img, xy0, 'local_max_bright' );
      fprintf('trackpoint OK!\n');
      xy0_comp   = xy0;
      xy0        = xyF;
      
      % Generate the return bytes
      message = typecast( uint16([xyF(1),xyF(2)]),'uint8');
      sh=sfigure(1); imshow(img); hold on;
      
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
    if ~isempty(sh)
      close(sh);
    end
    break;
    
  catch                          %#ok<CTCH> 
    s = lasterror;               %#ok<LERR>
    bInterestingError = isempty(strfind(s.message,'java.net.SocketTimeoutException: Accept timed out'));
    if bInterestingError
      disp(['last error was: ']); disp(s.stack); disp(s.message);
    end
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
