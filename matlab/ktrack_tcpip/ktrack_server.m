%
% run this server first, then run something like ktrack client:
%         ./simulator_client -m ../../testdata/F22mod.3ds


if ~exist('JavaAndPathsAreSetUp','var')
  import java.net.ServerSocket
  import java.io.*
  javaaddpath([pwd]);
  JavaAndPathsAreSetUp = true();
else
  disp(['java and paths are already set, continuing...']);
end
addpath('~/source/chernobylite/matlab/display_helpers/');
addpath('~/source/chernobylite/matlab/util/');
addpath('~/source/chernobylite/matlab/display_helpers/');
addpath('~/source/chernobylite/matlab/LSMlibPK/');

opts = getKtrackOpts(); % Internal options set here 
global KOpts;
KOpts = opts;

results           = [];

number_of_retries = opts.number_of_retries; % set to -1 for infinite
output_port       = opts.output_port;
retry             = 0;

!killall -9 simulator_client 
if ~exist('server_socket','var')
  server_socket  = [];
  io_socket      = [];
  sh             = [];
else
  server_socket.close(); %#ok<SUSENS>
  io_socket.close(); %#ok<SUSENS>
end
global TKR;  TKR = [];


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
    expected_bytes = headerLen + 640 * 480 * 3 + ... % header, rgb (8uc3),
                                 640 * 480 * 4;      %   and zbuffer (float32)
    RecvTimeout    = 5.0; % seconds in timeout per frame
    
    while gotFrame
      tA = tic(); tB = toc( tA );
      RecvBytes      = 0; % how many we've gotten so far
      data_raw       = [];
      bytes_available = input_stream.available;
      fprintf(1,'Available bytes from client = %d\n', bytes_available);
      while (RecvBytes < expected_bytes) && (tB < RecvTimeout) %(bytes_available == 0)
        bytes_available = input_stream.available;
        if bytes_available > 0
          print_Bytes = false;
          if print_Bytes; fprintf(1,'bytes = %d, total= %d .. ', bytes_available, RecvBytes); 
          end;
          
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
      
      % Unpack the stream of bytes, "custom protocol" 
      [img, g_WC, f, ...
           true_xy,true_Nframe, Zbuffer] = unpack_ktrack_data( data_raw, headerLen);
      fprintf('unpacked OK!\n');
      
      if true_Nframe > opts.max_input_frames
        fprintf('finished requested # of input frames!\n');
        fprintf('attempting to save results as %s\n', ... 
                                          KOpts.result_filename);
        save(opts.result_filename,'results'); gotFrame = false();
      end
      
      xy0prev=xy0;
      if true % Doh, need to record to TKR here.. % if opts.compensation
        [xy0 g_prv g_f2f] = getCompensation( g_WC, g_prv, xy0, f );
        fprintf('compensated OK!\n');
      end
      
      if opts.horizon
        [horizonU]    = getMetaHorizon( g_WC, img, f );
        fprintf('got horizon OK!\n');
      else
        horizonU = zeros( size(img,1), size(img,2) );
      end
      
      if opts.getPsiTru
        psi = getTargetTrueSDF( Zbuffer, true_xy );
        fprintf('got "true" sdf OK!\n');
      end
         
      % Run tracker. Inputs should be: {img, g_WC, U, f, tracker_type}
      localMaxBright = 'local_max_bright';
      levelsetMeans  = 'mean_align_levelset';
      trackerType    = levelsetMeans;
      xyF = getTrackPoint( img, xy0, horizonU, trackerType);
      fprintf('trackpoint OK!\n');
      xy0        = xyF;
     
      % Store into results for later evaluation
      TKR.true_xy     = true_xy;
      TKR.true_Nframe = true_Nframe;
      TKR.Zbuffer     = Zbuffer;
      TKR.xyF         = xyF;
      results         = push_to_results( results );
      
      % Generate the return bytes
      if true_Nframe > opts.max_input_frames 
        xyF = [666;666]; % send kill code to simulator
      end
      message = typecast( uint16([xyF(1),xyF(2)]),'uint8');
     
      sfigure(1); 
      title([ sprintf('Received image %05d, x=%3.2f,y=%3.2f', ...
         frameIdx,xyF(1),xyF(2)),', Server img#: ' num2str(true_Nframe) ] );
      
      
      % Return data to client via stream
      output_stream   = io_socket.getOutputStream();
      d_output_stream = DataOutputStream(output_stream);
      
      % output the data over the DataOutputStream
      fprintf(1, 'Writing %d bytes\n', length(message));
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
      s_ = s;
      disp(['last error was: ']); 
      s.stack.file
        s.stack.line
        s.stack.name
      %s.message
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
% 
%  if ~strcmp( levelsetMeans, trackerType ) 
%         sh=sfigure(1); imshow(img); hold on;
% 
%         plot( [xy0_comp(1), xyF(1)], [xy0_comp(2), xyF(2)],...
%           'c-o', 'MarkerSize',10,'LineWidth',2);
%         plot( [xy0prev(1), xy0_comp(1)], [xy0prev(2), xy0_comp(2)],...
%           '-mo','MarkerSize',4, ...
%           'MarkerFaceColor','m','LineWidth',2);
% 
%         hold off; 
%       end
