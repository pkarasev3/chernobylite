function server_socket = ktrack_server(message, output_port, number_of_retries)
    %
    % run this server first, then run something like ktrack client:
    %         ./simulator_client -m ../../testdata/F22mod.3ds
    
    import java.net.ServerSocket
    import java.io.*
    javaaddpath([pwd]);
    addpath('~/source/chernobylite/matlab/display_helpers/');
    addpath('~/source/chernobylite/matlab/util/');
    
    if (nargin < 3)
        number_of_retries = 1000; % set to -1 for infinite
    end
    if nargin < 2 
        output_port = 5001;
    end
    if nargin < 1
        message = 'xyzabc';
    end
    retry             = 0;
    
    server_socket  = [];
    io_socket      = [];

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
            
            gotFrame = 1;
            headerLen      = 184; % sizeof( meta_data ) in cpp
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
                  fprintf(1,'Available bytes from client = %d\n', bytes_available);
                  data_reader = DataReader(d_input_stream);
                  recv        = data_reader.readBuffer(bytes_available);
                  data_raw    = [data_raw, recv(:)'];
                  RecvBytes   = RecvBytes + bytes_available;
                end
                pause(0.05); tB = toc( tA );
              end
              if RecvBytes == expected_bytes
                gotFrame = true;
              else
                gotFrame = false;
                continue;
              end

              % Should be done writing from client now, so read it.
              fprintf('reply(1:8) as raw char array is:\n%s\n',num2str(data_raw(1:8)));
              meta_data  = typecast(data_raw(1:headerLen),'double'); 
              g_WC       = reshape(meta_data(7:7+15),[4,4])'; 
              disp('g_WC = '); disp(g_WC);
              img_raw    = typecast(data_raw(headerLen+1:end),'uint8');
              B=reshape( img_raw(1:3:end),[640,480])';
              G=reshape( img_raw(2:3:end),[640,480])';
              R=reshape( img_raw(3:3:end),[640,480])';
              img        = uint8(zeros(480,640,3)); img(:,:,1)=R; img(:,:,2)=G; img(:,:,3)=B;

              message = typecast( uint16([320,240]),'uint8');
              sfigure(1); imagesc(img); title(['Received image from client. ' ...
                                                   'Server sends the reply:  ' message] );




              output_stream   = io_socket.getOutputStream();
              d_output_stream = DataOutputStream(output_stream);

              % output the data over the DataOutputStream
              fprintf(1, 'Writing %d bytes\n', length(message))
              d_output_stream.write( uint8(message) );
              d_output_stream.flush;
              pause(0.05);
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
end
