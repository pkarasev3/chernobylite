% SERVER Write a message over the specified port
% 
% Usage - server(message, output_port, number_of_retries)
%
% Example: 
% server('abcdef0123456789alfghjldfgjls;dgj',5001)
%   should write 33 bytes to the client
function server(message, output_port, number_of_retries)

    import java.net.ServerSocket
    import java.io.*
    javaaddpath([pwd]);
    
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
                      
            % wait for 1 second for client to connect server socket
            server_socket = ServerSocket(output_port);
            server_socket.setSoTimeout(5000);

            io_socket = server_socket.accept(); 
            fprintf(1, 'Client connected\n');
           
            input_stream    = io_socket.getInputStream();
            d_input_stream  = DataInputStream(input_stream);
            
            tA = tic(); tB = toc( tA );
            bytes_available = input_stream.available;
            fprintf(1,'Available bytes from client = %d\n', bytes_available);
            while (bytes_available == 0) && (tB < 10.0)
              bytes_available = input_stream.available;
              fprintf(1,'Available bytes from client = %d\n', bytes_available);
              pause(0.1); tB = toc( tA );
            end
            
            % Should be done writing from client now, so read it.
            data_reader = DataReader(d_input_stream);
            reply       = data_reader.readBuffer(bytes_available);
            whos reply
            fprintf('reply(1:8) as raw char array is:\n%s\n',num2str(reply(1:8)));
            data_raw = typecast(reply(:),'double');
            img      = reshape( data_raw, [32 32] );
            figure(1); imagesc(img); title(['Received image from client. ' ...
                                                 'Server sends the reply:  ' message] );
            
            output_stream   = io_socket.getOutputStream();
            d_output_stream = DataOutputStream(output_stream);

            % output the data over the DataOutputStream
            fprintf(1, 'Writing %d bytes\n', length(message))
            d_output_stream.write( uint8(message) );
            d_output_stream.flush;
            
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
