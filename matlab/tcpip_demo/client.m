% CLIENT connect to a server and read a message
%
% Usage - message = client(host, port, number_of_retries)
function message = client(host, port, number_of_retries)

    import java.net.Socket
    import java.io.*
    javaaddpath([pwd]);

    % testreader = DataReader([]);
    
    if (nargin < 3)
        number_of_retries = 1000; % set to -1 for infinite
    end
    if nargin < 2
        port = 5001;
    end
    if nargin < 1
        host = 'localhost';
    end
    
    retry        = 0;
    io_socket    = [];
    message      = [];

    while true

        retry = retry + 1;
        if ((number_of_retries > 0) && (retry > number_of_retries))
            fprintf(1, 'Too many retries\n');
            break;
        end
        
        try
            fprintf(1, 'Retry %d connecting to %s:%d\n', ...
                    retry, host, port);

            % throws if unable to connect
            io_socket = Socket(host, port);
            fprintf(1, 'Connected to server %s port %d\n', host, port);
            pause(0.5);
            
            % get a buffered data input stream from the socket
            client2server_stream   = io_socket.getOutputStream();
            d_client2server_stream = DataOutputStream(client2server_stream);
           
            img   = zeros(32,32); img(2:2:end,2:2:end) = 169;
            bytes = typecast(img(:),'uint8');
            d_client2server_stream.write( bytes );
            d_client2server_stream.flush;
            
            
            % read data from the socket - wait a short time first
            pause(1);
            server2client_stream   = io_socket.getInputStream();
            d_server2client_stream = DataInputStream(server2client_stream);
                        
            tA = tic(); tB = toc( tA );
            bytes_available = server2client_stream.available;
            fprintf(1,'Available bytes from client = %d\n', bytes_available);
            while (bytes_available == 0) && (tB < 10.0)
              bytes_available = server2client_stream.available;
              fprintf(1,'Available bytes from client = %d\n', bytes_available);
              pause(0.1); tB = toc( tA );
            end            
            
            fprintf(1, 'Client is reading %d bytes from server\n', bytes_available);         
            data_reader = DataReader(d_server2client_stream);
            message     = data_reader.readBuffer(bytes_available);
            
            message = char(message');
            
            % cleanup
            io_socket.close;
            break;
            
        catch
            if ~isempty(io_socket)
                io_socket.close;
            end
            
            s = lasterror
            s.stack 
            % pause before retrying
            pause(0.01);
        end
    end
end
