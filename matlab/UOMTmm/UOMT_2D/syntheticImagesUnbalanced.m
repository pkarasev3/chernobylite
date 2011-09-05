function [mu0 mu1] = syntheticImagesUnbalanced(type, n)


    if mod(n,2)
        n=n+1;
    end
    %filter size
    fsize = round((n+1)/10);
    if ~mod(fsize,2)
        fsize = fsize+1;
    end
    if fsize < 3
        fsize = 3;
    end
    mu0 = ones(n+1);
    mu1 = ones(n+1);
    [X Y] = meshgrid(-n/2:n/2, -n/2:n/2);
    if strcmp(type, 'diffuseCirc')
    
        %DIFFUSE
        circ0 = X.^2+Y.^2 <= (n/8)^2;
        circ1 = X.^2+Y.^2 <= (n/4)^2;
        mu0(circ0) = 6;
        mu1(circ1) = 3;
        mu0 = conv2(mu0, 1/(fsize)^2*ones(fsize));
        mu1 = conv2(mu1, 1/(fsize)^2*ones(fsize));
    elseif strcmp(type, 'shiftCirc')
    
        %SHIFT
        circ0 = X.^2+Y.^2 <= (n/8)^2;
        circ1 = (X-n/20).^2+(Y+n/20).^2 <= (n/8)^2;
        mu0(circ0) = 6;
        mu1(circ1) = 6;
        mu0 = conv2(mu0, 1/(fsize)^2*ones(fsize));
        mu1 = conv2(mu1, 1/(fsize)^2*ones(fsize));
    end
    %crop to remove filter padding
    marg = (fsize-1);
    mu0 = mu0(1+marg:end-marg, 1+marg:end-marg);
    mu1 = mu1(1+marg:end-marg, 1+marg:end-marg);