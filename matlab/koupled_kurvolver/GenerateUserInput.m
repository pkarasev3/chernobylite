function [U deltaU num_inputs] = GenerateUserInput( phi_star, phi, U, ...
                                                               imgForU, Umax )

% Generate and accumulate user inputs
  num_inputs = 3;
  
  k = 1;
  U_= U;
  
  while k <= num_inputs   
    % User is the only place that reference phi_star exists !
    % Simulate their input after some time. 
    
    idx_u = find( abs( (phi_star > 0).*(  0 > phi ) - ...
      (phi_star < 0).*(  0 < phi ) ) > 0 );
    
    if( (numel(idx_u) < k ) )
      px = 1; py = 1;
    else
      idx_u   = idx_u( randperm(numel(idx_u)) );
      [py px] = ind2sub( size( phi ),idx_u(k) );
    end
    U_ = updateU( U_, phi_star,phi,px,py,imgForU,Umax);
    diffU=norm( U(:)-U_(:) );
    
    if( k<= 1);      fprintf('diffU = ');    end;
    fprintf(' %6.2g ,  ',diffU);
    k  = k+1;
    if( k==num_inputs);      fprintf('\n');    end;
    
  end
  
  % Update U
  U               = U_;
  deltaU          = U-U_;
end
