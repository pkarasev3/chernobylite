function [U deltaU num_inputs] = GenerateUserInput( phi_star, phi, prev_err, U, ...
                                                               imgForU, Umax, num_inputs )

% Generate and accumulate user inputs
  %num_inputs = 3;
    
  k = 1;
  U_= U;
  epsilon   = 2;
  Heavi     = @(z)  1 * (z >= epsilon) + (abs(z) < epsilon).*(1+z/epsilon+1/pi * sin(pi*z/epsilon))/2.0;
  
  while k <= num_inputs   
    % User is the only place that reference phi_star exists !
    % Simulate their input after some time. 
  
    % Find where there is disagreement with truth.
    % No input where error is naturally decreasing.
    idx_u = find( (abs( (phi_star > 0).*(  0 > phi ) - ...
      (phi_star < 0).*(  0 < phi ) )) ...
       .*(abs((phi)-(phi_star))>=abs(prev_err))          );
  
%     
%     currErr=phi(idx_u)-phi_star(idx_u);
%     prevErr=prev_err(idx_u);
%     idxDrop = find(abs( currErr(:)) < abs(prevErr(:)));
%     idx_u( abs( currErr(:)) < abs(prevErr(:)) ) = []; 
    
    
    idx_u2 = idx_u; 
    idx_u2( abs(phi(idx_u)) >= 4 ) = [];
    if ~isempty(idx_u2) 
      idx_u=idx_u2;
    end
    
    if( (numel(idx_u) < k ) )
      px = 1; py = 1;
    else
      idx_u   = idx_u( randperm(numel(idx_u)) );
      [py px] = ind2sub( size( phi ),idx_u(k) );
    end
    U_ = updateU( U_, phi_star,phi,px,py,imgForU,Umax);
    diffU=norm( U(:)-U_(:) );
    
    if( k<= 1);      fprintf('diffU = ');    end;
    %fprintf(' %6.2g ,  ',diffU);
    k  = k+1;
    if( k==num_inputs);      fprintf('\n');    end;
    
  end
  
  % Update U
  U               = U_;
  deltaU          = U-U_;
end
