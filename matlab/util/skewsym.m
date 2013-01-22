function out = skewsym( u )
   if( numel(u) == 3 )
   out = [0 -u(3) u(2);
           u(3) 0 -u(1);
           -u(2) u(1) 0 ];
   elseif(numel(u) == 16 )
     z        = real(logm(u));
     out      = [z(3,2); -z(3,1); z(2,1)];
   end
end
