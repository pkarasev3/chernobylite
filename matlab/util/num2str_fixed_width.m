function str_num = num2str_fixed_width( k, num_chars )

if( nargin < 2 )
  num_chars = 4;
end

str_num = num2str(k);
if( norm( k - round(k) ) < 1e-9 )
  while( numel(str_num) < num_chars )
    str_num = ['0' str_num]; %#ok<AGROW>
  end
else % if it is floating point ...
  while( numel(str_num) < num_chars )
    str_num = [str_num '0']; %#ok<AGROW>
  end
  while( numel(str_num) > num_chars )
    str_num = [str_num(1:end-1)]; %#ok<AGROW>
  end
end

end
