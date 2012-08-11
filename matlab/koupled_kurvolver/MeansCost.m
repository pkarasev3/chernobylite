function E = MeansCost(Img,phi,lambda, Heavi)
   % Compute mean-separation cost
   %
   
   if ndims(phi)==3
     [mu_i, mu_o] = compute_means3(Img,phi);
     term1 = trapz(trapz(trapz(  Heavi(phi).*(Img-mu_i).^2 ) ) );
     term2 = trapz(trapz(trapz(  (1-Heavi(phi)).*(Img-mu_o).^2 ) ) );
     [gx gy gz] = gradient( Heavi(phi) );
     term3 = lambda*trapz(trapz(trapz(  sqrt(gx.^2+gy.^2+gz.^2) ) ) );
   else
     [mu_i, mu_o] = compute_means2(Img,phi);
     term1 = (trapz(trapz(  Heavi(phi).*(Img-mu_i).^2 ) ) );
     term2 = (trapz(trapz(  (1-Heavi(phi)).*(Img-mu_o).^2 ) ) );
     [gx gy] = gradient( Heavi(phi) );
     term3 = lambda*(trapz(trapz(  sqrt(gx.^2+gy.^2) ) ) );
   end
   E = term1+term2+term3;

  function  [mu_i mu_o] = compute_means3( Img,phi )
     mu_i = trapz(trapz(trapz(Heavi( phi ) .* Img))) / trapz(trapz(trapz(Heavi( phi ) ) ));
     mu_o = trapz(trapz(trapz( (1-Heavi( phi )) .* Img))) / trapz(trapz(trapz( (1-Heavi( phi )) ) ));
  end

  function  [mu_i mu_o] = compute_means2( Img,phi )
       mu_i = (trapz(trapz(Heavi( phi ) .* Img))) / (trapz(trapz(Heavi( phi ) ) ));
       mu_o = (trapz(trapz( (1-Heavi( phi )) .* Img))) / (trapz(trapz( (1-Heavi( phi )) ) ));
  end

end


