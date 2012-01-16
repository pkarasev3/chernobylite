%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% UPWIND_HJ_ENO3_2D() computes the third-order upwind HJ ENO
% approximation to grad(phi)
%
% Usage: [phi_x, phi_y] = ...
%        UPWIND_HJ_ENO3_2D(phi, vel_x, vel_y, ghostcell_width, dX)
%
% Arguments:
% - phi:              function for which to compute upwind
%                       derivative
% - vel_x:            x-component of velocity to use in upwinding
% - vel_y:            y-component of velocity to use in upwinding
% - ghostcell_width:  number of ghostcells at boundary of
%                       computational domain
% - dX:               array containing the grid spacing
%                       in coordinate directions
%
% Return values:
% - phi_x:            x-component of third-order, upwind 
%                       HJ ENO derivative
% - phi_y:            y-component of third-order, upwind 
%                       HJ ENO derivative
%
% NOTES:
% - The vel_x and vel_y arrays are assumed to be the same size.
%
% - phi_x and phi_y have the same ghostcell width as phi.
%
% - All data arrays are assumed to be in the order generated by the
%   MATLAB meshgrid() function.  That is, data corresponding to the
%   point (x_i,y_j) is stored at index (j,i).  The output data arrays
%   will be returned with the same ordering as the input data arrays.
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Copyrights: (c) 2005 The Trustees of Princeton University and Board of
%                 Regents of the University of Texas.  All rights reserved.
%             (c) 2009 Kevin T. Chu.  All rights reserved.
% Revision:   $Revision: 149 $
% Modified:   $Date: 2009-01-18 00:31:09 -0800 (Sun, 18 Jan 2009) $
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
