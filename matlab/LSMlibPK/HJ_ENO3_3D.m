%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% HJ_ENO3_3D() computes the third-order plus and minus HJ ENO
% approximation to grad(phi)
%
% Usage: [phi_x_plus, phi_y_plus, phi_z_plus, ...
%         phi_x_minus, phi_y_minus, phi_z_minus] = ...
%        HJ_ENO3_3D(phi, ghostcell_width, dX)
%
% Arguments:
% - phi:              function for which to compute plus and minus
%                       spatial derivatives
% - ghostcell_width:  number of ghostcells at boundary of
%                       computational domain
% - dX:               array containing the grid spacing
%                       in coordinate directions
%
% Return values:
% - phi_x_plus:       x-component of third-order, plus
%                       HJ ENO derivative
% - phi_y_plus:       y-component of third-order, plus
%                       HJ ENO derivative
% - phi_z_plus:       z-component of third-order, plus
%                       HJ ENO derivative
% - phi_x_minus:      x-component of third-order, minus
%                       HJ ENO derivative
% - phi_y_minus:      y-component of third-order, minus
%                       HJ ENO derivative
% - phi_z_minus:      z-component of third-order, minus
%                       HJ ENO derivative
%
% NOTES:
% - phi_x_plus, phi_y_plus, phi_z_plus, phi_x_minus, phi_y_minus and
%   phi_z_minus have the same ghostcell width as phi.
%
% - All data arrays are assumed to be in the order generated by the
%   MATLAB meshgrid() function.  That is, data corresponding to the
%   point (x_i,y_j,z_k) is stored at index (j,i,k).  The output data 
%   arrays will be returned with the same ordering as the input data 
%   arrays.
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
