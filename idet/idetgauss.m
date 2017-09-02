function [iD] = idetgauss (iA)
%BEGINDOC=================================================================
% .Author.
%
%  Josef Matejka
%
%-------------------------------------------------------------------------
% .Description.
%
%   Get the interval enclosure using Gauss Elimination.
%
%
%-------------------------------------------------------------------------
% .Input parameters.
%   
%   iA ... an interval matrix
%
%------------------------------------------------------------------------
% .Output parameters.
%
%   iD ... interval enclosure of determinant
%
%------------------------------------------------------------------------
% .Implementation details. 
%
%
%------------------------------------------------------------------------
% .License.
%    Copyright (C) 2011  Charles University in Prague, Czech Republic
%
%   LIME 1.0 is free for private use and for purely academic purposes. 
%   It would be very kind from the future user of LIME 1.0 to give 
%   reference that this software package has been developed 
%   by at Charles University, Czech Republic.
%
%   For any other use of LIME 1.0 a license is required.
%
%   THIS SOFTWARE IS PROVIDED AS IS AND WITHOUT ANY EXPRESS OR IMPLIED 
%   WARRANTIES, INCLUDING, WITHOUT LIMITATIONS, THE IMPLIED WARRANTIES 
%   OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE.   
%    
%-------------------------------------------------------------------------
% .History.
%       
%   2017-07-10    first version
%
%------------------------------------------------------------------------
% .Todo.
%
%
%ENDDOC===================================================================
  [iQ, ~, E, S]= igauss(iA);
  if strcmp(S, 'sin')
    iD = entire;
  else
    iD = E * prod(diag(iQ));
  end
endfunction