function [iD] = idethad (iA)
%BEGINDOC=================================================================
% .Author.
%
%  Josef Matejka
%
%-------------------------------------------------------------------------
% .Description.
%
%   Get the interval enclosure using Hadamard's inequality.
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
  iD = infsupdec(0);
  v = zeros(size(iA, 1), 1);
  for i = 1:size(iA, 1)
    q = 0;
    for j = 1:size(iA, 2)
      q = q + mag(iA(i, j))^2;
    end
    v(i) = sqrt(q);
  end
  v = prod(v);
  iD = infsupdec(-v, v);
endfunction