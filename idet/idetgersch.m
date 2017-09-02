function [iD] = idetgersch (iA)
%BEGINDOC=================================================================
% .Author.
%
%  Josef Matejka
%
%-------------------------------------------------------------------------
% .Description.
%
%   Get the interval enclosure using Gerchgoring's circles.
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
  il = infsupdec(zeros(size(iA, 1), 1));
  for i = 1:size(iA, 1)
    il(i, 1) = iA(i, i);
    iA(i, i) = 0;
    irad = sum(mag(iA(i, :))) + wid(il(i, 1));
    il(i, 1) = infsupdec(mid(il(i)) - irad, mid(il(i)) + irad);
  end
  il = lowerboundsort(il');
  il = il';
  cur = 1;
  for i = 2:size(iA, 1)
    if intersect(il(cur, 1), il(i, 1)) != empty
      il(cur, 1) = union(il(cur, 1), il(i, 1));
      il(i, 1) = il(cur, 1);
    else
      cur = i;
    end
  end
  iD = prod(il(:, 1));
endfunction