function [iD] = idetcram (iA)
%BEGINDOC=================================================================
% .Author.
%
%  Josef Matejka
%
%-------------------------------------------------------------------------
% .Description.
%
%   Get the interval enclosure using Cramer's rule.
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
%   Implemented without recursion.
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
  iD = infsupdec(zeros(1, size(iA, 1)));
  b = zeros(size(iA, 1), 1);
  iC = iA;
  for i = 1:(size(iA, 1) - 1)
    b(i) = 1;
    iC(:, i) = b;
    b(i) = 0;
  end
  iD(size(iA, 1)) = iC(size(iC, 1), size(iC, 2));
  for i = (size(iA, 1) -1):-1:1
    b(i) = 1;
    iC(:, i) = iA(:, i);
    try
      x = iC\b;
    catch
      x = zeros(size(iA, 1), 1);
    end
    if(x(i) == 0)
      break;
    end
    iD(i) = iD(i + 1) / x(i);
    b(i) = 0;
  end
  iD = iD(1);
  if(x(i) == 0)
    iD = entire;
  end
endfunction