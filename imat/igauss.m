function [iA, iB, d, state] = igauss (iA, iB)
%BEGINDOC=================================================================
% .Author.
%
%  Josef Matejka
%
%-------------------------------------------------------------------------
% .Description.
%
%   Gaussian elimination
%
%   Function provides Gaussian elimination on an interval matrix [A|B],
%   that represents system AX = B.
%
%-------------------------------------------------------------------------
% .Input parameters.
%   
%   iA ...  m x n matrix 
%   iB ...  n x k matrix - optional
%
%------------------------------------------------------------------------
% .Output parameters.
%
%   iA  ...  m x n upper diagonal matrix, eliminated A matrix.
%              on diagonal 
%   iB  ...  n x k matrix, if we suppose PA = C, then PB = D, where P is
%              a matrix of elimination steps.
%   d   ...  the sign of a determinant of given matrix
%   state ...  'ok' means no error occured, 'sin' means that A is either
%              singular, or it seems singular due to wideing of interval
%              arithmetic.
%
%------------------------------------------------------------------------
% .Implementation details. 
%
%   We choose pivot for every column as the biggest absolute value of
%   mid of an interval
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
%   2016-10-22    first version
%   2017-03-23    changed we things, code should be more readable
%   2017-04-07    added state
%
%------------------------------------------------------------------------
% .Todo.
%
%
%ENDDOC===================================================================
  % check input
  if nargin < 1
    iA = rand(3);
  end
  if nargin < 2
    iB = zeros(size(iA, 1), 1);
  end

  % initialization
  state = 'ok';
  iAB = horzcat(iA, iB);
  d = 1;
  k = 1;
  
  while ((k <= size(iAB, 1)) && (k < size(iA, 2)))
    for i = k:size(iAB, 1)
      if (not(ismember(0, iAB(i, k))) &&
        (abs(mid(iAB(i, k))) > abs(mid(iAB(k, k)))))
          iAB([i k], :) = iAB([k i], :);
          d = -d;
      end
    end
    
    if not(ismember(0, iAB(k, k)))
      for i = (k + 1):size(iAB, 1)
        alpha = iAB(i, k) / iAB(k, k);
        iAB(i, :) = iAB(i, :) - alpha * iAB(k, :);
        iAB(i, k) = 0;
      end
    else
      state = 'sin';
    end
    k = k + 1;
  end
  
  asize = size(iA);
  absize = size(iAB);
  iA = iAB(1:asize(1), 1:asize(2));
  iB = iAB(1:asize, (asize(2) + 1):absize(2));
endfunction