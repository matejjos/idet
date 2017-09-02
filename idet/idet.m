function [iD] = idet (iA, precondition, method)
%BEGINDOC=================================================================
% .Author.
%
%  Josef Matejka
%
%-------------------------------------------------------------------------
% .Description.
%
%   Get the interval enclosure of determinant.
%
%
%-------------------------------------------------------------------------
% .Input parameters.
%   
%   iA ... an interval matrix
%   precondition ... a wanted precondition
%   method ... a wanted method
%
%------------------------------------------------------------------------
% .Output parameters.
%
%   iD ... interval enclosure of determinant
%
%------------------------------------------------------------------------
% .Implementation details. 
%
%  If no other parameters given apart from iA then igauss with inverse
%  precondition is used for general matrix, for symmetric is used method
%  based on enclosures of eigenvalues.
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
  % check input
  if nargin < 1
    iA = imatrand(5, 5, -5, 1);
  end
  sym = (issymmetric(sup(iA)) && issymmetric(inf(iA)));
  if nargin < 2
    if sym
      precondition = 'none';
    else
      precondition = 'inverse';
    end
  end
  if nargin < 3
    if sym
      method = 'sym';
    else
      method = 'igauss';
    end
  end
  denom = 1;
  % preconditioning
  switch precondition
    case 'inverse'
      C = inv(mid(iA));
      iA = C * iA;
      denom = vdet(C);
    case 'hansen'
      [C, ~] = lu(mid(iA));
      C = inv(C);
      iA = C * iA;
    case 'perturbation'
      row = randi(size(iA, 1));
      C = mid(iA);
      C(row, :) = rand(1, size(iA, 1)) .* (sup(iA(row, :)) - inf(iA(row, :)));
      C(row, :) = C(row, :) + inf(iA(row, :));
    case 'ldl'
      [C, ~] = ldl(mid(iA));
      C = inv(C);
      iA = C * iA * C';
    otherwise
  end
  
  switch method
    case 'gerschgoring'
      iD = idetgersch(iA) / denom;
    case 'hadamard'
      iD = idethad(iA) / denom;
    case 'cramer'
      iD = idetcram(iA) / denom;
    case 'sym'
      iD = prod(eigsymenc(iA)) / denom;
    otherwise
      iD = idetgauss(iA) / denom;
  end
endfunction