function idet = vdet (A)
%BEGINDOC=================================================================
% .Author.
%
%  Josef Matejka
%
%-------------------------------------------------------------------------
% .Description.
%
%   Counts verified determinant of a non-interval matrix A.
%
%-------------------------------------------------------------------------
% .Input parameters.
%   
%   A ... an non-interval matrix
%
%------------------------------------------------------------------------
% .Output parameters.
%
%   idet ... lower and upper bound on determinat
%
%------------------------------------------------------------------------
% .Implementation details. 
%
%   It's a function copied from Rump's INTLAB, modified for OCTAVE.
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
%   2017-03-30    first version
%
%------------------------------------------------------------------------
% .Todo.
%
%  Add verified matrix multiplication, verified LU inverses to get
%  better results.
%
%ENDDOC===================================================================
  if nargin != 1 
    error("Error! Number of input parameters: " + nargin + ", but only 1 is supported.")
  end
  validateattributes(A, {'numeric'},{'square'});
  [L, U, p] = lu(A, 'vector');
  I = speye(size(A));
  XL = I / L;
  XU = U \ I;
  iB = XL*infsupdec(A(p, :))*XU;
  c = mid(diag(iB));
  r = mag(sum(iB - diag(c), 2));
  ig = midrad(c, r);
  iu = prod(infsupdec(diag(XU)));
  idet = det(I(p, :))*prod(ig)/iu;
endfunction