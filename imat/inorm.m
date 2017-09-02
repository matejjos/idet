function [n, ver] = inorm (iA, type)
%BEGINDOC=================================================================
% .Author.
%
%  Jaroslav Horacek
%
%-------------------------------------------------------------------------
% .Description.
%
%   computes the norm of an interval matrix
%
%-------------------------------------------------------------------------
% .Input parameters.
%   
%   iA ...  interval matrix 
%   type ... string characterizing norm type
%       1 - column norm (default) 
%       2 - Euclid norm
%       'inf' - row norm
%       'col' - column norm
%
%------------------------------------------------------------------------
% .Output parameters.
%
%   mn ...  matrix norm according to its type
%   ver ... 1 if the result is verified, 0 otherwise
%
%------------------------------------------------------------------------
% .Implementation details. 
%
%
%------------------------------------------------------------------------
% .Licence.
%  
%   Copyright (C) 2017  Charles University in Prague, Czech Republic
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
%   2017-04-06    first version
%   2017-04-06    max and frobenius norm added
%
%------------------------------------------------------------------------
% .Todo.
%   
%   Finish 'fro', 'u' norms, Hamming norm
%
%ENDDOC===================================================================

if nargin < 2
  type = 1;
end


switch(type)

  case {1, 'col'}
    mags = mag(iA);
    sums = sum(mags);
    n = max(sums);   
    
  case {Inf, 'inf'}
    mags = mag(transpose(iA));
    sums = sum(mags);
    n = max(sums);  
    
  case 'max'
    mags = mag(iA);  
    n = max(max(mags));
      
    
  otherwise 
    error("interval norm: wrong parameter type value");

endswitch

ver = 1;


end