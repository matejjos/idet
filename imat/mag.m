function m = mag(iA)
%BEGINDOC=================================================================
% .Author.
%
%  Jaroslav Horacek
%
%-------------------------------------------------------------------------
% .Description.
%
%   computes the magnitude of an interval matrix (if iA is real returns abs)
%
%-------------------------------------------------------------------------
% .Input parameters.
%   
%   iA ...  interval or real matrix 
%
%------------------------------------------------------------------------
% .Output parameters.
%
%   m ...  matrix magnitude 
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
%
%------------------------------------------------------------------------
% .Todo.
%   
%
%ENDDOC===================================================================


  if isa(iA, 'infsup')
    m = mag(iA);
  else
    m = abs(iA);
  end

end