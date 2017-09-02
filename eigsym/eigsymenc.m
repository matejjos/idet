function ienc = eigsymenc (iA) 


%BEGINDOC=================================================================
% .Author.
%
%  Jaroslav Horacek
%
%-------------------------------------------------------------------------
% .Description.
%
%   Computes enclosure of eigenvalues of iA symmetric n x n matrix
%
%-------------------------------------------------------------------------
% .Input parameters.
%   
%   iA ... n x n interval matrix
%
%------------------------------------------------------------------------
% .Output parameters.
%
%   ienc ... enclosures of eigenvalues of iA
%
%------------------------------------------------------------------------
% .Implementation details. 
%
%------------------------------------------------------------------------
%  .Licence.
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
%   2017-03-30    first version
%   2017-07-10    due to renaming of some method, we had to change them
%                 here as well.
%
%------------------------------------------------------------------------
% .Todo.
%
%ENDDOC===================================================================


ienc1 = eigsymencdirect(iA, 'both');
ienc2 = eigsymencindirect(iA);

ienc = infsup(zeros(size(ienc1, 1), 1));
for i = 1:size(ienc1, 1)
  ienc(i) = intersect(ienc1(i), ienc2(i));
end

end