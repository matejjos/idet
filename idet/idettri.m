function id = idettri(iA)
%BEGINDOC=================================================================
% .Author.
%
%  Jaroslav Horacek
%
%-------------------------------------------------------------------------
% .Description.
%
%   computes hull of interval determinant of a tridiagonal matrix in P
%
%-------------------------------------------------------------------------
% .Input parameters.
%   
%   iA ...  square interval or real matrix 
%
%------------------------------------------------------------------------
% .Output parameters.
%
%   id ... determinant
%
%------------------------------------------------------------------------
% .Implementation details. 
%
%   Brute force according to Edge theorem by Rohn is used, therefore it may 
%   take a lot of time  
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
%   2017-06-26    first version
%
%------------------------------------------------------------------------
% .Todo.
%ENDDOC===================================================================


[m, n] = size(iA);
if m != n
  error "idethull: Matrix must be of a square shape";
end


minors = repmat(infsup(0,0), 1, n);
minors(1) = iA(1,1);

if n >= 2
  minors(2) = iA(1,1)*iA(2,2) - iA(1,2)*iA(2,1);
  
  
  for i=3:n
    
   
    d1 = iA(i,i)*minors(i-1) - iA(i, i-1)* iA(i-1, i) * sup(minors(i-2));
    d2 = iA(i,i)*minors(i-1) - iA(i, i-1)* iA(i-1, i) * inf(minors(i-2));
  
    list = [inf(d1), inf(d2), sup(d1), sup(d2)];
  
    minors(i) = infsup( min(list), max(list) );
   
   
   
  end
  
  end

id = minors(n);

end




























