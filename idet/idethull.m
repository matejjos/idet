function id = idethull(iA)
%BEGINDOC=================================================================
% .Author.
%
%  Jaroslav Horacek
%
%-------------------------------------------------------------------------
% .Description.
%
%   computes hull of interval determinant of matrix by brute force
%
%-------------------------------------------------------------------------
% .Input parameters.
%   
%   iA ...  square interval or real matrix 
%
%------------------------------------------------------------------------
% .Output parameters.
%
%   id ...  determinant
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
%  Verified real det computation 
%
%ENDDOC===================================================================


[m, n] = size(iA);
if m != n
  error "idethull: Matrix must be of a square shape";
end

%counting number of nonzero elements + creating bitmap
bitmap = zeros(n,n);
c = 0;
for i=1:n
  for j=1:n
    if (inf(iA(i,j)) < 0 || sup(iA(i,j)) > 0)
      c = c + 1;
      bitmap(i, j) = 1;
    end   
  end
end




gen = zeros(1, c);
detmin = Inf;
detmax = -Inf;

while any(gen == 0)

  % mapping gen into M
  dummyc = 1;
  M = zeros(n,n);
  for i=1:n
    for j=1:n
      if (bitmap(i,j) == 1)
        M(i, j) = gen(dummyc);
        dummyc = dummyc + 1;
      end   
    end
  end
  
  

  % generating 0/1 matrix from gen vector 
  U = M .* sup(iA); 
  M = M - 1;
  M = abs(M);
  L = M .* inf(iA);
 
  % computing determinant of real matrix 
  % !!! verification
  edgeM = U + L; 
  d = det(edgeM);
  
  if(d > detmax)
    detmax = d;
  end
  
  if(d < detmin)
    detmin = d;
  end
  
  % generating next 0/1 vector
  i = c;
  while gen(i) == 1 
  i = i - 1;  
  end
  
  if i >= 1
    gen(i) = 1;   
  end
  
  for j=i+1:c
    gen(j) = 0;
  end
 
  
end 

  id = infsup(detmin, detmax);

end




























