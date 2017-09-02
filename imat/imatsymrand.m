function iA = imatsymrand (n, radexp, density, midbounds)
%BEGINDOC=================================================================
% .Author.
%
%  Josef Matejka
%
%-------------------------------------------------------------------------
% .Description.
%
%   Generated a symmetrix interval matrix.
%
%-------------------------------------------------------------------------
% .Input parameters.
%   
%   n  ... size of n x n matrix
%   radexp ... the bound on interval rad is 10^radexp
%   density ... how many non-zero columns should matrix have
%   midbounds ... 2-array - [a,b] a <= b interval for random midpoints 
%                 if a=b  we get the same midpoints everywhere
%
%------------------------------------------------------------------------
% .Output parameters.
%
%   iA ... symmetric interval matrix
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
%   2017-03-30    first version
%   2017-04-08    renamed
%
%------------------------------------------------------------------------
% .Todo.
%
%
%ENDDOC===================================================================
  %% input checking
  if nargin < 1
    n = 3;
  end
  
  if nargin < 2
    radexp = 0;
  end
  
  if nargin < 3
    density = 1;
  end
  
  if nargin < 4
    midbounds = [-20 20];
  end
  
  %% input validation
  validateattributes(radexp, {'numeric'}, {'integer'});
  validateattributes(n, {'numeric'}, {'integer', 'positive'});
  validateattributes(density, {'numeric'},
    {'positive', '>=', 0, '<=', 1});
  validateattributes(midbounds,  {'numeric'},
    {'size', [1 2]});
  
  if midbounds(1) > midbounds(2)
     error("Wrong [a,b] bounds for midpoints a >= b");   
  end
  
  %% number of distinct values
  s = (n^2 + n) / 2;
  
  Amid = zeros(n);
  Arad = zeros(n);
  %% vector of those distinct values
  amid = randi((midbounds(2) - midbounds(1)) + 1, s, 1);
  amid = amid + midbounds(1) - 1;
  %% vector of radii
  arad = rand(s, 1) * 10^radexp;
  
  k = 1;
  for i = 1:n
    for j = i:n
      Amid(i, j) = Amid(i, j) + amid(k);
      Arad(i, j) = Arad(i, j) + arad(k);
      if i != j
        Arad(j, i) = Arad(j, i) + arad(k);
        Amid(j, i) = Amid(j, i) + amid(k);
      end
      k = k + 1;
    end
  end
  
  A = full(sprandsym(n, density) != 0);
  iA = midrad(Amid .* A, Arad .* A);
endfunction