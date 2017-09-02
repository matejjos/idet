function iA = imatrand (m, n, radexp, density, midbounds)
%BEGINDOC=================================================================
% .Author.
%
%  Josef Matejka
%
%-------------------------------------------------------------------------
% .Description.
%
%   This function will generate random sparse interval matrix, with
%   given width, density and sizes.
%
%-------------------------------------------------------------------------
% .Input parameters.
%   
%   m ... number of rows
%   n ... number of columns
%   radexp ... the bound on interval rad is 10^radexp
%   density ... how many non-zero columns should matrix have
%   midbounds ... 2-array - [a,b] a <= b interval for random midpoints 
%                 if a=b  we get the same midpoints everywhere
%
%------------------------------------------------------------------------
% .Output parameters.
%
%   iA ... a sparse interval matrix
% 
%------------------------------------------------------------------------
% .Implementation details. 
%
%
%------------------------------------------------------------------------
% .License.
%    Copyright (C) 2017  Charles University in Prague, Czech Republic
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
%   2017-03-14    first version
%   2017-03-22    new name, second version
%   2017-04-04    renamed, switch of parameters m,n
%   2017-04-04    added new parameter midbounds
%
%------------------------------------------------------------------------
% .Todo.
%
%
%ENDDOC===================================================================
  %% input checking
  if nargin < 1
     m = 3;
  end   
     
  if nargin < 2
    n = m;
  end

  if nargin < 3
    radexp = 0;
  end
  
  if nargin < 4
    density = 1;
  end
  
  if nargin < 5
    midbounds = [-20, 20];
  end
  
  %% input validation
  validateattributes(m, {'numeric'}, {'integer', 'positive'});
  validateattributes(n, {'numeric'}, {'integer', 'positive'});
  validateattributes(radexp, {'numeric'}, {'integer'});
  validateattributes(density, {'numeric'},
    {'positive', '>=', 0, '<=', 1});
  validateattributes(midbounds,  {'numeric'},
    {'size', [1 2]});

    
  if midbounds(1) > midbounds(2)
     error("Wrong [a,b] bounds for midpoints a >= b");   
  end
  
  
  % computing sparse matrix with given density (Octave functions)
  A = full(sprand(m, n, density) > 0);
  
  % computing midpoint Amid
  
  if(midbounds(1) == midbounds(2))
    Arand = repmat(midbounds(1), m, n);
  else 
    Arand = randi((midbounds(2) - midbounds(1)) + 1,m, n);
    Arand = Arand + midbounds(1) - 1;  
  end
  
  Amid = A .* Arand;
  
  % computing matrix of radii
  Arad = rand(m, n) .* (10^radexp);
  Arad = A .* Arad;
  
  iA = midrad(Amid, Arad);
  
end
