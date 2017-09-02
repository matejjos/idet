function ilambdas = vereigsym(A)


%BEGINDOC=================================================================
% .Author.
%
%  Jaroslav Horacek
%
%-------------------------------------------------------------------------
% .Description.
%
%   Computes verified enclosures of eigenvalues of real symmetric matrix A
%
%-------------------------------------------------------------------------
% .Input parameters.
%   
%   A ... n x n real matrix
%
%------------------------------------------------------------------------
% .Output parameters.
%
%   ilambdas ...  enclosures of eigenvalues, eigenvalues are sorted in 
%                 increasing order
%
%------------------------------------------------------------------------
% .Implementation details. 
%
%  It uses Rump's algorithm described e.g. in 
%  Fast Verification for Respective Eigenvalues of Symmetric Matrix
%  Shinya Miyajima,Takeshi Ogita, and Shin’ichi Oishi 
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
%
%------------------------------------------------------------------------
% .Todo.
%
%  incorporate more tight verification in the paper mentioned above
%
%ENDDOC===================================================================

[V, L] = eig(A);

lambdas = sort(diag(L));
n = length(lambdas);


iA = infsup(A);
iL = infsup(L);
iV = infsup(V);

iR = iA * iV - iV * iL;
iG = eye(n) - transpose(iV)*iV;

delta = (sqrt( inorm(iR, 'col') * inorm(iR, 'inf') )) / (1 - inorm(iG, 'inf'));

ilambdas = midrad(lambdas, delta);


end