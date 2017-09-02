function [ ib ] = eigsymencindirect( iA, midEigs, upperboundmethod )
%BEGINDOC==================================================================
% .Author.
%
%  Nikolay Stoyanov Kaleyski
%
%--------------------------------------------------------------------------
% .Description.
%
%  Indirect interlacing method for computing an outer approximation for the
%  bounds of the eigenvalue intervals of a real symmetric interval matrix.
%
%--------------------------------------------------------------------------
% .Input parameters.
%
%  iA ... a real symmetric interval matrix
%  midEigs ... precomputd verified eigenvalue intervals for the central
%  matrix of iA; will be computed with vereig if not supplied
%  upperboundmethod ... optionally, a method for computing an outer bound
%  for all eigenvalues of iA; a version of the direct interlacing algorithm
%  is used by default
%
%--------------------------------------------------------------------------
% .Output parameters.
%
%  ib ... an interval vector with the computed approximation
%
%--------------------------------------------------------------------------
% .Implementation details.
%
%  This is a wrapper function which utilizes the implementing function
%  eigsymencouterindirectsingleupper.m for computing bounds for individual
%  eigenvalues and then combining them into a complete approximation.
%
%   More details in Bounds on Real Eigenvalues and Singular Values of 
%   Interval Matrices, Milan Hladík, David Daney, Elis Tsigaridas, 2010s 
%
%--------------------------------------------------------------------------
% .License
%
%  [license goes here]
%
%--------------------------------------------------------------------------
% .History.
%
%  2014-07-08   first version
%  2017-03-30   renamed eigsymencindirect
%  2017-03-30   linked with our own vereigsym
%
%--------------------------------------------------------------------------
% .Todo.
%
%
%ENDDOC====================================================================

%Default values for input parameters
if nargin < 2
    midEigs = vereigsym(mid(iA));
end

if nargin < 3
    upperboundmethod = @eigsymencouterindirectdefaultbound;
end

%Size of the input matrix
[n,~] = size(iA);

minEigs = flipud(-midEigs); %central matrix of -iA
radEigs = upperboundmethod(infsup(-rad(iA),rad(iA))); %bounds for the
%matrix infsup(-rad(iA),rad(iA)); we do not provide a possibility for
%passing precomputed approximations for it since it most probably will
%not have been precomputed

ib = infsup(-inf(n,1),inf(n,1));

for i = 1:n
    up = eigsymencouterindirectsingleupper(midEigs, radEigs, n-i+1);
    down = eigsymencouterindirectsingleupper(minEigs, radEigs, i);
    
    ib(i,1) = infsup(-down,up);
end

end

