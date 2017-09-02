function [ ienc, largest ] = eigsymencrohn( iA, iEigsAc, ...
    iEigsAd, spectral)
%BEGINDOC==================================================================
% .Author.
%
%  Nikolay Stoyanov Kaleyski, Jaroslav Horacek
%
%--------------------------------------------------------------------------
% .Description.
%
%  Very fast outer approximation of the eigenvalue interval bounds for a
%  real symmetric interval matrix
%
%--------------------------------------------------------------------------
% .Input parameters.
%
%  iA ........ a real symmetric interval matrix
%  iEigsAc ... optionally, a precomputed vector of eigenvalue bounds of the
%  central matrix of iA (must be verified)
%  iEigsAd ... similarly, a precomputed vector for the radial matrix of iA;
%  needs to be verified
%  spectral .. similarly, a precomputed spectral radius of the radial
%  matrix of iA; needs to be verified
%
%--------------------------------------------------------------------------
% .Output parameters.
%
%  ienc ... an interval vector containing an outer approximation of the 
%  eigenvalues of iA
%  largest ... the upper bound for the largest eigenvalue of iA
%
%--------------------------------------------------------------------------
% .Implementation details.
%
%  The function computes its results by trivially following the statement
%  of the theorem. Both sets of eigenvalues, 'eigsAc' and 'eigsAd' must al-
%  ways be computed, even if we're looking only for the largest eigen-
%  value's upper bound, and the computation of the upper and lower bounds
%  takes very little time (it only involves adding and substracting a given
%  vector from another), so it is always performed.
%
%  If the function is called repeatedly from a more complex algorithm,
%  chances are that the verified eigenvalues of mid(iA) and rad(iA) will
%  need to be computed for several different functions; therefore, we
%  allow the caller to provid pre-computed results so that we avoid
%  redunadant computations.
%
%  method vereigsym returns eigenvalues sorted l1 <= l2 <= ...
%
%--------------------------------------------------------------------------
% .License
%
%  [License goes here]
%
%--------------------------------------------------------------------------
% .History
%
%  2013-11-24   first version, adapted from the previous function "Rohn.m"
%  2014-02-20   added support for precomputed data
%  2017-03-30   linked to our own implemented method vereigsym
%  2017-03-30   renamed to eigsymencrohn
%
%--------------------------------------------------------------------------
% .Todo.
%
%
%ENDDOC====================================================================

%Compute eigenvalues of mid(iA) if none are provided
if nargin < 2
    iEigsAc = vereigsym(mid(iA));
end

%Compute eigenvalues of rad(iA) if none are provided
if nargin < 3
    iEigsAd = vereigsym(rad(iA));
end

%Compute the spectral radius of rad(iA) if not provided
if nargin < 4
    spectral = max( mag(iEigsAd) );
end

%Apply the theorem
ienc = infsup( inf(iEigsAc) - spectral, sup(iEigsAc) + spectral );
[n,~] = size(iA);
%Bound on largest eigenvalue as a separate output parameter for convenience
largest = sup(ienc(n,1));

end

