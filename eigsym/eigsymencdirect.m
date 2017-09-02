function [ ib ] = eigsymencdirect( ...
    iA, ...
    mode, ...
    upperboundmethod, ...
    indexselectionmethod ...
)
    
%BEGINDOC==================================================================
% .Author.
%
%  Nikolay Stoyanov Kaleyski
%
%--------------------------------------------------------------------------
% .Description.
%
%  Computation of an outer approximation of the bounds of the eigenvalue
%  intervals of the real symmetric matrix iA using the direct interlacing
%  algorithm.
%
%--------------------------------------------------------------------------
% .Input parameters.
%
%  iA ... a real symmetric interval matrix
%  mode ... one of 'FORWARD', 'REVERSE' or 'BOTH', specifying which
%           variant of the algorithm to run; default is 'FORWARD'
%  upperboundmethod ... a function receiving a symmetric interval matrix
%  and returning an upper outer bound for its largest eigenvalue interval;
%  default is a quick estimate using Rohn's theorem and the eigenvalues of
%  mag(iA)
%  indexselectionmethod ... a method for selecting indices to add or
%  remove in accordance with the implementation of the direct interlacing
%  method; see the description of eigsymencouterdirectsingleupper.m 
%  for more details; default is the first such method described in the 
%  bachelor thesis (p. 19)
%
%--------------------------------------------------------------------------
% .Output parameters.
%
%  ib ... the computed outer approximation of the eigenvalue interval
%  bounds
%
%--------------------------------------------------------------------------
% .Implementation details.
%
%  This is essentially a wrapper method which calls the function
%  eigsymencouterdirectsingleupper.m for the individual eigenvalues of iA
%  and -iA (for upper and lower bounds). A vector is used for keeping track
%  of the indices which have been added to/removed from the matrix so as
%  not to have eigsymencouterdirectsingleupper.m account for more than one
%  index at once.
%
%   More details in Bounds on Real Eigenvalues and Singular Values of 
%   Interval Matrices, Milan Hladík, David Daney, Elis Tsigaridas, 2010s 
%
%--------------------------------------------------------------------------
% .License.
%
%  [License goes here]
%
%--------------------------------------------------------------------------
% .History.
%
%  2014-07-07   first version
%  2017-03-30   renamed to eigsymencdirect
%
%--------------------------------------------------------------------------
% .Todo.
%
%  
%ENDDOC====================================================================

%Default parameter values
if nargin < 2
    mode = 'FORWARD';
end

if nargin < 3
    upperboundmethod = @eigsymencouterdefaultupperbound;
end

if nargin < 4
    indexselectionmethod = @eigsymencouterdefaultindexselector;
end

[n,~] = size(iA); %dimension of the input matrix

ib = infsup( -inf(n,1), inf(n,1) ); %initial approximation is the "worst
%possible", so that "better" approximations are always found with the first
%direction used

%The function "interlacing" is used for improving the existing
%approximation by sequentially processing all eigenvalue indices

if strcmpi(mode, 'forward') || strcmpi(mode, 'both')
    ib = interlacing(ib,iA,'FORWARD',upperboundmethod, ...
        indexselectionmethod);
end

if strcmpi(mode, 'reverse') || strcmpi(mode, 'both')
    ib = interlacing(ib,iA,'REVERSE',upperboundmethod, ...
        indexselectionmethod);
end

end

%Calls the single upper bound method successively
function b = interlacing(b,iA,direction, ...
    upperboundmethod, indexselectionmethod)

    [n,~] = size(iA);

    %Selected indices (added or removed, depending on the direction) for
    %iA and -iA respectively:
    plusindices = [];
    minusindices = [];
    
    %The order in which the indices are processed depends on the
    %algorithm's direction; the iterationset vector is then used by the
    %following "for" loop
    if strcmpi(direction, 'reverse')
        iterationset = n:-1:1;
    else
        iterationset = 1:n;
    end
    
    for i = iterationset
        %Calling eigsymencouterdirectsingleupper returns the computed upper
        %bound along with an updated set of selected indices
        
        [up,plusindices] = eigsymencouterdirectsingleupper(iA, n-i+1, ...
            'STEP', b, direction, plusindices, upperboundmethod, ...
            indexselectionmethod);
    
        %neither the index
        %strategy, nor the existing approximation play any role here, as
        %we're calculating the eigenvalue approximations one by one; we
        %still need to specify some values
        
        [down,minusindices] = eigsymencouterdirectsingleupper(-iA, ...
            n-i+1, 'STEP', b, direction, minusindices, ...
            upperboundmethod, indexselectionmethod);
        
        %Select better values (between the newly computed ones and the
        %existing approximation)
        up = min(up, sup(b(n-i+1,1)));
        down = max(-down, inf(b(i,1)));
        
        %Possibly improve existing approximation
        b(n-i+1,1) = infsup(inf(b(n-i+1,1)), up);
        b(i,1) = infsup(down, sup(b(i,1)));
    end

end