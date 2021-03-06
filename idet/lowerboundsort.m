function [isorted] = lowerboundsort (iv)
%BEGINDOC=================================================================
% .Author.
%
%  Josef Matejka
%
%-------------------------------------------------------------------------
% .Description.
%
%   Quicksort intervals by their lowerbound.s
%
%
%-------------------------------------------------------------------------
% .Input parameters.
%   
%   iv ...  array of intervals
%
%------------------------------------------------------------------------
% .Output parameters.
%
%   isorted ... sorted intervals
%
%------------------------------------------------------------------------
% .Implementation details. 
%
%   We call it recursively, as pivot is choosen the first one.
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
%   2017-07-10    first version
%
%------------------------------------------------------------------------
% .Todo.
%
%
%ENDDOC===================================================================
  if size(iv, 1) <= 1
    isorted = iv;
    return
  end
  pivot = iv(1);
  iv = iv(1:end);
  isorted = [lowerboundsort(iv(inf(iv) <= pivot)), pivot,
    lowerboundsort(iv(inf(iv) > pivot))];
endfunction