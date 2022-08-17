function tf = valVarCell(c,v)
% VALVARCELL validates c with validation function v for both cell and
% non-cell inputs of c. If c is a cell, v validates all elements of c
% individually. If c is not a cell, v validates c directly.
%
% tf = VALVARCELL(v, c)
%   Validate c directly or the elements of c using validation function v,
%   depending on whether c is a cell or not.

tf = (iscell(c) && all(cellfun(@(p) v(p),c))) || v(c);