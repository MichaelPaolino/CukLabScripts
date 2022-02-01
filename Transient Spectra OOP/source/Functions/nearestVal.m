function [valsOut, ind] = nearestVal(valsIn,targetVals)
%Find the nearest values targetVals in vector vals. This is useful if a 
%user wants to select a subset of delays or wavelengths from their data. 
%The function takes the raw delays or wavelengths and returns the nearest 
%values and indicies to the target values.
%inputs:
%   valsIn: vector of values to search in for targetVals
%   targetVals: vector or values to search for in valsIn
%outputs:
%   valsOut: vector of unique values that were closest to targetVals
%   ind: indicies of valsOut, i.e. valsOut = valsIn(ind);

    %initialize new delay vectors
    nTargetVals = length(targetVals);
    ind = zeros(nTargetVals,1);
    
    %find nearest delay to user specified target delay and record the index
    for ii = 1:nTargetVals
        [~,ind(ii)] = min(abs(valsIn - targetVals(ii)));
    end

    %update delay values based on index and remove any duplicate values
    [valsOut, ia] = unique(valsIn(ind),'stable');
    ind = ind(ia);
end