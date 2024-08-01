% This functions approximates a NURBS curve by computing enough number 
% of points on the ideal curve using control points (cPts), weights,
% knot vector, and the shape (or basis) function values.
%
% INPUTS:
%   knotVector - A normalised real-valued knot vector. 
%   shapeFuncsDg - 
%   weights - A real-valued vector range between 0.0 and 1.0.
%   cPts - A real-valued 3-d vector containing all the necessary control points.
% OUTPUTS:
%   bspln - The computed points on the B-spline curve.
%   usfv - A reduced 2d matrix of shape function values without zero rows.
%   uksi - A reduced 2d matrix of knot values without zero rows.
% Example:
%   [nurbsCurve, usfv, uksi] = getnurbscrv(knotVector, shapeFuncDg,...
%   weights, cPts);
function [nurbsCurve, usfv, uksi] = getnurbscrv(knotVector, shapeFuncDg,...
    weights, cPts)
	[bspSFVs, ksi] = computebspbfunctions(knotVector, shapeFuncDg); % Computes the basis functions.	
    [shpeFuncs, ~, usfv, uksi] = pairrngdmin(bspSFVs, ksi,...
    knotVector, shapeFuncDg);   % pairs the domain (knot values) with range (shapefunctions).
	denominator = weights * shpeFuncs;
    rSFVs = weights'.* shpeFuncs./denominator;  
    if size(cPts, 3) ~= size(shpeFuncs, 1)  % Checks if the number of control points is correct.
        error('NO NURBS CURVE FOR YOU!')
    else
    nurbsCurve = zeros(size(cPts, 1), size(shpeFuncs, 2));  % Initiliases the 2-d vector of the points on the curve.
    for ii = 1:size(cPts, 3)
        nurbsCurve = nurbsCurve + cPts(:, :, ii) .* rSFVs(ii, :);   % Estimates the NURBS curve.
    end
        for ii = 1:size(nurbsCurve, 1)
            for jj = 1:size(nurbsCurve, 2)
                if nurbsCurve(ii, jj) == 0 
                    nurbsCurve(ii, jj) = NaN;   % Substitutes all zeros in the vector with NaNs.
                end
            end
        end
    end
end
