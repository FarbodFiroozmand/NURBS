% This functions approximates a NURBS surface by computing sufficient number 
% of points on the ideal surface using control points (cPts), weights,
% knot vectors, and the shape (or basis) function values.
%
% INPUTS: 
%   cPts - This parameter is a real-valued 3-d vector containing all the necessary control points.
%   weights - This parameter is a real-valued vector ranging between 0.0 and 1.0.
%   nShapeFuncDg/mShapeFuncDg - This parameter is a real-valued vector containg the shape function values.
%   knotVectorKsi/knotVectorEta - A normalised real-valued knot vector.
% OUTPUTS:
%   srf - The computed points on the NURBS surface.
%   nUSFV/mUSFV - A reduced 2d matrix of shape function values without zero rows.
%   uksi/ueta - A reduced 2d matrix of knot values without zero rows.
% Example:
%   [srf, nUSFV, mUSFV, uksi, ueta] = getnurbsrf(cPts, weights,...
%       nShapeFuncDg, mShapeFuncDg, knotVectorKsi, knotVectorEta)
function [srf, nUSFV, mUSFV, uksi, ueta] = getnurbsrf(cPts, weights,...
    nShapeFuncDg, mShapeFuncDg, knotVectorKsi, knotVectorEta)
    [nSFV, ksi] = computebspbfunctions(knotVectorKsi, nShapeFuncDg);   % Computes the basis functions.
    [mSFV, eta] = computebspbfunctions(knotVectorEta, mShapeFuncDg);
    [nShpeFuncs, ~, nUSFV, uksi] = pairrngdmin(nSFV, ksi, knotVectorKsi,...
        nShapeFuncDg);  % pairs the domain (knot values) with range (shapefunctions).
    [mShpeFuncs, ~, mUSFV, ueta] = pairrngdmin(mSFV, eta, knotVectorEta,...
        mShapeFuncDg);
    denominator = nShpeFuncs' * weights * mShpeFuncs;
    if size(cPts, 1) + size(cPts, 2) ~= size(nShpeFuncs, 1) + size(mShpeFuncs, 1)   % Checks if the number of control points is correct.
        error('NO B-SPLINE SURFACE FOR YOU!')
    else
        srf = zeros(size(nShpeFuncs, 2), size(mShpeFuncs, 2), size(cPts, 3));   % Initiliases the 3-d vector of the points on the curve.
        for ii = 1:size(cPts, 3)
           srf(:,:,ii) = nShpeFuncs' * (weights .* cPts(:, :, ii)) *...
               mShpeFuncs ./denominator;        % Estimates the NURBS surface.
        end
        for ii = 1:size(srf, 1)
            for jj = 1:size(srf, 2)
                if srf(ii, jj, :) == 0 
                    srf(ii, jj, :) = NaN;    % Substitutes all zeros in the vector with NaNs.
                end
            end
        end
    end
end
