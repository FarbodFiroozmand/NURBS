% This script plots a random NURBS curve of degree 3 with a knot vector
% of knotVector = [0,0,0,0.33,0.66,1,1,1] and random control points.
clc; clear; close all;
%% Inputs
knotVector = [0,0,0,0.33,0.66,1,1,1];
shapeFuncDg = 3;
num = length(knotVector) - 1 - shapeFuncDg;
cps = zeros(2, 1, num);
rng("shuffle")
weights = rand(1, num);
for jj=1:num
    cps(:,:,jj) = 200 * rand(2,1) - 100;
end
%% Plots Nurbs Curves
figure(1)
for jj=1:num
    cps(:,:,jj) = 200 * rand(2,1) - 100;
end
[curve, usfv, uksi] = getnurbscrv(knotVector, shapeFuncDg, weights, cps);
cps = reshape(permute(cps, [2, 1, 3]), size(cps, 1)*size(cps, 2), size(cps, 3));
plot(cps(1,:), cps(2,:),"--o" , "color", "k",...
    "MarkerFaceColor", "k", 'LineWidth', 1.0)
hold on
plot(curve(1, :), curve(2, :), "k", 'LineWidth',1.5)
xlabel('\fontname{Courier}\fontsize{12} x') 
ylabel('\fontname{Courier}\fontsize{12} y')
subtitle({'\fontname{Courier}\fontsize{14} Polynomial Degree = 3'; ...
    '\fontname{Courier}\fontsize{14}\xi = [0,0,0,0.33,0.66,1,1,1]'})
title('\fontname{Courier}\fontsize{16} A Random NURBS Curve')
