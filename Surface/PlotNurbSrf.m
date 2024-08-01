% This script plots a random NURBS surface of degree 2 along
% the ksi-axis and degree 1 along the eta-axis, using knot vectors
% [0,0,0,0.5,1,1,1] and [0,0,0,1,1,1], and random control points. 
clc; clear; close all;
%% Inputs
knotVectorKsi = [0,0,0,0.5,1,1,1];
knotVectorEta = [0,0,1,1];
nShapeFuncDg = 2;
mShapeFuncDg = 1;
%% Plots NURBS Surface
nNum = length(knotVectorKsi) - 1 - nShapeFuncDg;
mNum = length(knotVectorEta) - 1 - mShapeFuncDg;
rng("shuffle")
weights = rand(nNum, mNum);
cPs = zeros(nNum, mNum, 3);
for ii=1:nNum
    for jj = 1:mNum
        cPs(ii, jj, :) = rand(3,1);
    end
end
[srf, nUSFV, mUSFV, uksi, ueta] = getnurbsrf(cPs, weights, nShapeFuncDg, mShapeFuncDg,...
    knotVectorKsi, knotVectorEta);
x=(nonzeros(srf(:,:,1)));
y=(nonzeros(srf(:,:,2)));
z=(nonzeros(srf(:,:,3)));

figure(1)
plot3(x,y,z,'Color','#48494B')
for ii = 1:size(cPs, 1)
    hold on
    plot3(cPs(ii,:,1), cPs(ii,:,2), cPs(ii,:,3), '--.', 'color', 'k',...
        'LineWidth', 0.75)
end
for ii = 1:size(cPs, 2)
    hold on
    plot3(cPs(:,ii,1), cPs(:,ii,2), cPs(:,ii,3), '--.', 'color', 'k',...
        'LineWidth', 0.75)
end
grid on
xlabel('\fontname{Courier}\fontsize{12} x') 
ylabel('\fontname{Courier}\fontsize{12} y')
zlabel('\fontname{Courier}\fontsize{12} z')
subtitle({'\fontname{Courier}\fontsize{14} Polynomial Degrees, \xi and \eta = 2 and 1'; ...
    '\fontname{Courier}\fontsize{14}\xi = [0,0,0,0.5,1,1,1], \eta = [0,0,0,1,1,1]'})
title('\fontname{Courier}\fontsize{16} A Random NURBS Surface')
