clear all; clc;

Delta = 1;
nPoints = 5e3;
dither = generate2D_dither(Delta,nPoints);

scatterplot(dither(:,1) + 1j*dither(:,2))
hist3(dither,'CDataMode','auto','FaceColor','interp')
title(strcat('Npoints = ',num2str(nPoints)));