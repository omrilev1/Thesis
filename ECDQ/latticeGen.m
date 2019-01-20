function [latticePoints] = latticeGen(G,nPoints)
% This function gets the generator matrix of 2D Lattice G , and outputs 
% n cell centers of the lattice
% the outputs x , y are vectors that contain the cordinates

latticePoints = [];

nPointsX = round(sqrt(nPoints));
nPointsY = round(sqrt(nPoints));

[i,j] = meshgrid(-nPointsX/2:1:nPointsX/2 , -nPointsY/2:1:nPointsY/2);

for idx1 = 1:size(i,1)
    for idx2 = 1:size(i,2)
        
        latticePoints = [latticePoints G(:,1) * i(idx1,idx2) + G(:,2) * j(idx1,idx2)];
        
    end
end
end

