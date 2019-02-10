function [Lqij,nextL] = statUpdate(numOfChuncks,Ms,newL,Lqij_init)
Lqij = Lqij_init;
nextL = [];

for k = 1:numOfChuncks - 1
    
    % Future statistics section
    nextL = [zeros(2*(Ms + 1) - 2, Ms + 1) newL(1:end - 2, 1:Ms, k + 1)];
    
    % Current statistics
    Lqij(:, :, k) = [newL(:, :, k); nextL];
    
end % for k
end

