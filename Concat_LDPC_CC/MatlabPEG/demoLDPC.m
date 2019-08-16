% demo for generating regular PEG LDPCs
% Date: May 10 2013
% Hatef Monajemi (monajemi@stanford.edu)

mex GenerateLDPC_mex.C BigGirth.C CyclesOfGraph.C Random.C GenerateLDPC.C
n = 252;
N = 504;
dc= 5;
A = buildLDPCmatrix(n,N,dc);

figure;
spy(A,7);
title(['PEG-LDPC: N = ',num2str(N), ...
        ', n = ', num2str(n), ...
        ', d = ', num2str(dc)], 'fontsize', 18)
set(gca, 'fontSize', 17)
