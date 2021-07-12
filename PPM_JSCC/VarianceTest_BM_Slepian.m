 clear all; clc;

nIters = 1e4;
ENR = -5:1:25; 
ENRlin = 10.^(ENR/10);
MaxVarBM = zeros(size(ENR)); 
MaxVarSlep = zeros(size(ENR)); 
dx = 1e-3; beta = 5; 
xMax = 25;
x = (1/beta):dx:xMax;

for i=1:length(ENR)
   currMaxBM = zeros(1,nIters); currMaxSlep = zeros(1,nIters);
    parfor iter=1:nIters
        
        n = randn(1,length(x));
        
        % create BM Minus Parabols
        BM = cumsum(n) - sqrt(1/(8*ENRlin(i)))*(x.^2);
        Slepian = filter(ones(1,beta/dx)/sqrt(beta/dx),1,n) - sqrt(1/(8*ENRlin(i)))*(x.^2);
%         figure;hold all;
%         plot(BM);plot(Slepian);legend('BM','Slepian');
        % finx Max
        [maxBM,maxIdxBm] = max(BM);
        [maxSlep,maxIdxSlep] = max(Slepian);
        
        currMaxBM(iter) = maxIdxBm*dx; currMaxSlep(iter) = maxIdxSlep*dx; 
    end
    
    MaxVarBM(i) = std(currMaxBM);
    MaxVarSlep(i) = std(currMaxSlep);
end

figure;hold all;plot(ENR,MaxVarBM,'LineWidth',2);
plot(ENR,MaxVarSlep,'--','LineWidth',2);
grid on; grid minor;
xlabel('ENR [dB]'); legend('Var BM','Var Slepian');
title(strcat('Variance for \beta = ',num2str(beta)));