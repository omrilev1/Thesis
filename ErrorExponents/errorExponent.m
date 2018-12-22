%% This script plots the F_0 function , in the erasure decoding case for BSC and BEC Channels
%% ErrprExponent = -log(F(rho,s))

% Code
clear ; close all; clc;

s = [0.5];
rho = [0.5];

epsilon_bsc = 0.00005:0.00005:0.5;
epsilon_bec = 0.0001:0.0001:1;

F_bsc = zeros(length(rho),length(s),length(epsilon_bsc));
F_bec = zeros(length(rho),length(s),length(epsilon_bec));

for i=1:length(rho)
    for j=1:length(s)
        F_bec(i,j,:) = epsilon_bec + (2^(-1*rho(i)))*(1-epsilon_bec);
        
        F_bsc(i,j,:) =  (2^(-1*rho(i)))*(epsilon_bsc.^(1-s(j)) + (1-epsilon_bsc).^(1-s(j))).*((epsilon_bsc.^(s(j)/rho(i)) + (1-epsilon_bsc).^(s(j)/rho(i)))).^(rho(i));
        
    end
end

resultPlot(s,rho,epsilon_bec,epsilon_bsc,F_bec,F_bsc);

% Simple Functions
function [] = resultPlot(sToPlot,rhoToPlot,epsilon_bec,epsilon_bsc,F_bec,F_bsc)

% translate plots to capacity
C = 0.0001:0.0001:1;
C_bec = 1-epsilon_bec;
C_bsc = 1-Hb(epsilon_bsc);

[~,bec_indices] = min(abs(C_bec - C(:)),[],1);
[~,bsc_indices] = min(abs(C_bsc - C(:)),[],1);

% plot
currLegend = [];
figure; hold all
for i=1:length(rhoToPlot)
    for j=1:length(sToPlot)
        temp_bec = F_bec(i,j,bec_indices);
        temp_bsc = F_bsc(i,j,bsc_indices);
        plot(C,temp_bec(:),'--','Linewidth',2)
        plot(C,temp_bsc(:),'-','Linewidth',2)
        currLegend = [currLegend; strcat('BEC : S =',num2str(sToPlot(i)),' \rho =',num2str(rhoToPlot(j)));....
            strcat('BSC : S =',num2str(sToPlot(i)),' \rho =',num2str(rhoToPlot(j)))];
    end
end
grid on; grid minor;
xlabel('C [bits/Tx]'); ylabel('F_0')
legend(currLegend);
title('F_0 (s,\rho) : BEC and BSC')

end

function [h] = Hb(x)
h = zeros(size(x));
h(x == 0) = 0;
h(x>0) = -1*x(x>0).*log2(x(x>0)) - (1-x(x>0)).*log2(1 - x(x>0));
end


