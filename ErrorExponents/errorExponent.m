clear ; close all; clc;

s = [0.5];
rho = [0.5];

% [mesh_rho,mesh_s] = meshgrid(rho,s);

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

% translate plots to capacity
C = 0.0001:0.0001:1;
C_bec = 1-epsilon_bec; 
C_bsc = 1-Hb(epsilon_bsc);

[~,bec_indices] = min(abs(C_bec - C(:)),[],1);
[~,bsc_indices] = min(abs(C_bsc - C(:)),[],1);

%2D plot
figure; hold all
for i=1:length(rho)
    for j=1:length(s)
        temp_bec = F_bec(i,j,bec_indices);
        temp_bsc = F_bsc(i,j,bsc_indices);
        plot(C,temp_bec(:),'--','Linewidth',2)
        plot(C,temp_bsc(:),'-','Linewidth',2)
    end
end
grid on; grid minor;
xlabel('C [bits/Tx]'); ylabel('F_0')
legend('BEC : s = 0.1 , \rho = 0.1','BSC : s = 0.1 , \rho = 0.1',...
    'BEC : s = 0.5 , \rho = 0.1','BSC : s = 0.5 , \rho = 0.1',...
    'BEC : s = 0.1 , \rho = 0.5','BSC : s = 0.1 , \rho = 0.5',...
    'BEC : s = 0.5 , \rho = 0.5','BSC : s = 0.5 , \rho = 0.5');
title('F_0 (s,\rho) : BEC and BSC')
function [h] = Hb(x)
    h = zeros(size(x));
    h(x == 0) = 0;
    h(x>0) = -1*x(x>0).*log2(x(x>0)) - (1-x(x>0)).*log2(1 - x(x>0));
end