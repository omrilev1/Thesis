% compaer SDR Of PPM and Linear
clear all; close all; clc;
maxEN0 = 20;
optEN0 = 1;% [1 5 9 13 17];

for i=1:length(optEN0)
    EN0 = 0:0.01:maxEN0;
    EN0_lin = 10.^(EN0/10);
    beta = (13/8)^(1/3) * (EN0_lin.^(-5/6)) .* exp(EN0_lin/6);
    D_S = (13/8 + 2*(sqrt(2*beta.*EN0_lin) - 1).*exp(-EN0_lin.*(1 - 1./sqrt(2*beta.*EN0_lin)).^2))./((sqrt(beta.*EN0_lin) - 1/sqrt(2)).^4)...
        + exp(-beta.*EN0_lin)./(beta.^2);
    D_L = 2*beta.*sqrt(EN0_lin).*exp(-1*EN0_lin/2).*(1 + 3.*sqrt(2*pi./EN0_lin) + 12*exp(-1)./(beta.*sqrt(EN0_lin)) + ...
        8*exp(-1)./(sqrt(8*pi)*beta) + sqrt(8./(pi*EN0_lin)) + 12^(3/2) * exp(-3/2) ./(beta.*sqrt(32*pi*EN0_lin))) ...
        + sqrt(8*pi)*beta.*exp(-EN0_lin) .* (1 + 4*exp(-1) ./ (beta*sqrt(2*pi)));
    SDR_PPM = 1./(D_S + D_L);
    % D_S = 13/8 ./ (beta.*EN0_lin).^2;
    % D_L = 2*beta.*sqrt(EN0_lin).*exp(-1*EN0_lin/2);
    
    optEN0_Lin = 10^(optEN0(i)/10);
    beta = (13/8)^(1/3) * (optEN0_Lin.^(-5/6)) .* exp(optEN0_Lin/6);
    D_S_opt = (13/8 + 2*(sqrt(2*beta.*EN0_lin) - 1).*exp(-EN0_lin.*(1 - 1./sqrt(2*beta.*EN0_lin)).^2))./((sqrt(beta.*EN0_lin) - 1/sqrt(2)).^4)...
        + exp(-beta.*EN0_lin)./(beta.^2);
    D_L_opt = 2*beta.*sqrt(EN0_lin).*exp(-1*EN0_lin/2).*(1 + 3.*sqrt(2*pi./EN0_lin) + 12*exp(-1)./(beta.*sqrt(EN0_lin)) + ...
        8*exp(-1)./(sqrt(8*pi)*beta) + sqrt(8./(pi*EN0_lin)) + 12^(3/2) * exp(-3/2) ./(beta.*sqrt(32*pi*EN0_lin))) ...
        + sqrt(8*pi)*beta.*exp(-EN0_lin) .* (1 + 4*exp(-1) ./ (beta*sqrt(2*pi)));
    SDR_PPM_Opt = 1./(D_S_opt + D_L_opt);
    SDR_PPM_Opt_Approx = 0.283 * exp(EN0_lin/3) .* (EN0_lin.^(1/3));
    
    figure;semilogy(EN0,max(1,1./(D_S + D_L)),'LineWidth',2); hold on;
    semilogy(EN0,1 + 2*EN0_lin,'-.','LineWidth',2); hold on;
    semilogy(EN0,max(1,SDR_PPM_Opt),'-','LineWidth',2); hold on;
    semilogy(EN0,max(1,SDR_PPM_Opt_Approx),'.-','LineWidth',2); hold on;
    semilogy(EN0,(1 + EN0_lin .^ 2),'--','LineWidth',2)
    grid on; grid minor;
    xlabel('ENR [dB]'); ylabel('SDR');
    lgd = legend({'PPM','Linear','PPM Opt','PPM Opt Approx ','Squared Profile'},'FontSize',12);
    %     legend(lgd);
    title(strcat('Opt EN0 = ' ,num2str(optEN0(i)),' [dB]'));
end


%% M Linear Layers and 1 PPM
alpha = 0.005:0.005:10;
x = 0.005:0.005:10;
M = 2:1:30;
[xGrid,alphaGrid,Mgrid] = meshgrid(x,alpha,M);

cost_linear = (exp(alphaGrid)./xGrid) + ...
    0.5*(xGrid).*((exp(-2*alphaGrid).*(1 - exp(-alphaGrid.*Mgrid))./(1 - exp(-1*alphaGrid)))).*...
    (exp(alphaGrid*2) - 1).*...
    (1 + sqrt((1 + 4*exp(alphaGrid.*(2 + 1)))./(1 - exp(alphaGrid * 2)).^2));
E_PPM = 10^(8.8/10) * xGrid.*exp(-1*alphaGrid.*Mgrid);

[v,loc] = min(0.5*cost_linear(:) + E_PPM(:));
[ii,jj,k] = ind2sub(size(E_PPM),loc);
x_opt = x(jj); 
alpha_opt = alpha(ii); 
M_opt = 15; 
opt_cost_verify = 0.5 * (exp(alpha_opt)./x_opt) + ...
    0.25*(x_opt).*((exp(-2*alpha_opt).*(1 - exp(-alpha_opt.*(M_opt - 1)))./(1 - exp(-1*alpha_opt)))).*...
    (exp(alpha_opt*2) - 1).*...
    (1 + sqrt((1 + 4*exp(alpha_opt.*(2 + 1)))./(1 - exp(alpha_opt * 2)).^2)) + 10^(8.8/10) * x_opt.*exp(-1*alpha_opt.*M_opt);

[minEnergy_linear] = min(min(min(0.5*cost_linear + E_PPM)));


%% Second Option
% optimal PPM Energy for squared profile:
% First layer is linear and its energy is x1/2
% The second layer is with PPM - we optimize it for N0, and demand it to be greater than 1 on N_1

N1 = 0.005:0.005:300;

bestEnergy = 1e3;
useApprox = false;
for i=1:length(N1)
    N0 = N1(i):0.05:200;
    
    %% first layer energy
    E0 = 0.5/N1(i);
    
    %% PPM Energy
    
    if useApprox
        % approximation
        kappa = 0.283;
        E1 = N0 .* Lambert_W(((N1(i)./N0).^2 / kappa).^3);
    else
        % exact
        E1 = N0 .* inv_ppm_sdr(N0./N1(i),1);
    end
    
    % final optimization
    E_tot = E0(:) + E1(:);
    if min(E_tot) < bestEnergy
        bestEnergy = min(E_tot);
    end
    
    if mod(i,200) == 0
        disp(strcat('Minimal Energy is = ',num2str(bestEnergy)));
    end
end

disp(strcat('Minimal Energy is = ',num2str(bestEnergy)));


function EN0 = inv_ppm_sdr(ratio,val)

EN_Opt = -20:0.05:13;
EN_Opt_lin = 10.^(EN_Opt/10);

[EN_Opt_Mat,ratio_mat] = meshgrid(EN_Opt_lin(:),ratio(:));
EN_Eval = EN_Opt_Mat .* ratio_mat;

beta = (13/8)^(1/3) * (EN_Opt_Mat.^(-5/6)) .* exp(EN_Opt_Mat/6);

if min(min(beta .* EN_Eval)) < 1/2
    disp('The Approximation do not holds anymore!');
end
D_S_opt = (13/8 + 2*(sqrt(2*beta.*EN_Eval) - 1).*exp(-EN_Eval.*(1 - 1./sqrt(2*beta.*EN_Eval)).^2))./((sqrt(beta.*EN_Eval) - 1/sqrt(2)).^4)...
    + exp(-beta.*EN_Eval)./(beta.^2);
D_L_opt = 2*beta.*sqrt(EN_Eval).*exp(-1*EN_Eval/2).*(1 + 3.*sqrt(2*pi./EN_Eval) + 12*exp(-1)./(beta.*sqrt(EN_Eval)) + ...
    8*exp(-1)./(sqrt(8*pi)*beta) + sqrt(8./(pi*EN_Eval)) + 12^(3/2) * exp(-3/2) ./(beta.*sqrt(32*pi*EN_Eval))) ...
    + sqrt(8*pi)*beta.*exp(-EN_Eval) .* (1 + 4*exp(-1) ./ (beta*sqrt(2*pi)));
SDR_PPM_Opt = 1./(D_S_opt + D_L_opt);
% SDR_PPM_Opt = max(1,SDR_PPM_Opt);

[minVal,valIdx] = min(abs(SDR_PPM_Opt - val),[],2);
EN0 = EN_Opt_lin(valIdx);

end