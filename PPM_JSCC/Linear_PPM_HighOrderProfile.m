close all; clear all; clc;

L =3:2:100;
alpha = 0.01:0.01:20;
x = 0.01:0.01:20;
M = 2:1:35;
[xGrid,alphaGrid,Mgrid] = meshgrid(x,alpha,M);

kappa = 0.283;
minEnergy_linear = zeros(1,length(L));
minEnergy_ppm = zeros(1,length(L));

parfor l=1:length(L)
    
    %% PPM
    cost_ppm = 0.5*(exp(alphaGrid)./xGrid).^(L(l) -1) + ...
        0.25*(xGrid).*((exp(-2*alphaGrid).*(1 - exp(-alphaGrid.*Mgrid))./(1 - exp(-1*alphaGrid)))).*...
        (exp(alphaGrid*(L(l))) - 1).*...
        (1 + sqrt((1 + 4*exp(alphaGrid.*(L(l) + 1)))./(1 - exp(alphaGrid * L(l))).^2)) + ...
        (xGrid).* Lambert_W(((exp(alphaGrid * (L(l) - 2)))/(kappa)).^3) .*exp(-alphaGrid.*Mgrid)./ (exp(alphaGrid) - 1);% (exp(alphaGrid) - 1);
    minEnergy_ppm(l) = min(min(min(cost_ppm)));
    
    %% Linear 
    cost_linear = (exp(alphaGrid)./xGrid).^(L(l) -1) + ...
        0.5*(xGrid).*((exp(-2*alphaGrid)./(1 - exp(-1*alphaGrid)))).*...
        (exp(alphaGrid*(L(l))) - 1).*...
        (1 + sqrt((1 + 4*exp(alphaGrid.*(L(l) + 1)))./(1 - exp(alphaGrid * L(l))).^2));
    
    minEnergy_linear(l) = 0.5 * min(min(min(cost_linear)));
    
    disp(strcat('Finished L = ',num2str(L(l))));
end

figure;semilogy(L,minEnergy_linear,'-o','LineWidth',2);
hold on; semilogy(L,minEnergy_ppm,'-*','LineWidth',2);
xlabel('Profile Order - L'); ylabel('Minimum Energy');
legend('Linear Layers','PPM Layers');
grid on; grid minor;
