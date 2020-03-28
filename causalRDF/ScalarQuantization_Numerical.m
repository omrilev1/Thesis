% Simulation of RDF in the low rate regime, using scalar quantizer
close all; clear all; clc;

Delta = 5 + (0.45:0.25:5.5);
sigma_n = 1/4;
sigma_x_sigma_n = sigma_n^2/(1 + sigma_n^2);

niters = 1e7;
iterLen = 1;


dx = 1e-3;

R_temp = 0.5*Delta.*exp(-Delta.^2/4);

D = zeros(size(Delta)); R = zeros(size(Delta));
for i=1:length(Delta)
    
    grid_quant = (-5:1:5)*Delta(i);
    err = 0;
    
    R(i) = calc_rate(Delta(i));
    for iter=1:niters
        
        x = randn(iterLen,1);
        y = x + sigma_n * randn(size(x));
        
        % scalar quantization
        [~,minIdx] = min(abs(x - grid_quant),[],2);
        minVal = grid_quant(minIdx);minVal = reshape(minVal,[],1);
        
        if abs(minVal) > 0
            disp('something');
        end
        % estimation
        x_y_mmse = y/(1 + sigma_n^2);
        
%         
        dt = linspace(0,1,1e4);
        x_axis = (minVal - 0.5*Delta(i) - x_y_mmse) + dt.*Delta(i);
        
        
%         x_axis2 = (minVal - 0.5*Delta(i) - x_y_mmse) : dx : (minVal + 0.5*Delta(i) - x_y_mmse);
        func = exp(-x_axis.^2 / (2*sigma_x_sigma_n));
        
        x_mmse = x_y_mmse + sum(x_axis.*func*dx,2)./sum(func*dx,2);
        
        err = err + sum((x - x_mmse).^2);
    end
    
    D(i) = err/(niters*iterLen);
    display(strcat('Finished Delta = ',num2str(Delta(i)),' D = ',num2str(D(i))));
    
end

D_temp = min(D) : 1e-5 : sigma_x_sigma_n;
figure;hold all;
plot(D,R,'-o','LineWidth',2);
plot(D_temp,0.5*log(sigma_x_sigma_n./D_temp),'--','LineWidth',2);
plot(D_temp,0.5*log(1./D_temp),'--','LineWidth',2);
grid on; grid minor;
xlabel('Distortion'); ylabel('Rate');
title('Low Rate Approximation');
legend('Scalar Quantization - empirical curve','Optimal Curve - Two-Sided SI','No SI');


function R = calc_rate(Delta)

grid_quant = (-25:1:25)*Delta;
R = 0;
dx = 1e-3;
x = -12:dx:12;
pdf = (1/sqrt(2*pi)) * exp(-x.^2 / 2);

for idx = 1 : length(grid_quant)
    
    valid_x = (abs(x - grid_quant(idx)) < Delta/2);
    
    if sum(valid_x) == 0
        continue
    end
    prob = sum(pdf(valid_x))*dx;
    
    if prob > 0
        R = R - prob*log(prob);
    end
    
end
end
