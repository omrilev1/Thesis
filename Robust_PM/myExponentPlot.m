close all; clear all; clc;
p = 0:0.01:0.3;
epsilon = [0 0.05 0.1];
lambda = 0 : 0.01 : 1;
N = 1:50;
exponent = zeros(length(epsilon),length(p));
beta = 0 : 0.007 : 1;

% My Exponent plot
for n=1:length(N)
    for i=1:length(epsilon)
        for j=1:length(p)
            
            conv1 = binaryConvolutions(p(j),epsilon(i));
            conv2 = binaryConvolutions(p(j),1 - epsilon(i));
            
            conv = max(conv1,conv2);% max(conv1.^lambda,conv2.^lambda);
            
            exponent(i,j,n) = -log2(conv); % -1 - log2(conv)
            if exponent(i,j,n)  < 1/n
                exponent(i,j,n) = 0;
            end
        end
    end
end

% Anusha's Exponent
opt_beta = zeros(length(p),length(N));
for n=1:length(N)
    for j=1:length(p)
        opt_Psi = zeros(size(beta));
        
        for k=1:length(beta)
            
            curr_beta = beta(k);
            curr_p = p(j);
            
            curr_psi = @(lambda) lambda*curr_beta + log2(((2*curr_p)^(lambda) + (2*(1-curr_p))^(lambda))/2);
            lambda_opt = fminbnd(curr_psi,0,1);
            opt_Psi(k) = -log2(((2*curr_p).^(lambda_opt) + (2*(1-curr_p)).^(lambda_opt))/2) - lambda_opt*curr_beta;
        end
        
        [~,beta_opt_idx] = min(abs(beta - opt_Psi + 1/n));
        curr_beta_opt = beta(beta_opt_idx);
        
        opt_beta(j,n) = curr_beta_opt;
    end
    display(strcat('Finished n = ',num2str(n)));
end
% Plot : vs n
p_for_plot = [0.05,0.1,0.15,0.2,0.25];
for i=1:length(p_for_plot)
    figure; hold all
    
    [~,p_idx] = min(abs(p - p_for_plot(i)));
    plot(1./N,opt_beta(p_idx,:),'LineWidth',2);
    for j=1:length(epsilon)
        plot(1./N,reshape(exponent(j,p_idx,:),1,[]),'LineWidth',2);
    end
    grid on; grid minor;
    xlabel('1/n'); ylabel('Exponent');
    legend('Anusha','\epsilon = 0','\epsilon = 0.05','\epsilon = 0.1');
    title(strcat('Bound for p = ',num2str(p_for_plot(i))));
end

% Plot : vs p
figure; hold all

plot(p,opt_beta(:,length(N)),'LineWidth',2);
for j=1:length(epsilon)
    plot(p,1 - reshape(exponent(j,:,length(N)),1,[]),'LineWidth',2);
end
grid on; grid minor;
xlabel('p'); ylabel('Exponent');
legend('Anusha','\epsilon = 0','\epsilon = 0.05','\epsilon = 0.1');
title(strcat('Bound for n_{max} = ',num2str(length(N))));



function [conv] = binaryConvolutions(p,eps)
% This function calculate the convolution between 2 numbers
conv = p*(1-eps) + (1-p)*eps;
end

function [R] = binaryEntropyCalc(p)
% This function calculate the convolution between 2 numbers
H = zeros(size(p));
H(p > 0) = -p(p > 0).*log2(p(p > 0)) - (1-p(p > 0)).*log2(1-p(p > 0));
R = 1 - H;
end