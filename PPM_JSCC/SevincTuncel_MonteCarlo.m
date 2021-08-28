% Digital PPM: Simulation of Burnashev scheme
% The simulation is for uniform source
clear all; close all; clc;
addpath(genpath('fast_linear_interpolation'));

ENR = (4:1:16);
ENRlin = 10.^(ENR/10);

Nrun = 6e6; %

delta = 1/2 + 1/(2*sqrt(2));
dist = {'Uniform'};

final_MSE = zeros(length(dist),length(ENR));
final_Beta = zeros(length(dist),length(ENR));

for distIdx = 1:length(dist)
    disp('===============================');
    disp(strcat('Current Distribution: ',dist{distIdx}));
    for i=1:length(ENR)
        % pick beta to find the optimal point of the curve
        
        switch dist{distIdx}
            case 'Uniform'
                beta_opt = 3 * exp(ENRlin(i)/6);
                beta = beta_opt; % max(1,beta_opt*(0.75:(1/4):1.55));
            case 'Gaussian'
                beta_opt = 2*exp(ENRlin(i)/6);% 1.6*exp(ENRlin(i)/12);
                beta = beta_opt; % max(1,beta_opt*(0.75:(1/4):1.55));
        end
        MSE = zeros(size(beta));
        currENRlin = ENRlin(i);
        
        %calculate optimal Quantization for the Gaussian case
        if strcmp(dist(distIdx),'Gaussian')
            % equations from Sevinc' Tuncel paper, appendix A, for the calculation
            % of the compander
            dx = 0.5e-5; x = -6.75:dx:6.75;
            c = 2.2219;delta = 1/2 + 1/(2*sqrt(2));
            lambda_x = (1/(6^(1/3) * c^(2/3))) * exp(-x.^2 / 6) ./ ((2*c*x.^2 + 0.6573).^(1/3));
            G = cumsum(lambda_x*dx);
            
        else
            c = 0.8281;dx = 0.5e-5; x = -1/2:dx:1/2;
            lambda_x = (1/(6^(1/3) * c^(2/3))) * 1 ./ ((2*c*x.^2 + 0.1385).^(1/3));
            G = cumsum(lambda_x*dx);
        end
        
        for beta_idx = 1:length(beta)
            curr_beta = beta(beta_idx);
            currMSE = calc_perf_Tuncel(curr_beta,Nrun,currENRlin,dist{distIdx},G(:),x(:),dx);
            MSE(beta_idx) = currMSE;
        end
        %     [optVal,optIdx] = min(abs(MSE - 13/8 ./ (beta.*currENRlin).^2));
        [optVal,optIdx] = min(MSE(1:beta_idx));
        final_MSE(distIdx,i) = MSE(optIdx);
        final_Beta(distIdx,i) = beta(optIdx);
        disp(strcat('Finished ENR = ',num2str(ENR(i))));
    end
    
    switch dist{distIdx}
        case 'Uniform'
            optSDR_Analytic_Exact = ((delta^(2/3))/4) .* exp(-ENRlin/3);
            opt_beta = 2 * exp(ENRlin/6);
            save('UniformTuncelOptSDR.mat','ENR','final_MSE','optSDR_Analytic_Exact','final_Beta','opt_beta');
            
            figure;subplot(211);
            hold all;
            plot(ENR,10*log10(1/12./final_MSE(distIdx,:)),'-.ko','LineWidth',1.5);
            plot(ENR,10*log10(4/12./optSDR_Analytic_Exact),'--ms','LineWidth',1.5);
            xlabel('ENR [dB]'); ylabel('SDR [dB]'); legend({'Empirical','Approximation'},'Location','northwest','FontSize',12);
            grid on; grid minor;
            subplot(212);
            semilogy(ENR,final_Beta(distIdx,:),'-.ko','LineWidth',1.5);hold on;
            semilogy(ENR,opt_beta,'--ms','LineWidth',1.5);
            xlabel('ENR [dB]'); ylabel('Optimal Beta');
            grid on; grid minor;
            
        case 'Gaussian'
            optSDR_Analytic_Exact = exp(-ENRlin/3);
            opt_beta = 2 * exp(ENRlin/6);
            save('GaussianTuncelOptSDR.mat','ENR','final_MSE','optSDR_Analytic_Exact','final_Beta','opt_beta');
            
            figure;subplot(211);
            hold all;
            plot(ENR,10*log10(1./final_MSE(distIdx,:)),'-.ko','LineWidth',1.5);
            plot(ENR,10*log10(1./optSDR_Analytic_Exact),'--ms','LineWidth',1.5);
            xlabel('ENR [dB]'); ylabel('SDR [dB]'); legend({'Empirical','Approximation'},'Location','northwest','FontSize',12);
            grid on; grid minor;
            subplot(212);
            semilogy(ENR,final_Beta(distIdx,:),'-.ko','LineWidth',1.5);hold on;
            semilogy(ENR,opt_beta,'--ms','LineWidth',1.5);
            xlabel('ENR [dB]'); ylabel('Optimal Beta');
            grid on; grid minor;
            
    end
end


function MSE = calc_perf_Tuncel(beta,Nrun,ENRlin,distType,G,x,dx)

beta_round = round(beta);
orthMat = sparse(eye(beta_round));

currMSE = zeros(1,Nrun);
currENR = ENRlin;

if (beta_round == 1)
    switch distType
        case 'Uniform'
            MSE = 1/12;
        case 'Gaussian'
            MSE = 1;
    end
    return
end
% Generate Quantization Table

parfor n=1:Nrun
    
    % generate source and quantize
    switch distType
        case 'Uniform'
            S = rand - 0.5;  % [-0.5,0.5]
            
        case 'Gaussian'
            S = randn;
            
    end
    
    % companding: we take the ivalue from G, and then
    % linearly interpolate with neighbours to get finer resolution
    
%     if round((S - min(x))/dx) == 0
%         S_Companding = G(1);
%     elseif round((S - min(x))/dx) >= (length(G) - 1)
%         S_Companding = G(end);
%     else
%         S_0  = G(round((S - min(x))/dx) + 1);
%         S_1  = G(round((S - min(x))/dx) + 2);
%         S_m1 = G(round((S - min(x))/dx));
%         
%         S_Companding = interp1q([-1;0;1],[S_m1,;S_0;S_1],(S - min(x))/dx - round((S - min(x))/dx));
%     end
    S_Companding = G(round((S - min(x))/dx) + 1);
    
    % Scalar uniform quantization
    S_q = round((beta_round - 1)* S_Companding);
    
    % Modulation
    TxPulse = sparse(orthMat(S_q + 1,:));
    TxPulse = sqrt(currENR)*TxPulse;
    
    % AWGN
    noise = randn(size(TxPulse));
    noise = sqrt(1/2) * noise;% sqrt(2*Fs/W)*noise;
    r = TxPulse + noise;
    
    % Correlator receiver: multiplication with eye matrix 
    % is equivalent to pick the maximal element of the vector 
%     corr = orthMat*r(:);
    
    [~,maxIdx] = max(r(:));
    
    % De-quantization + De-companding
    iHat = (maxIdx - 1 + 0.5)/beta_round;
    [~,idx] = min(abs(G - iHat));
    sHat  = x(idx);
    
    currMSE(n) = (S - sHat)^2;
end

MSE = sum(currMSE)/Nrun;
end