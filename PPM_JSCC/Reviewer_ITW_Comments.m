clear clc close all

ENR_dB = 4:0.125:25; 
K = 0.058; 
K_reviewer = 0.073; 

ENR = 10.^(ENR_dB/10) ;
beta = (312 * pi^0.5)^(1/3) .* ENR.^(-5/6) .* exp(ENR./6) ;

% from the paper 
D_s = ((13/8)./(beta .* ENR).^2) .* (1+16/13 * sqrt(beta .* ENR) .* exp(-ENR));

% from the paper 
D_s_reviewer = ((13/8)./(beta .* ENR).^2) .* (1+16/13 * sqrt(ENR/2) .* exp(-ENR/4));

% what I think is correct 
P_L = beta .* sqrt(ENR) .* exp(-ENR/2)/16/sqrt(pi);

% from the paper 
D_L = (1/6)*(1 + 2./beta + 4./beta.^2); 
D = D_s + P_L .* D_L; 

D_reviewer = D_s_reviewer + P_L .* D_L; 

% from the paper 

D_lower_bound = K * exp(-ENR./3) .* ENR.^(-1/3) ; 
D_lower_bound_reviewer = K_reviewer * exp(-ENR./3) .* ENR.^(-1/3) ;

% from the paper 
SDR = 1/12 .* (1./D) ;
SDR_dB = 10 * log10(SDR) ;
SDR_lower_bound = 1/12 .* (1./D_lower_bound) ;
SDR_dB_lower_bound = 10 * log10(SDR_lower_bound) ; 

SDR_reviewer = 1/12 .* (1./D_reviewer) ; 
SDR_dB_reviewer = 10 * log10(SDR_reviewer) ;
SDR_lower_bound_reviewer = 1/12 .* (1./D_lower_bound_reviewer) ;
SDR_dB_lower_bound_reviewer = 10 * log10(SDR_lower_bound_reviewer);

figure; 
plot(ENR_dB, SDR_dB, '-ok', ...
    ENR_dB, SDR_dB_reviewer, '-.+c',...
    ENR_dB, SDR_dB_lower_bound,'--sm', ...
    ENR_dB, SDR_dB_lower_bound_reviewer,':*b');
grid on; box on; ylim([0, 80]); xlabel('ENR [dB]'); ylabel('SDR [dB]');
legend('SDR from (11)', 'SDR from (11), corrected by reviewer', 'Lower Bound from (14)', 'Lower Bound from (14), corrected by reviewer')


% Uniform Optimization 
alpha = 1:1e-3:16; 
[minVal,minIdx] = min(13/8 ./ (alpha.^2) + alpha/(96*sqrt(pi)));
disp(strcat('Optimal K is ',num2str(minVal)));

% Gaussian Optimization 
alpha = 1:1e-3:16; 
[minVal_Gauss,minIdx_Gauss] = min(13/8 ./ (alpha.^2) + 2*alpha);
disp(strcat('Optimal K Gauss is ',num2str(minVal_Gauss)));
