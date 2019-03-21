function [MSE,power_scale,beta,gamma] = calc_feedback_coeffs(P,Pin,sigma_n,sigma_n_feedback,Nfeedback,p_alias)
%This function calculate the feedback shceme coefficents, i.e the
%estimaiton coefficients, the zooming factor, the input power scale, and
%the MSE

MSE = zeros(1,Nfeedback);
power_scale = zeros(1,Nfeedback);
beta = zeros(1,Nfeedback);
gamma = zeros(1,Nfeedback);


snr = P/sigma_n;
% initialization

for n = 1 : Nfeedback
    if n==1
        MSE(1) = Pin/snr;
    else
        MSE(n) = (1 - power_scale(n-1)^2 * MSE(n-1) / (P + sigma_n^2))^2 * MSE(n-1) + ...
            ((power_scale(n-1) * MSE(n-1) / (P + sigma_n^2))^2) * ...
            ((power_scale(n-1)*sigma_n_feedback/gamma(n-1))^2 + sigma_n^2);
    end
    gamma(n) = sqrt((3*P/(qfuncinv(p_alias/2))^2) - sigma_n_feedback^2)/sqrt(MSE(n));
    power_scale(n) = sqrt(P/(MSE(n) + (sigma_n_feedback/gamma(n))^2));
    beta(n) = (power_scale(n)*MSE(n))/(P + sigma_n^2);
    
    
end
end

