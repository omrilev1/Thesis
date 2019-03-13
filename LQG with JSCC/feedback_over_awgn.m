function [theta_final] = feedback_over_awgn(theta,P,m,snr,delta_snr,N_feedback)
% This function simulate the communication with noisy feedback scheme of
% Assaf and Ofer []. 
% Inputs : 
%   o theta : information signal. Inside the power limits (needs to be
%   normalized before and multiplied after by appropriate scale factor) 
%   o P : input power constraint
%   o m : PAM order
%   o snr : channel snr
%   o delta_snr : snr + delta_snr is the snr of the feedback path
X = zeros(1,N_feedback);
theta_hat = zeros(1,N_feedback);

% Initialize noise parameters
sigma_n = sqrt(P/10^(snr/10));
sigma_n_feedback = sqrt(P/10^((snr + delta_snr)/10));

% generate feedback constants 

d = sqrt(12*P);
delta = d/2;
dither = d*(rand(1,N_feedback) - 0.5); % uni[-d/2,d/2]

% calculate feedback estimation constants
p_m = 1e-6/(2*N_feedback);
lambda = 3/(qfuncinv(p_m/2))^2;
gamma = sqrt((lambda*d^2/12 - sigma_n_feedback^2)/


%% Initialization : generate the PAM symbol
% PAM constellation scaling factor
R = log2(m);
eta = sqrt(P)*sqrt(3/(2^(2*R) - 1));

% generate PAM points 
pam_symbol = pammod(0:1:(m-1),m,0,'gray');

% quantize theta and generate the transmitted symbol
[~,quant_idx] = min(abs(theta - eta*pam_symbol));

% first iteration - simple uncoded transmission
X(1) = eta*pam_symbol(quant_idx);
Y = X(1) + sigma_n*randn;
theta_hat(1) = Y;

for i=1:N_feedback
    
    %% Terminal B : estimate the new estimation error and sends back to Terminal A 
    % feedback signal generation 
    X_tilde = mod(gamma(i)*theta_hat(i) + dither(i) + delta,2*delta) - delta;
    % feedback channel 
    Y_tilde = X_tilde + sigma_n_feedback*randn;
    
    %% Terminal A : estimate the new Rx estimation error and sends to Terminal B 
    epsilon_tilde = mod(Y_tilde - gamma(i)*X(1) - dither(i) + delta,2*delta) - delta;
    X(i+1) = alpha(i)*epsilon_tilde;
    Y = X(i+1) + sigma_n*randn;
    
    %% Terminal B : estimate the new estimation error
    theta_hat(i+1) = theta_hat(i) - beta(i+1)*Y;
    
end

%% Final Stage : information signal decoding 
decod_idx = min(abs(theta_hat(N_feedback+1) - eta*pam_symbol));
theta_final = pam_symbol(decod_idx);

end


