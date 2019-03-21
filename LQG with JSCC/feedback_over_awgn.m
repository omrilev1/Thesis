function [theta_final,MSE_final] = feedback_over_awgn(theta,P,Pin,snr,delta_snr,N_feedback,p_alias)
% This function simulate the communication with noisy feedback scheme of
% Assaf and Ofer []. 
% Inputs : 
%   o theta : information signal. Inside the power limits (needs to be
%   normalized before and multiplied after by appropriate scale factor) 
%   o P : input power constraint
%   o snr : channel snr
%   o delta_snr : snr + delta_snr is the snr of the feedback path
X = zeros(1,N_feedback);
theta_hat = zeros(1,N_feedback);

% Initialize noise parameters
sigma_n = sqrt(P/10^(snr/10));
sigma_n_feedback = sqrt(P/10^((snr + delta_snr)/10));

% generate feedback constants 

d = sqrt(3*P);
delta = d;
dither = 2*d*(rand(1,N_feedback) - 0.5); % uni[-delta,delta]

% calculate feedback estimation constants
[MSE,power_scale,beta,gamma] = calc_feedback_coeffs(P,Pin,sigma_n,sigma_n_feedback,N_feedback,p_alias);

%% Initialization : generate the dirst estimator

% first iteration - simple uncoded transmission
X(1) = sqrt(P/Pin)*theta;
Y = X(1) + sigma_n*randn;
theta_hat(1) = sqrt(Pin/P)*Y;

for i=1:N_feedback
    
    %% Terminal B : estimate the new estimation error and sends back to Terminal A 
    % feedback signal generation 
    X_tilde = mod(gamma(i)*theta_hat(i) + dither(i) + delta,2*delta) - delta;
    % feedback channel 
    Y_tilde = X_tilde + sigma_n_feedback*randn;
    
    %% Terminal A : estimate the new Rx estimation error and sends to Terminal B 
    epsilon_tilde = mod(Y_tilde - gamma(i)*theta - dither(i) + delta,2*delta) - delta;
    X(i+1) = power_scale(i)*epsilon_tilde/gamma(i);
    Y = X(i+1) + sigma_n*randn;
    
    %% Terminal B : estimate the new estimation error
    theta_hat(i+1) = theta_hat(i) - beta(i)*Y;
    
end

%% Final Stage : information signal is the signal in the last iteration 
theta_final = theta_hat(N_feedback+1);
MSE_final = MSE(N_feedback);

end


