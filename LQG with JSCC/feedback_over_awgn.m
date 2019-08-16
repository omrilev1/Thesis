function [theta_final,MSE_final] = feedback_over_awgn(theta,P,Pin,snr,delta_snr,N_feedback,p_alias,type)
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

snr_F = 10^(snr/10);
snr_B = 10^((snr + delta_snr)/10);

% generate feedback constants

d = sqrt(3*P);
delta = d;
dither = 2*d*(rand(1,N_feedback) - 0.5); % uni[-delta,delta]

% calculate feedback estimation constants
[MSE,power_scale,beta,gamma] = calc_feedback_coeffs(P,Pin,sigma_n,sigma_n_feedback,N_feedback,p_alias,type);

%% Initialization : generate the first estimator

% first iteration - simple uncoded transmission
X(1) = sqrt(P/Pin)*theta;
Y = X(1) + sigma_n*randn;
theta_hat(1) = sqrt(Pin/P)*Y;

for i=1:N_feedback
    
    
    switch type
        case 'modulo'
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
        case 'linear'
            %% Terminal B : estimate the new estimation error and sends back to Terminal A
            % feedback signal generation
            X_tilde = gamma(i)*theta_hat(i);
            % feedback channel
            Y_tilde = X_tilde + sigma_n_feedback*randn;
            
            %% Terminal A : estimate the new Rx estimation error and sends to Terminal B
            epsilon_tilde = Y_tilde - gamma(i)*theta;
            X(i+1) = power_scale(i)*epsilon_tilde/gamma(i);
            Y = X(i+1) + sigma_n*randn;
            
            %% Terminal B : estimate the new estimation error
            theta_hat(i+1) = theta_hat(i) - beta(i)*Y;
            
        case 'linear optimized'
            
            if N_feedback == 1
                break;
            else                
                
                % we simplify the whole transmission scheme to matrix
                % operations, and then perform the matrix operation 1 time
                % (we can write it down by interactive feedback scheme)
                a = snr_F/snr_B;
                b = N_feedback/(N_feedback - 1);
                c = b*snr_F*(1 + a);
                
                % find the optimal power coefficent \gamma
                % for now we take the nominal value from the paper
%                 gamma = 1/sqrt(N_feedback);% (N_feedback - 1)/N_feedback;

                [gamma] = calc_opt_gamma(a,b,c,N_feedback);
                theta_norm = sqrt((1-gamma)*N_feedback*snr_F/Pin)*theta;

                % find the optimal \beta, we need to find the smallest
                % positive root, and check for some regularity conditions
                
                % optimal constant beta
                %                 beta = sqrt((N_feedback - 1)/(N_feedback + (1 + 10^(-1*(delta_snr/10))*N_feedback*gamma*10^(snr/10))));
                p = zeros(1,2*N_feedback+1);
                p(1) = 1;
                p(end-2) = -1*(N_feedback + (1 + snr_F/snr_B)*N_feedback*gamma*snr_F);
                p(end) = N_feedback - 1;
                beta = roots(p);
                beta = beta((beta > 0) & (beta < 1));
                beta = real(beta(imag(beta) < 1e-9));
                
                for k = 1:length(beta)
                    [curr_beta,curr_idx] = min(beta);
                    if curr_beta < 1e-2
                        beta(curr_idx) = [];
                    else
                        beta = curr_beta;
                        break
                    end
                end
                
                q_vec = sqrt((1-beta^2)/(1 - beta^(2*N_feedback))) * beta.^(0:1:(N_feedback-1));
                
                % construct the precoding matrix F
                F = zeros(N_feedback,N_feedback);
                for j=2:N_feedback
                    curr_row = -1*((1-beta^2)/(1 + snr_F/snr_B))*beta.^((1:(j-1)) - 2);
                    F(j,1:(j-1)) = fliplr(curr_row);
                end
                
                %% The actual estimation is via the next equations
                % y_vec = F(z_vec + n_vec) + g_vec*\theta + z_vec
                % \hat{\theta} = (g_vec)^T * y_vec
                
                z_vec = randn(N_feedback,1);
                n_vec = sqrt(snr_F/snr_B)*randn(N_feedback,1);
                theta_hat(N_feedback + 1) = q_vec * (F*(z_vec + n_vec) + q_vec(:) * theta_norm + z_vec);
                
                theta_hat(N_feedback + 1) = theta_hat(N_feedback + 1) * sqrt(Pin/((1-gamma)*N_feedback*snr_F)); 
            end
            break
    end
    
end

%% Final Stage : information signal is the signal in the last iteration
theta_final = theta_hat(N_feedback+1);
MSE_final = MSE(N_feedback);

end

function [gamma_opt] = calc_opt_gamma(a,b,c,N)

    p = 10; 
    
    gamma = 1/sqrt(N);
    for k=1:30
        gamma = 1 - (a*(b + c*gamma)^N + b + c)/(N*c);
    end
    
    if ((gamma > 1) || (gamma < 0) || abs(gamma - 1/sqrt(N)) > p/sqrt(N))
        gamma_opt = 1/sqrt(N);
    else
        gamma_opt = gamma;
    end
end

