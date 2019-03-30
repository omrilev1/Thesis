%% Simulation of LQG over power limited awgn channel
% we simulate 3 cases :
% - controller and observer has full access to the control signal u_t
% - only controller has the control signal
% - noisy feedback between both

% The sensor has access to the plant, and needs to transmit the current
% state to the controller who controls the plant. The sensor has acces to
% x_t and limited access to u_t, while the controller has access to u_t but no
% access to x_t

% The evolution of x_t is given by :
% x_(t+1) = alpha*x_t + w_(t+1) + u_t
% where w_(t+1) is wgn with variance W

% This script simulate only linear schemes , i.e no modulo operations

close all; clear all; clc;
debug = 0;
%% coefficeints initialization
% model coefficients
alpha = 2;
W = 1; % input gaussian variance
T = 2^7; % horizon
delta = 3;
P = delta^2/3;

% cost parameters
Q = 1;
R = 1;
ratio_up = 3; % upper limit of control coeffs ratio
ratio_down = 1/3; % lower limit of control coeffs ratio


SNR = [20 15 10 8 6];
deltaSNR = 15;
snrLin = 10.^(SNR/10);
snrFeedback_Lin = 10.^((SNR + deltaSNR)/10);

N_avg = 100;

initArrays;

P_alias = 1e-2;
N_feedback = 16;
% set timesharing parameters - we need to scale the powers to the overall
% avergae power will be P
uncoded_timesharing_factor = 150;

cnt = 0;
for i=1:length(SNR)
    sigma_z = sqrt(P/snrLin(i));
    sigma_v = sqrt(P/snrFeedback_Lin(i));
    
    % calculate LQG parameters for current run
    k_vec = calcLQG(Q,R,T,alpha);
    k_vec_up = calcLQG(Q,Q*ratio_up,T,alpha);
    k_vec_down = calcLQG(Q,Q*ratio_down,T,alpha);
    
    P_X_linear = zeros(T,N_avg);
    P_X_hat_linear = zeros(T,N_avg);

    P_X_partial = zeros(T,N_avg);
    P_X_hat_partial = zeros(T,N_avg);
    
    P_X_modulo = zeros(T,N_avg);
    P_X_hat_modulo = zeros(T,N_avg);
    
    P_X_feedback = zeros(T,N_avg);
    P_X_hat_feedback = zeros(T,N_avg);
    
    for n = 1 : N_avg
        for t=1:T
            
            %% plant
            if t==1
                % initialization
                curr_w = sqrt(W)*randn;
                x_fullAccess(i,t,n) = curr_w;
                x_zeroAccess(i,t,n) = curr_w;                
                x_partialAccess(i,t,n) = curr_w;
                x_zeroAccess_modulo(i,t,n) = curr_w;
                x_partialAccess_feedback(i,t,n) = curr_w;
                x_partialAccess_feedback_linear(i,t,n) = curr_w;
                
            elseif t==T
                curr_w = sqrt(W)*randn;
                % no control in the last stage - evolve state and break
                x_fullAccess(i,t,n) = alpha*x_fullAccess(i,t-1,n) + curr_w + u_fullAccess(i,t-1,n);
                J_fullAccess(i,t,n) = (J_fullAccess(i,t-1,n)*(T-1) + (Q*(x_fullAccess(i,t,n))^2))/T;
                
                % no control in the last stage - evolve state and break
                x_zeroAccess(i,t,n) = alpha*x_zeroAccess(i,t-1,n) + curr_w + u_zeroAccess(i,t-1,n);
                J_zeroAccess(i,t,n) = (J_zeroAccess(i,t-1,n)*(T-1) + (Q*(x_zeroAccess(i,t,n))^2))/T;   
                
                % no control in the last stage - evolve state and break
                x_partialAccess(i,t,n) = alpha*x_partialAccess(i,t-1,n) + curr_w + u_partialAccess(i,t-1,n);
                J_partialAccess(i,t,n) = (J_partialAccess(i,t-1,n)*(T-1) + (Q*(x_partialAccess(i,t,n))^2))/T;
                
                % no control in the last stage - evolve state and break
                x_zeroAccess_modulo(i,t,n) = alpha*x_zeroAccess_modulo(i,t-1,n) + curr_w + u_zeroAccess_modulo(i,t-1,n);
                J_zeroAccess_modulo(i,t,n) = (J_zeroAccess_modulo(i,t-1,n)*(T-1) + (Q*(x_zeroAccess_modulo(i,t,n))^2))/T;
                
                % no control in the last stage - evolve state and break
                x_partialAccess_feedback(i,t,n) = alpha*x_partialAccess_feedback(i,t-1,n) + curr_w + u_partialAccess_feedback(i,t-1,n);
                J_partialAccess_feedback(i,t,n) = (J_partialAccess_feedback(i,t-1,n)*(T-1) + (Q*(x_partialAccess_feedback(i,t,n))^2))/T;

                % no control in the last stage - evolve state and break
                x_partialAccess_feedback_linear(i,t,n) = alpha*x_partialAccess_feedback_linear(i,t-1,n) + curr_w + u_partialAccess_feedback_linear(i,t-1,n);
                J_partialAccess_feedback_linear(i,t,n) = (J_partialAccess_feedback_linear(i,t-1,n)*(T-1) + (Q*(x_partialAccess_feedback_linear(i,t,n))^2))/T;
                break
                
                
            else
                % evolution
                curr_w = sqrt(W)*randn;
                x_zeroAccess(i,t,n) = alpha*x_zeroAccess(i,t-1,n) + curr_w + u_zeroAccess(i,t-1,n);                
                x_fullAccess(i,t,n) = alpha*x_fullAccess(i,t-1,n) + curr_w + u_fullAccess(i,t-1,n);
                x_partialAccess(i,t,n) = alpha*x_partialAccess(i,t-1,n) + curr_w + u_partialAccess(i,t-1,n);
                x_zeroAccess_modulo(i,t,n) = alpha*x_zeroAccess_modulo(i,t-1,n) + curr_w + u_zeroAccess_modulo(i,t-1,n);
                x_partialAccess_feedback(i,t,n) = alpha*x_partialAccess_feedback(i,t-1,n) + curr_w + u_partialAccess_feedback(i,t-1,n);
                x_partialAccess_feedback_linear(i,t,n) = alpha*x_partialAccess_feedback_linear(i,t-1,n) + curr_w + u_partialAccess_feedback_linear(i,t-1,n);
                
            end
            
            %% observer
            if t==1
                
                % in the first stage, we just normalize the current instant and
                % send to the controller
                a_zeroAccess = sqrt(P/W)*x_zeroAccess(i,t,n);
                a_fullAccess = sqrt(P/W)*x_fullAccess(i,t,n);
                a_partialAccess = sqrt(P/W)*x_partialAccess(i,t,n);
                a_zeroAccess_modulo = sqrt(P/W)*x_zeroAccess_modulo(i,t,n);
                a_partialAccess_feedback = sqrt(P/W)*x_partialAccess_feedback(i,t,n);
                a_partialAccess_feedback_linear = sqrt(P/W)*x_partialAccess_feedback_linear(i,t,n);
                
                P_X_linear(t,n) = W;
                P_X_modulo(t,n) = W;
                P_X_feedback(t,n) = W;
                P_X_partial(t,n) = W;
            else
                %% full access :
                % subtract the last stage receive side estimation, and normalize
                a_fullAccess = sqrt(P/P_error_predict(i,t,n))*(x_fullAccess(i,t,n) - x_hat_predict_fullAccess(i,t,n));
                
                %% Zero access , linear scheme : normalize to worst case power, rest is same as Kochman-Zamir
                P_X_linear(t,n) = W + P_X_linear(t-1,n) * (alpha - k_vec_up(t-1))^2 + (alpha^2 + (alpha - k_vec_up(t-1))^2) * P_error_estim_zeroAccess(i,t-1,n);                
                a_zeroAccess = sqrt(P/P_X_linear(t,n)) * x_zeroAccess(i,t,n);
                
                %% zero access with modulo : scale the unknown source part using "zooming" factor, and dither and send through the channel
                
                if mod(t-1,uncoded_timesharing_factor) == 0
                    P_X_modulo(t,n) = W + (k_vec(t-1))^2 * P_error_estim_zeroAccess_modulo(i,t-1,n);
                    P_X_hat_modulo(t,n) = P_X_modulo(t,n) + P_error_predict_zeroAccess_modulo(i,t-1,n);
                    curr_zoom = calc_zoom(P_alias,P,P_X_modulo(t,n),snrLin(i));
                else
                    curr_zoom = calc_zoom(P_alias,P,P_error_predict_zeroAccess_modulo(i,t,n),snrLin(i));
                end
                curr_dither = 2*delta*(rand - 0.5);
                
                % combine uncoded transmission sometime
                if mod(t,uncoded_timesharing_factor) == 0
                    a_zeroAccess_modulo = x_zeroAccess_modulo(i,t,n);
                else
                    a_zeroAccess_modulo = mod(curr_zoom*x_zeroAccess_modulo(i,t,n) + curr_dither + delta,2*delta) - delta;
                end
                
                %% partial access : estimate the noisy control, and then subtract noisy control estimate
                P_X_partial(t,n) = W + P_X_partial(t-1,n) * (alpha - k_vec_up(t-1))^2 + (alpha^2 + (alpha - k_vec_up(t-1))^2) * P_error_estim_partialAccess(i,t-1,n);
                P_X_hat_partial(t,n) = P_X_partial(t,n) + P_error_predict_partialAccess(i,t,n);
                
                r = sqrt(P/P_X_hat_partial(t,n))*x_hat_predict_partialAccess(i,t,n) + sigma_v*randn;
                estim_coeff = sqrt(P_X_hat_partial(t,n)/P);
                x_estim = estim_coeff*r;
                
                norm_partialAccess = P_error_predict_partialAccess(i,t,n) + estim_coeff^2 * sigma_v^2;
                a_partialAccess = sqrt(P/norm_partialAccess) * (x_partialAccess(i,t,n) - x_estim);
                
                %% partial access with feedback : estimate the noisy control after feedback iterations,then subtract noisy control estimate
                %  feedback is with modulo
                
                [x_final,MSE_final] = feedback_over_awgn(x_hat_predict_partialAccess_feedback(i,t,n),P,P_X_hat_partial(t,n),SNR(i) + deltaSNR,...
                    -1*deltaSNR,N_feedback,P_alias,'modulo');
                
                norm_partialAccess_feedback = P_error_predict_partialAccess_feedback(i,t,n) + MSE_final;
                a_partialAccess_feedback = sqrt(P/norm_partialAccess_feedback) * (x_partialAccess_feedback(i,t,n) - x_final);
                
                %% partial access with feedback : estimate the noisy control after feedback iterations,then subtract noisy control estimate
                %  feedback is with linear scaling
                
                [x_final_linear,MSE_final_linear] = feedback_over_awgn(x_hat_predict_partialAccess_feedback_linear(i,t,n),P,P_X_hat_partial(t,n),SNR(i) + deltaSNR,...
                    -1*deltaSNR,N_feedback,P_alias,'linear');
                
                norm_partialAccess_feedback_linear = P_error_predict_partialAccess_feedback_linear(i,t,n) + MSE_final_linear;
                a_partialAccess_feedback_linear = sqrt(P/norm_partialAccess_feedback_linear) * (x_partialAccess_feedback_linear(i,t,n) - x_final_linear);
            end
            
            %% channel
            z = sigma_z*randn;
            b_fullAccess = a_fullAccess + z;
            b_partialAccess = a_partialAccess + z;
            b_zeroAccess_modulo = a_zeroAccess_modulo + z;
            b_partialAccess_feedback = a_partialAccess_feedback + z;
            b_zeroAccess = a_zeroAccess + z;
            b_partialAccess_feedback_linear = a_partialAccess_feedback_linear + z;
            
            if debug
                % check the SNR
                snr_fullAccess(i,t) = a_fullAccess^2/sigma_z^2;
                snr_partialAccess(i,t) = a_partialAccess^2/sigma_z^2;
                snr_zeroAccess_modulo(i,t) = a_zeroAccess_modulo^2/sigma_z^2;
                snr_partialAccess_feedback(i,t) = a_partialAccess_feedback^2/sigma_z^2;
            end
            
            %% controller
            
            % estimation
            if t==1
                x_hat_estim_zeroAccess(i,t,n) = snrLin(i)/(1+snrLin(i))*sqrt(W/P)*b_zeroAccess;
                x_hat_estim_fullAccess(i,t,n) = snrLin(i)/(1+snrLin(i)) * sqrt(W/P) * b_fullAccess;
                x_hat_estim_partialAccess(i,t,n) = snrLin(i)/(1+snrLin(i)) * sqrt(W/P) * b_partialAccess;
                x_hat_estim_zeroAccess_modulo(i,t,n) = snrLin(i)/(1+snrLin(i)) * sqrt(W/P) * b_zeroAccess_modulo;
                x_hat_estim_partialAccess_feedback(i,t,n) = snrLin(i)/(1+snrLin(i)) * sqrt(W/P) * b_partialAccess_feedback;
                x_hat_estim_partialAccess_feedback_linear(i,t,n) = snrLin(i)/(1+snrLin(i)) * sqrt(W/P) * b_partialAccess_feedback_linear;
                
                P_error_estim(i,t,n) = W/(1+snrLin(i));
                P_error_estim_zeroAccess(i,t,n) = W/(1+snrLin(i));
                P_error_estim_partialAccess(i,t,n) = W/(1+snrLin(i));
                P_error_estim_zeroAccess_modulo(i,t,n) = W/(1+snrLin(i));
                P_error_estim_partialAccess_feedback(i,t,n) = W/(1+snrLin(i));
                P_error_estim_partialAccess_feedback_linear(i,t,n) = W/(1+snrLin(i));
            else
                
                %% full access : update current state estimation
                x_hat_estim_fullAccess(i,t,n) = x_hat_predict_fullAccess(i,t,n) + ...
                    snrLin(i)/(1+snrLin(i)) * sqrt(P_error_predict(i,t,n)/P) * b_fullAccess;
                P_error_estim(i,t,n) = P_error_predict(i,t,n)/(1+snrLin(i));
                
                %% zero access linear: recover a CUBE and update estimate
                zeta = (sqrt(P/P_X_linear(t,n)) / (P + sigma_z^2))*P_error_predict_zeroAccess(i,t,n);
                y_tilde = zeta*b_zeroAccess;
                x_hat_estim_zeroAccess(i,t,n) = x_hat_predict_zeroAccess(i,t,n) + y_tilde ;
                P_error_estim_zeroAccess(i,t,n) = (1 - zeta*sqrt(P/P_X_linear(t,n)))^2*P_error_predict_zeroAccess(i,t,n) + ...
                    zeta^2 * (sigma_z)^2; 
                
                %% partial access : update current state estimation
                gamma = sqrt(P/norm_partialAccess);
                partial_zeta = gamma*P_error_predict_partialAccess(i,t,n)/(P + (sigma_z)^2);
                x_hat_estim_partialAccess(i,t,n) = x_hat_predict_partialAccess(i,t,n) + partial_zeta * b_partialAccess;
                
                P_error_estim_partialAccess(i,t,n) = (1-gamma*partial_zeta)^2*P_error_predict_partialAccess(i,t,n) + ...
                    (partial_zeta*gamma*estim_coeff)^2 * (sigma_v)^2 + partial_zeta^2*(sigma_z)^2;
                %% zero access : modulo receiver
                
                % handle uncoded transmissions
                if mod(t,uncoded_timesharing_factor) == 0
                    b_tilde = b_zeroAccess_modulo;
                    
                    % update state estimate
                    x_hat_estim_zeroAccess_modulo(i,t,n) = b_tilde;
                    P_error_estim_zeroAccess_modulo(i,t,n) = sigma_z^2;
                    
                elseif mod(t-1,uncoded_timesharing_factor) == 0
                    alpha_mmse = snrLin(i)/(1 + snrLin(i));
                    b_tilde = mod(delta + alpha_mmse*b_zeroAccess_modulo - curr_dither,2*delta) - delta;
                    
                    % update state estimate
                    x_hat_estim_zeroAccess_modulo(i,t,n) = b_tilde/curr_zoom;
                    
                    P_error_estim_zeroAccess_modulo(i,t,n) = (P/snrLin(i))/curr_zoom^2;
                else
                    
                    % simple lmmse, assume there is no modulo-aliasing
                    b_tilde = b_zeroAccess_modulo - curr_dither - curr_zoom*(alpha - k_vec(t-1))*x_hat_estim_zeroAccess_modulo(i,t-1,n);
                    b_tilde_mod = mod(delta + b_tilde,2*delta) - delta;
                    
                    r_estim = curr_zoom*P_error_predict_zeroAccess_modulo(i,t,n)/(curr_zoom^2*P_error_predict_zeroAccess_modulo(i,t,n) + P/snrLin(i));
                    e_hat = r_estim * b_tilde_mod;
                    
                    %                     % update state estimate
                    x_hat_estim_zeroAccess_modulo(i,t,n) = x_hat_predict_zeroAccess_modulo(i,t,n) + e_hat;
                    
                    P_error_estim_zeroAccess_modulo(i,t,n) = (1 - curr_zoom^2*P_error_predict_zeroAccess_modulo(i,t,n)/(curr_zoom^2*P_error_predict_zeroAccess_modulo(i,t,n) + P/snrLin(i)))^2*P_error_predict_zeroAccess_modulo(i,t,n) + ...
                        (curr_zoom*P_error_predict_zeroAccess_modulo(i,t,n)/(curr_zoom^2*P_error_predict_zeroAccess_modulo(i,t,n) + P/snrLin(i)))^2 * P/(1 + snrLin(i));
                end
                
                %% partial access with feedback : we assume that the feedback output is perfect, so we just treat it
                %% as in the full access scheme
                x_hat_estim_partialAccess_feedback(i,t,n) = x_hat_predict_partialAccess_feedback(i,t,n) + ...
                    snrLin(i)/(1+snrLin(i)) * sqrt(norm_partialAccess_feedback/P) * b_partialAccess_feedback;
                P_error_estim_partialAccess_feedback(i,t,n) = P_error_predict_partialAccess_feedback(i,t,n)/(1+snrLin(i));

                %% partial access with feedback : we assume that the feedback output is perfect, so we just treat it
                %% as in the full access scheme
                x_hat_estim_partialAccess_feedback_linear(i,t,n) = x_hat_predict_partialAccess_feedback_linear(i,t,n) + ...
                    snrLin(i)/(1+snrLin(i)) * sqrt(norm_partialAccess_feedback_linear/P) * b_partialAccess_feedback_linear;
                P_error_estim_partialAccess_feedback_linear(i,t,n) = P_error_predict_partialAccess_feedback_linear(i,t,n)/(1+snrLin(i));
            end
            
            %% control generation
            u_zeroAccess(i,t,n) = -k_vec(t)*x_hat_estim_zeroAccess(i,t,n);
            u_partialAccess(i,t,n) = -k_vec(t)*x_hat_estim_partialAccess(i,t,n);
            u_fullAccess(i,t,n) = -k_vec(t)*x_hat_estim_fullAccess(i,t,n);
            u_partialAccess_feedback(i,t,n) = -k_vec(t)*x_hat_estim_partialAccess_feedback(i,t,n);
            u_partialAccess_feedback_linear(i,t,n) = -k_vec(t)*x_hat_estim_partialAccess_feedback_linear(i,t,n);
            
            if mod(t,uncoded_timesharing_factor) == 0
                u_zeroAccess_modulo(i,t,n) = -alpha*x_hat_estim_zeroAccess_modulo(i,t,n);
            else
                u_zeroAccess_modulo(i,t,n) = -k_vec(t)*x_hat_estim_zeroAccess_modulo(i,t,n);
            end
            %% update estimates
            x_hat_predict_fullAccess(i,t+1,n) = alpha*x_hat_estim_fullAccess(i,t,n) + u_fullAccess(i,t,n);
            P_error_predict(i,t+1,n) = alpha^2 * P_error_estim(i,t,n) + W;

            x_hat_predict_zeroAccess(i,t+1,n) = alpha*x_hat_estim_zeroAccess(i,t,n) + u_zeroAccess(i,t,n);
            P_error_predict_zeroAccess(i,t+1,n) = alpha^2*P_error_estim_zeroAccess(i,t,n) + W;
            
            x_hat_predict_partialAccess(i,t+1,n) = alpha*x_hat_estim_partialAccess(i,t,n) + u_partialAccess(i,t,n);
            P_error_predict_partialAccess(i,t+1,n) = alpha^2*P_error_estim_partialAccess(i,t,n) + W;
            
            x_hat_predict_zeroAccess_modulo(i,t+1,n) = alpha*x_hat_estim_zeroAccess_modulo(i,t,n) + u_zeroAccess_modulo(i,t,n);
            P_error_predict_zeroAccess_modulo(i,t+1,n) = alpha^2*P_error_estim_zeroAccess_modulo(i,t,n) + W;
            
            x_hat_predict_partialAccess_feedback(i,t+1,n) = alpha*x_hat_estim_partialAccess_feedback(i,t,n) + u_partialAccess_feedback(i,t,n);
            P_error_predict_partialAccess_feedback(i,t+1,n) = alpha^2*P_error_estim_partialAccess_feedback(i,t,n) + W;
            
            x_hat_predict_partialAccess_feedback_linear(i,t+1,n) = alpha*x_hat_estim_partialAccess_feedback_linear(i,t,n) + u_partialAccess_feedback_linear(i,t,n);
            P_error_predict_partialAccess_feedback_linear(i,t+1,n) = alpha^2*P_error_estim_partialAccess_feedback_linear(i,t,n) + W;
            
            %% Calculate current cost
            if t==1
                J_fullAccess(i,t,n) = (Q*(x_fullAccess(i,t,n))^2 + R*(u_fullAccess(i,t,n))^2);
                J_zeroAccess(i,t,n) = (Q*(x_zeroAccess(i,t,n))^2 + R*(u_zeroAccess(i,t,n))^2);
                J_partialAccess(i,t,n) = (Q*(x_partialAccess(i,t,n))^2 + R*(u_partialAccess(i,t,n))^2);
                J_zeroAccess_modulo(i,t,n) = (Q*(x_zeroAccess_modulo(i,t,n))^2 + R*(u_zeroAccess_modulo(i,t,n))^2);
                J_partialAccess_feedback(i,t,n) = (Q*(x_partialAccess_feedback(i,t,n))^2 + R*(u_partialAccess_feedback(i,t,n))^2);
                J_partialAccess_feedback_linear(i,t,n) = (Q*(x_partialAccess_feedback_linear(i,t,n))^2 + R*(u_partialAccess_feedback_linear(i,t,n))^2);
            else
                J_fullAccess(i,t,n) = (J_fullAccess(i,t-1,n)*(t-1) + (Q*(x_fullAccess(i,t,n))^2 + R*(u_fullAccess(i,t,n))^2))/t;
                J_zeroAccess(i,t,n) = (J_zeroAccess(i,t-1,n)*(t-1) + (Q*(x_zeroAccess(i,t,n))^2 + R*(u_zeroAccess(i,t,n))^2))/t;
                J_partialAccess(i,t,n) = (J_partialAccess(i,t-1,n)*(t-1) + (Q*(x_partialAccess(i,t,n))^2 + R*(u_partialAccess(i,t,n))^2))/t;
                J_zeroAccess_modulo(i,t,n) = (J_zeroAccess_modulo(i,t-1,n)*(t-1) + (Q*(x_zeroAccess_modulo(i,t,n))^2 + R*(u_zeroAccess_modulo(i,t,n))^2))/t;
                J_partialAccess_feedback(i,t,n) = (J_partialAccess_feedback(i,t-1,n)*(t-1) + (Q*(x_partialAccess_feedback(i,t,n))^2 + R*(u_partialAccess_feedback(i,t,n))^2))/t;
                J_partialAccess_feedback_linear(i,t,n) = (J_partialAccess_feedback_linear(i,t-1,n)*(t-1) + (Q*(x_partialAccess_feedback_linear(i,t,n))^2 + R*(u_partialAccess_feedback_linear(i,t,n))^2))/t;
            end
        end
    end
    
    figure;
    % generate sample path plots, and cost
    hold all
    plot(1:T,mean(x_fullAccess(i,:,:),3),'-g','LineWidth',1.5)
    plot(1:T,mean(x_zeroAccess(i,:,:),3),'-k','LineWidth',1.5)
    plot(1:T,mean(x_partialAccess(i,:,:),3),'-r','LineWidth',1.5)
    plot(1:T,mean(x_zeroAccess_modulo(i,:,:),3),'--c','LineWidth',1.5)
    plot(1:T,mean(x_partialAccess_feedback(i,:,:),3),'-b','LineWidth',1.5)
    plot(1:T,mean(x_partialAccess_feedback_linear(i,:,:),3),'-m','LineWidth',1.5)
    
    xlabel('t'); ylabel('state')
    legend('x(t) for full control access','x(t) for zero control access','x(t) for partial control access',...
        'x(t) for zero control access, with Kochman-Zamir',...
        strcat('x(t) for partial control access, with feedback :',num2str(N_feedback),' iterations'),...
        strcat('x(t) for partial control access, with linear feedback :',num2str(N_feedback),' iterations'))
    title(strcat('state , snr = ',num2str(SNR(i)),'[dB]',' \Deltasnr = ',num2str(deltaSNR)))
    grid on; grid minor;
    
    figure; hold all
    plot(1:T,10*log10(mean(J_fullAccess(i,:,:),3)),'-g','LineWidth',2)
    plot(1:T,10*log10(mean(J_zeroAccess(i,:,:),3)),'--k','LineWidth',2)
    plot(1:T,10*log10(mean(J_partialAccess(i,:,:),3)),'-.r','LineWidth',2)
%     plot(1:T,10*log10(mean(J_zeroAccess_modulo(i,:,:),3)),':c','LineWidth',2)
    plot(1:T,10*log10(mean(J_partialAccess_feedback(i,:,:),3)),'b','LineWidth',2)
    plot(1:T,10*log10(mean(J_partialAccess_feedback_linear(i,:,:),3)),':m','LineWidth',2)
    grid on; grid minor;
    xlabel('t'); ylabel('cost [dB]');
    title(strcat('Cost in [dB] over time , snr = ',num2str(SNR(i)),'[dB]',' \Deltasnr = ',...
        num2str(deltaSNR),'Q = ',num2str(Q),'R = ',num2str(R),' \alpha = ',num2str(alpha)));
    legend('Full Access','Zero Access','Partial Access',...
        'Zero Access - kochman zamir',strcat('Partial Access - feedback : ',num2str(N_feedback),' iterations'),...
        strcat('Partial Access - linear feedback : ',num2str(N_feedback),' iterations'));
%     legend('Full Access','Zero Access','Partial Access',...
%      strcat('Partial Access - feedback : ',num2str(N_feedback),' iterations'),...
%         strcat('Partial Access - linear feedback : ',num2str(N_feedback),' iterations'));
%     

end
figure;
subplot(121);hold all
plot(1:T,10*log10(mean(J_fullAccess(1,:,:),3)),'-gs','LineWidth',1.5)
plot(1:T,10*log10(mean(J_partialAccess(1,:,:),3)),'-r*','LineWidth',1.5)
plot(1:T,10*log10(mean(J_zeroAccess_modulo(1,:,:),3)),'-cd','LineWidth',1.5)
plot(1:T,10*log10(mean(J_partialAccess_feedback(1,:,:),3)),'-b>','LineWidth',1.5)
grid on; grid minor;
xlabel('t'); ylabel('cost [dB]');
title(strcat('J_{t} [dB] , snr = ',num2str(SNR(1)),'[dB]','Q = ',num2str(Q),', R = ',num2str(R),'\alpha = ',num2str(alpha)));
legend('Full Access','Partial Access',...
    'Zero Access - kochman zamir',strcat('Partial Access - N_{feedback} : ',num2str(N_feedback)));

subplot(122);hold all
plot(1:T,10*log10(mean(J_fullAccess(4,:,:),3)),'-gs','LineWidth',1.5)
plot(1:T,10*log10(mean(J_partialAccess(4,:,:),3)),'-r*','LineWidth',1.5)
plot(1:T,10*log10(mean(J_zeroAccess_modulo(4,:,:),3)),'-cd','LineWidth',1.5)
plot(1:T,10*log10(mean(J_partialAccess_feedback(4,:,:),3)),'-b>','LineWidth',1.5)
grid on; grid minor;
xlabel('t'); ylabel('cost [dB]');
title(strcat('J_{t} [dB] , snr = ',num2str(SNR(4)),'[dB]','Q = ',num2str(Q),', R = ',num2str(R),'\alpha = ',num2str(alpha)));
legend('Full Access','Partial Access',...
    'Zero Access - kochman zamir',strcat('Partial Access - N_{feedback} : ',num2str(N_feedback)));
function [k_vec] = calcLQG(Q,R,T,alpha)

% calculate k
k_vec = zeros(1,T-1);
s_vec = zeros(1,T-1);
for t=T:-1:2
    if t==T
        s_vec(t) = Q;
        k_vec(t-1) = alpha*s_vec(t)/(s_vec(t) + R);
    else
        s_vec(t) = ((alpha^2)*R*s_vec(t+1))/(s_vec(t+1) + R) + Q;
        k_vec(t-1) = alpha*s_vec(t)/(s_vec(t) + R);
    end
end
end

function [curr_zoom] = calc_zoom(P_alias,P,P_predict,snrLin)

a = 3*P / ((qfuncinv(P_alias/2))^2);
b = P/snrLin;
c = P_predict;

curr_zoom = sqrt((a - b)/c);

end

function [k_hat] = ML_interval_index(b_zeroAccess_modulo,r_estim,e_hat_frac,curr_zoom,curr_dither,delta,alpha,L_t,x_hat_estim_zeroAccess_modulo)
k_grid = -3:1:3;
metric = abs(b_zeroAccess_modulo - (mod(delta + curr_zoom*(e_hat_frac - 2*delta*r_estim*k_grid + ...
    (alpha-L_t)*x_hat_estim_zeroAccess_modulo) + curr_dither,2*delta) - delta)).^2;

[~,k_ind] = min(metric);
k_hat = k_grid(k_ind);
end