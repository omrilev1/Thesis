%% Simulation of Markov Source tracking over power limited awgn channel
% we simulate 3 cases :
% - sensor and observer has full access to the last state estimation \hat{s}_t|t
% - only sensor has the estimation
% - noisy feedback between both

% The evolution of s_t is given by :
% s_(t+1) = lambda*s_t + w_(t+1)
% where w_(t+1) is wgn with variance W

% The goal is to minimize the estimation error ||s_{t} - \hat{s}_{t|t}||^2
clear all; clc; close all

%% coefficeints initialization
% model coefficients
lambda = 0.98;
W = 1; % input gaussian variance
T = 2^9; % horizon
P = 1; % power constraint normalized to 1
P_ss = W/(1-lambda^2);
P_timesharing = P;

delta = sqrt(3*P);

sign_cheat = 0; % to transmit in the feedback path only the absolute value, with half of the power

SNR = 14; % [20 15 10 8 6];
deltaSNR = 8;
snrLin = 10.^(SNR/10);
snrFeedback_Lin = 10.^((SNR + deltaSNR)/10);

N_avg = 256;

initArrays;

P_alias = 9*1e-3;
N_feedback = 4;
% set timesharing parameters - we need to scale the powers to the overall
% avergae power will be P
uncoded_timesharing_factor = 2^8;

% init Tuncel parameters
delta_Tuncel = 3.6;% 2.6
partition = (-6.5:1:6.5)*delta_Tuncel;
codebook = generateGaussianCodeBook(partition);
alpha = 0.6;% 0.6;
beta = -1.165;% -1.165;
centers = [partition(1) (partition(1:end-1) + partition(2:end))/2 partition(end)];


for i=1:length(SNR)
    sigma_z = sqrt(P/snrLin(i));
    sigma_v = sqrt(P/snrFeedback_Lin(i));
    
    P_X_linear = zeros(T,N_avg);
    P_X_hat_linear = zeros(T,N_avg);
    
    P_X_partial = zeros(T,N_avg);
    P_X_hat_partial = zeros(T,N_avg);
    
    P_X_modulo_feedback = zeros(T,N_avg);
    P_X_modulo_feedback_linear = zeros(T,N_avg);
    
    P_X_full_debug = zeros(T,N_avg);
    P_X_linear_debug = zeros(T,N_avg);
    P_X_partial_debug = zeros(T,N_avg);
    P_X_partial_feedback_linear_debug = zeros(T,N_avg);
    P_X_partial_feedback_debug = zeros(T,N_avg);
    
    P_X_modulo = zeros(T,N_avg);
    P_X_modulo_debug = zeros(T,N_avg);
    P_X_hat_modulo = zeros(T,N_avg);
    
    P_X_Tuncel = zeros(T,N_avg);
    P_X_Tuncel_debug = zeros(T,N_avg);
    P_X_Tuncel_modulo = zeros(T,N_avg);
    
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
                x_zeroAccess_modulo(i,t,n) = curr_w;
                x_zeroAccess_Tuncel(i,t,n) = curr_w;
                x_partialAccess_feedback(i,t,n) = curr_w;
                x_partialAccess_feedback_linear(i,t,n) = curr_w;
                
            elseif t==T
                curr_w = sqrt(W)*randn;
                % no control in the last stage - evolve state and break
                x_fullAccess(i,t,n) = lambda*x_fullAccess(i,t-1,n) + curr_w;
                J_fullAccess(i,t,n) = ((x_fullAccess(i,t,n) - x_hat_predict_fullAccess(i,t,n))^2);
                
                % no control in the last stage - evolve state and break
                x_zeroAccess(i,t,n) = lambda*x_zeroAccess(i,t-1,n) + curr_w;
                J_zeroAccess(i,t,n) = (x_zeroAccess(i,t,n) - x_hat_predict_zeroAccess(i,t,n))^2;
                
                % no control in the last stage - evolve state and break
                x_zeroAccess_modulo(i,t,n) = lambda*x_zeroAccess_modulo(i,t-1,n) + curr_w;
                J_zeroAccess_modulo(i,t,n) = (x_zeroAccess_modulo(i,t,n) - x_hat_predict_zeroAccess_modulo(i,t,n))^2;
                
                % no control in the last stage - evolve state and break
                x_zeroAccess_Tuncel(i,t,n) = lambda*x_zeroAccess_Tuncel(i,t-1,n) + curr_w;
                J_zeroAccess_Tuncel(i,t,n) = (x_zeroAccess_Tuncel(i,t,n) - x_hat_predict_zeroAccess_Tuncel(i,t,n))^2;
                
                % no control in the last stage - evolve state and break
                x_partialAccess_feedback(i,t,n) = lambda*x_partialAccess_feedback(i,t-1,n) + curr_w;
                J_partialAccess_feedback(i,t,n) = (x_partialAccess_feedback(i,t,n) - x_hat_predict_partialAccess_feedback(i,t,n))^2;
                
                
                % no control in the last stage - evolve state and break
                x_partialAccess_feedback_linear(i,t,n) = lambda*x_partialAccess_feedback_linear(i,t-1,n) + curr_w;
                J_partialAccess_feedback_linear(i,t,n) = (x_partialAccess_feedback_linear(i,t,n) - x_hat_predict_partialAccess_feedback_linear(i,t,n))^2;
                break
                
                
            else
                % evolution
                curr_w = sqrt(W)*randn;
                x_zeroAccess(i,t,n) = lambda*x_zeroAccess(i,t-1,n) + curr_w;
                x_fullAccess(i,t,n) = lambda*x_fullAccess(i,t-1,n) + curr_w;
                x_zeroAccess_modulo(i,t,n) = lambda*x_zeroAccess_modulo(i,t-1,n) + curr_w;
                x_zeroAccess_Tuncel(i,t,n) = lambda*x_zeroAccess_Tuncel(i,t-1,n) + curr_w;
                x_partialAccess_feedback(i,t,n) = lambda*x_partialAccess_feedback(i,t-1,n) + curr_w;
                x_partialAccess_feedback_linear(i,t,n) = lambda*x_partialAccess_feedback_linear(i,t-1,n) + curr_w;
                
            end
            
            %% observer
            if t==1
                
                % in the first stage, we just normalize the current instant and
                % send to the controller
                a_zeroAccess = sqrt(P/W)*x_zeroAccess(i,t,n);
                a_fullAccess = sqrt(P/W)*x_fullAccess(i,t,n);
                a_zeroAccess_modulo = sqrt(P/W)*x_zeroAccess_modulo(i,t,n);
                a_zeroAccess_Tuncel = sqrt(P/W)*x_zeroAccess_Tuncel(i,t,n);
                a_partialAccess_feedback = sqrt(P/W)*x_partialAccess_feedback(i,t,n);
                a_partialAccess_feedback_linear = sqrt(P/W)*x_partialAccess_feedback_linear(i,t,n);
                
                P_X_linear(t,n) = W;
                P_X_modulo(t,n) = W;
                P_X_feedback(t,n) = W;
                P_X_Tuncel(t,n) = W;
            else
                %% full access :
                % subtract the last stage receive side estimation, and normalize
                a_fullAccess = sqrt(P/P_error_predict(i,t,n))*(x_fullAccess(i,t,n) - x_hat_predict_fullAccess(i,t,n));
                P_X_full_debug(t,n) = a_fullAccess.^2;
                
                %% Zero access , linear scheme : normalize according to some strategy, rest is same as Kochman-Zamir
                P_X_linear(t,n) = W + lambda^2 * P_X_linear(t-1,n);
                a_zeroAccess = sqrt(P/P_X_linear(t,n)) * x_zeroAccess(i,t,n);
                P_X_linear_debug(t,n) = a_zeroAccess.^2;
                
                %% zero access with modulo : scale the unknown source part using "zooming" factor, and dither and send through the channel
                
                P_X_modulo(t,n) = W + lambda^2*P_X_modulo(t-1,n);
                %                 curr_zoom = 1.15*sqrt(P/P_X_modulo(t,n));
                curr_zoom = calc_zoom(P_alias,P,P_error_predict_zeroAccess_modulo(i,t,n),snrLin(i));
                curr_dither = 0; % 2*delta*(rand - 0.5);
                
                % combine uncoded transmission sometime
                if mod(t,uncoded_timesharing_factor) == 0
                    a_zeroAccess_modulo = sqrt(P_timesharing/P_X_modulo(t,n))*x_zeroAccess_modulo(i,t,n);
                else
                    %                     a_zeroAccess_modulo = mod(curr_zoom*x_zeroAccess_modulo(i,t,n) + delta,2*delta) - delta;
                    a_zeroAccess_modulo = mod(curr_zoom*x_zeroAccess_modulo(i,t,n) + curr_dither + delta,2*delta) - delta;
                    
                end
                P_X_modulo_debug(t,n) = a_zeroAccess_modulo.^2;
                
                %% Tuncel scheme : quantize and send the quantized value + quantization value after proper scaling
                
                P_X_Tuncel(t,n) = W + lambda^2*P_X_Tuncel(t-1,n);
                
                % Tuncel coding
                normX_Tuncel = sqrt(P/P_X_Tuncel(t,n)) * x_zeroAccess_Tuncel(i,t,n);
                currDistance = abs(normX_Tuncel - centers(:));
                [~,minIdx] = min(currDistance,[],1);
                T_Tuncel = codebook(minIdx);
                S_Tuncel = normX_Tuncel - T_Tuncel;
                
                a_zeroAccess_Tuncel = alpha*T_Tuncel + beta*S_Tuncel;
                P_X_Tuncel_debug(t,n) = a_zeroAccess_Tuncel.^2;
                
                %% partial access with feedback : estimate the noisy control after feedback iterations,then subtract noisy control estimate
                %  feedback is with modulo
                P_X_modulo_feedback(t,n) = W + lambda*P_X_modulo_feedback(t-1,n);
                
                [a_partialAccess_feedback,MSE_final] = feedback_over_awgn(x_partialAccess_feedback(i,t,n),P,P_X_modulo_feedback(t,n),SNR(i),...
                    deltaSNR,N_feedback,1e-6,'modulo');
                
                %% partial access with feedback : estimate the noisy control after feedback iterations,then subtract noisy control estimate
                %  feedback is with linear scaling
                P_X_modulo_feedback_linear(t,n) = W + lambda*P_X_modulo_feedback_linear(t-1,n);
                
                [a_partialAccess_feedback_linear,MSE_final_linear] = feedback_over_awgn(x_partialAccess_feedback_linear(i,t,n),P,P_X_modulo_feedback_linear(t,n),SNR(i),...
                    deltaSNR,N_feedback,1e-6,'linear');
            end
            
            %% channel
            z = sigma_z*randn;
            b_fullAccess = a_fullAccess + z;
            b_zeroAccess_modulo = a_zeroAccess_modulo + z;
            b_zeroAccess = a_zeroAccess + z;
            b_zeroAccess_Tuncel = a_zeroAccess_Tuncel + z;
            
            if t==1
                b_partialAccess_feedback = a_partialAccess_feedback + z;
                b_partialAccess_feedback_linear = a_partialAccess_feedback_linear + z;
            else
                b_partialAccess_feedback = a_partialAccess_feedback;
                b_partialAccess_feedback_linear = a_partialAccess_feedback_linear;
            end
            %% sensor
            
            % estimation
            if t==1
                x_hat_estim_zeroAccess(i,t,n) = snrLin(i)/(1+snrLin(i))*sqrt(W/P)*b_zeroAccess;
                x_hat_estim_fullAccess(i,t,n) = snrLin(i)/(1+snrLin(i)) * sqrt(W/P) * b_fullAccess;
                x_hat_estim_zeroAccess_modulo(i,t,n) = snrLin(i)/(1+snrLin(i)) * sqrt(W/P) * b_zeroAccess_modulo;
                x_hat_estim_zeroAccess_Tuncel(i,t,n) = snrLin(i)/(1+snrLin(i)) * sqrt(W/P) * b_zeroAccess_Tuncel;
                x_hat_estim_partialAccess_feedback(i,t,n) = snrLin(i)/(1+snrLin(i)) * sqrt(W/P) * b_partialAccess_feedback;
                x_hat_estim_partialAccess_feedback_linear(i,t,n) = snrLin(i)/(1+snrLin(i)) * sqrt(W/P) * b_partialAccess_feedback_linear;
                
                P_error_estim(i,t,n) = W/(1+snrLin(i));
                P_error_estim_zeroAccess(i,t,n) = W/(1+snrLin(i));
                P_error_estim_zeroAccess_modulo(i,t,n) = W/(1+snrLin(i));
                P_error_estim_zeroAccess_Tuncel(i,t,n) = W/(1+snrLin(i));
                P_error_estim_partialAccess_feedback(i,t,n) = W/(1+snrLin(i));
                P_error_estim_partialAccess_feedback_linear(i,t,n) = W/(1+snrLin(i));
            else
                
                %% full access : update current state estimation
                x_hat_estim_fullAccess(i,t,n) = x_hat_predict_fullAccess(i,t,n) + ...
                    snrLin(i)/(1+snrLin(i)) * sqrt(P_error_predict(i,t,n)/P) * b_fullAccess;
                P_error_estim(i,t,n) = P_error_predict(i,t,n)/(1+snrLin(i));
                
                %% zero access linear: recover a CUBE and update estimate
                
                zeta = sqrt(P/P_X_linear(t,n));
                xi = zeta *P_error_predict_zeroAccess(i,t,n) / (zeta^2 * P_error_predict_zeroAccess(i,t,n) + sigma_z^2);
                y_tilde = xi*(b_zeroAccess - zeta * lambda * x_hat_estim_zeroAccess(i,t-1,n));
                x_hat_estim_zeroAccess(i,t,n) = x_hat_predict_zeroAccess(i,t,n) + y_tilde ;
                P_error_estim_zeroAccess(i,t,n) = (1 - zeta*xi)^2*P_error_predict_zeroAccess(i,t,n) + ...
                    xi^2 * (sigma_z)^2;
                
                %% zero access : modulo receiver
                
                % handle uncoded transmissions
                if mod(t,uncoded_timesharing_factor) == 0
                    
                    zeta_modulo = sqrt(P_timesharing/P_X_modulo(t,n));
                    xi_modulo = zeta_modulo *P_error_predict_zeroAccess_modulo(i,t,n) / (zeta_modulo^2 * P_error_predict_zeroAccess_modulo(i,t,n) + sigma_z^2);
                    y_tilde_modulo = xi_modulo*(b_zeroAccess_modulo - zeta_modulo * lambda * x_hat_estim_zeroAccess_modulo(i,t-1,n));
                    x_hat_estim_zeroAccess_modulo(i,t,n) = x_hat_predict_zeroAccess_modulo(i,t,n) + y_tilde_modulo ;
                    P_error_estim_zeroAccess_modulo(i,t,n) = (1 - zeta_modulo*xi_modulo)^2*P_error_predict_zeroAccess_modulo(i,t,n) + ...
                        xi_modulo^2 * (sigma_z)^2;
                    
                    
                else
                    
                    b_tilde = b_zeroAccess_modulo - curr_dither - curr_zoom*x_hat_predict_zeroAccess_modulo(i,t,n);
                    b_tilde_mod = mod(delta + b_tilde,2*delta) - delta;
                    
                    r_estim = curr_zoom*P_error_predict_zeroAccess_modulo(i,t,n)/(curr_zoom^2*P_error_predict_zeroAccess_modulo(i,t,n) + sigma_z^2);
                    e_hat = r_estim * b_tilde_mod;
                    
                    % update state estimate
                    x_hat_estim_zeroAccess_modulo(i,t,n) = x_hat_predict_zeroAccess_modulo(i,t,n) + e_hat;
                    
                    P_error_estim_zeroAccess_modulo(i,t,n) = (1 - curr_zoom^2*P_error_predict_zeroAccess_modulo(i,t,n)/(curr_zoom^2*P_error_predict_zeroAccess_modulo(i,t,n) + sigma_z^2))^2*P_error_predict_zeroAccess_modulo(i,t,n) + ...
                        (r_estim)^2 * sigma_z^2;
                    
                    
                    %                     zeta_modulo = curr_zoom;
                    %                     b_tilde = b_zeroAccess_modulo - zeta_modulo*lambda*x_hat_estim_zeroAccess_modulo(i,t-1,n);
                    %                     b_tilde_mod = mod(delta + b_tilde,2*delta) - delta;
                    %
                    %                     r_estim = zeta_modulo*P_error_predict_zeroAccess_modulo(i,t,n)/(zeta_modulo^2*P_error_predict_zeroAccess_modulo(i,t,n) + sigma_z^2);
                    %                     e_hat = r_estim * b_tilde_mod;
                    %
                    %                     % update state estimate
                    %                     x_hat_estim_zeroAccess_modulo(i,t,n) = x_hat_predict_zeroAccess_modulo(i,t,n) + e_hat;
                    %
                    %                     P_error_estim_zeroAccess_modulo(i,t,n) = (1 - zeta_modulo^2*P_error_predict_zeroAccess_modulo(i,t,n)/(zeta_modulo^2*P_error_predict_zeroAccess_modulo(i,t,n) + sigma_z^2))^2*P_error_predict_zeroAccess_modulo(i,t,n) + ...
                    %                         (zeta_modulo*P_error_predict_zeroAccess_modulo(i,t,n)/(zeta_modulo^2*P_error_predict_zeroAccess_modulo(i,t,n) + sigma_z^2))^2 * sigma_z^2;
                end
                
                %% Zero Access : Tuncel Receiver
                
                % start with parameter estimation
                rho = lambda^2*P_X_Tuncel(t,n) / (lambda^2*P_X_Tuncel(t,n) + W);
                y = sqrt(P/P_X_Tuncel(t,n))*x_hat_predict_zeroAccess_Tuncel(i,t,n);

                T_hat = MMSE_decoder_Int(b_zeroAccess_Tuncel,y,...
                    delta_Tuncel,codebook,alpha,beta,rho,sigma_z);
                S_hat = MMSE_decoder_Frac(b_zeroAccess_Tuncel,y,...
                    delta_Tuncel,alpha,beta,rho,sigma_z,codebook);
                x_hat_estim_zeroAccess_Tuncel(i,t,n) = sqrt(P_X_Tuncel(t,n)/P) * (T_hat + S_hat);
                %% partial access with feedback : we assume that the feedback output is perfect, so we just treat it
                %% as in the full access scheme
                
                xi = P_error_predict_partialAccess_feedback(i,t,n) / (P_error_predict_partialAccess_feedback(i,t,n) + MSE_final);
                y_tilde = xi*(b_partialAccess_feedback - lambda*x_hat_estim_partialAccess_feedback(i,t-1,n));
                x_hat_estim_partialAccess_feedback(i,t,n) = x_hat_predict_partialAccess_feedback(i,t,n) + y_tilde ;
                P_error_estim_partialAccess_feedback(i,t,n) = (1 - xi)^2*P_error_predict_partialAccess_feedback(i,t,n) + ...
                    xi^2 * MSE_final;
                
                %% partial access with feedback : we assume that the feedback output is perfect, so we just treat it
                %% as in the full access scheme
                
                
                xi = P_error_predict_partialAccess_feedback_linear(i,t,n) / (P_error_predict_partialAccess_feedback_linear(i,t,n) + MSE_final_linear);
                y_tilde = xi*(b_partialAccess_feedback_linear - lambda*x_hat_estim_partialAccess_feedback_linear(i,t-1,n));
                x_hat_estim_partialAccess_feedback_linear(i,t,n) = x_hat_predict_partialAccess_feedback_linear(i,t,n) + y_tilde ;
                P_error_estim_partialAccess_feedback_linear(i,t,n) = (1 - xi)^2*P_error_predict_partialAccess_feedback_linear(i,t,n) + ...
                    xi^2 * MSE_final_linear;
                
            end
            
            %% update estimates
            x_hat_predict_fullAccess(i,t+1,n) = lambda*x_hat_estim_fullAccess(i,t,n);
            P_error_predict(i,t+1,n) = lambda^2 * P_error_estim(i,t,n) + W;
            
            x_hat_predict_zeroAccess(i,t+1,n) = lambda*x_hat_estim_zeroAccess(i,t,n);
            P_error_predict_zeroAccess(i,t+1,n) = lambda^2*P_error_estim_zeroAccess(i,t,n) + W;
            
            x_hat_predict_zeroAccess_modulo(i,t+1,n) = lambda*x_hat_estim_zeroAccess_modulo(i,t,n);
            P_error_predict_zeroAccess_modulo(i,t+1,n) = lambda^2*P_error_estim_zeroAccess_modulo(i,t,n) + W;
            
            x_hat_predict_zeroAccess_Tuncel(i,t+1,n) = lambda*x_hat_estim_zeroAccess_Tuncel(i,t,n);
            
            x_hat_predict_partialAccess_feedback(i,t+1,n) = lambda*x_hat_estim_partialAccess_feedback(i,t,n);
            P_error_predict_partialAccess_feedback(i,t+1,n) = lambda^2*P_error_estim_partialAccess_feedback(i,t,n) + W;
            
            x_hat_predict_partialAccess_feedback_linear(i,t+1,n) = lambda*x_hat_estim_partialAccess_feedback_linear(i,t,n);
            P_error_predict_partialAccess_feedback_linear(i,t+1,n) = lambda^2*P_error_estim_partialAccess_feedback_linear(i,t,n) + W;
            
            %% Calculate current cost
            J_fullAccess(i,t,n) = (x_fullAccess(i,t,n) - x_hat_estim_fullAccess(i,t,n))^2 ;
            J_zeroAccess(i,t,n) = (x_zeroAccess(i,t,n) - x_hat_estim_zeroAccess(i,t,n))^2;
            J_zeroAccess_modulo(i,t,n) = (x_zeroAccess_modulo(i,t,n) - x_hat_estim_zeroAccess_modulo(i,t,n))^2;
            J_zeroAccess_Tuncel(i,t,n) = (x_zeroAccess_Tuncel(i,t,n) - x_hat_estim_zeroAccess_Tuncel(i,t,n))^2;
            J_partialAccess_feedback(i,t,n) = (x_partialAccess_feedback(i,t,n) - x_hat_estim_partialAccess_feedback(i,t,n))^2;
            J_partialAccess_feedback_linear(i,t,n) = (x_partialAccess_feedback_linear(i,t,n) - x_hat_estim_partialAccess_feedback_linear(i,t,n))^2;
            
        end
    end
    
    figure; hold all
    plot(1:(T-1),10*log10(mean(J_fullAccess(i,1:end-1,:),3)),'-g','LineWidth',2)
    plot(1:(T-1),10*log10(mean(J_zeroAccess(i,1:end-1,:),3)),'--k','LineWidth',2)
    plot(1:(T-1),10*log10(mean(J_zeroAccess_modulo(i,1:end-1,:),3)),':c','LineWidth',2)
    plot(1:(T-1),10*log10(mean(J_zeroAccess_Tuncel(i,1:end-1,:),3)),':g','LineWidth',2)
    plot(1:(T-1),10*log10(mean(J_partialAccess_feedback(i,1:end-1,:),3)),'b','LineWidth',2)
    plot(1:(T-1),10*log10(mean(J_partialAccess_feedback_linear(i,1:end-1,:),3)),':m','LineWidth',2)
    grid on; grid minor;
    xlabel('t'); ylabel('cost [dB]');
    %     title(strcat('Cost in [dB] over time , snr = ',num2str(SNR(i)),'[dB]',' \Deltasnr = ',...
    %         num2str(deltaSNR),'[dB] Q = ',num2str(Q),'R = ',num2str(R),' \alpha = ',num2str(alpha)));
    legend('Full Access','Zero Access','Zero Access - kochman zamir','Zero Access - Tuncel',...
        strcat('Modulo feedback : ',num2str(N_feedback),' iterations'),...
        strcat('Linear feedback : ',num2str(N_feedback),' iterations'));
    %     ylim([9 20.5])
    powerAnalysis(P_X_Tuncel_debug,N_avg,'Tuncel',sigma_z)
    powerAnalysis(P_X_modulo_debug,N_avg,'Kochman',sigma_z)
    
end

clear all;

function [curr_zoom] = calc_zoom(P_alias,P,P_predict,snrLin)

a = 3*P / ((qfuncinv(P_alias/2))^2);
b = P/snrLin;
c = P_predict;

curr_zoom = sqrt((a - b)/c);

end

function [] = powerAnalysis(P_x,N_avg,run,sigma)

inputPower = sum(P_x,2)/N_avg;

figure;hold all
plot(inputPower)
title([run strcat('measured SNR = ',num2str(10*log10(mean(inputPower(2:end-1))/sigma^2)))])

end

function [codebook] = generateGaussianCodeBook(partition)

dx = 0.001;
x = (-3 + partition(1)):dx:(partition(end) + 3);

pdf = (1/sqrt(2*pi))*exp(-0.5*x.^2);
codebook = zeros(1,length(partition) - 1);
for i=1:(length(partition) - 1)
    
    curIdx = find((x >= partition(i)) & (x <= partition(i+1)));
    curPDF = pdf(curIdx) / sum(pdf(curIdx));
    curX = x(curIdx);
    codebook(i) = sum(curX.*curPDF);
end

% generate first and last codewords

% first
curIdx = find(x <= partition(1));
curPDF = pdf(curIdx) / sum(pdf(curIdx));
curX = x(curIdx);
codebook = [sum(curX.*curPDF) codebook];

% last codeword
curIdx = find(x >= partition(end));
curPDF = pdf(curIdx) / sum(pdf(curIdx));
curX = x(curIdx);
codebook = [codebook sum(curX.*curPDF)];
end

function [T_hat] = MMSE_decoder_Int(RxSig_tuncel,y,delta,codebook,alpha,beta,rho,sigmaW)

ds = 0.001;
s = -delta/2:ds:delta/2;

% start with 2D calculations 
calc1 = codebook(:) + s;
calc2 = alpha*codebook(:) + beta*s;

expTerm1 = -1*(0.5*(calc1).^2);
expTerm2 = - (0.5/(1-rho^2))*(y - rho*calc1).^2;
expTerm3 = 0.5*(1/sigmaW^2)*(RxSig_tuncel - calc2).^2; 

Integrand = exp(expTerm1 + expTerm2 - expTerm3);

num = sum(Integrand,2);
T_hat = sum(codebook(:).*num(:),1)./sum(num(:));

% quantize to the nearset r_k
% [~,minIdx] = min(abs(T_hat - codebook));
% T_hat = codebook(minIdx);
end

function [S_hat] = MMSE_decoder_Frac(RxSig_tuncel,y,delta,alpha,beta,rho,sigmaW,codebook)


ds = 0.001;
s = ones(length(codebook),1) * (-delta/2:ds:delta/2);

% start with 2D calculations 
calc1 = codebook(:) + s;
calc2 = alpha*codebook(:) + beta*s;

expTerm1 = -1*(0.5*(calc1).^2);
expTerm2 = - (0.5/(1-rho^2))*(y - rho*calc1).^2;
expTerm3 = 0.5*(1/sigmaW^2)*(RxSig_tuncel - calc2).^2; 
Integrand = exp(expTerm1 + expTerm2 - expTerm3);


num = sum(sum(Integrand.*s,2),1);
den = sum(sum(Integrand));
S_hat = num/den;

end
