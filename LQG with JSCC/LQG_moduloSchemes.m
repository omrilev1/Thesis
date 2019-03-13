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
%% coefficeints initialization
% model coefficients
alpha = 2;
W = 1;
T = 100; % horizon
delta = 3;
P = delta^2/3;

p_alias = 1e-3;
% cost parameters
Q = 1;
R = 1;


SNR = [20 15 10 5];
snrLin = 10.^(SNR/10);

initArrays; 

cnt = 0;
for i=1:length(SNR)
    
    sigma_z = sqrt(P/10^(SNR(i)/10));
    
    % calculate LQG parameters for current run
    [gamma_vec,k_vec] = calcLQG(sigma_z,Q,R,W,T,delta,alpha,p_alias);
    
    for t=1:T
        
        %% plant
        if t==1
            % initialization
            curr_w = sqrt(W)*randn;
            x_fullAccess(i,t) = curr_w;
            x_zeroAccess(i,t) = curr_w;
            x_zeroAccess_lin(i,t) = curr_w;
            
        elseif t==T
            curr_w = sqrt(W)*randn;
            % no control in the last stage - evolve state and break
            x_fullAccess(i,t) = alpha*x_fullAccess(i,t-1) + curr_w + u_fullAccess(i,t-1);
            J_fullAccess(i,t) = J_fullAccess(i,t-1) + (Q*(x_fullAccess(i,t))^2)/T;
            
            % no control in the last stage - evolve state and break
            x_zeroAccess(i,t) = alpha*x_zeroAccess(i,t-1) + curr_w + u_zeroAccess(i,t-1);
            J_zeroAccess(i,t) = J_zeroAccess(i,t-1) + (Q*(x_zeroAccess(i,t))^2)/T;
            
            % no control in the last stage - evolve state and break
            x_zeroAccess_lin(i,t) = alpha*x_zeroAccess_lin(i,t-1) + curr_w + u_zeroAccess_lin(i,t-1);
            J_zeroAccess_lin(i,t) = J_zeroAccess_lin(i,t-1) + (Q*(x_zeroAccess_lin(i,t))^2)/T;
            break
            
            
        else
            % evolution
            curr_w = sqrt(W)*randn;
            x_zeroAccess(i,t) = alpha*x_zeroAccess(i,t-1) + curr_w + u_zeroAccess(i,t-1);
            x_zeroAccess_lin(i,t) = alpha*x_zeroAccess_lin(i,t-1) + curr_w + u_zeroAccess_lin(i,t-1);
            x_fullAccess(i,t) = alpha*x_fullAccess(i,t-1) + curr_w + u_fullAccess(i,t-1);
        end
        
        %% observer
        if t==1
            
            % in the first stage, we just normalize the current instant and
            % send to the controller
            a_fullAccess = sqrt(P/W)*x_fullAccess(i,t);
            a_zeroAccess = sqrt(P/W)*x_zeroAccess(i,t);
            a_zeroAccess_lin = sqrt(P/W)*x_zeroAccess_lin(i,t);
        else
            % full access : subtract the last stage receive side
            % estimation, and normalize
            a_fullAccess = sqrt(P/P_error_predict(i,t))*(x_fullAccess(i,t) - x_hat_predict_fullAccess(i,t));
            
            % zero access : subtract the current state, add dither and take
            % modulo 
            curr_dither = 2*delta*(rand - 0.5);
            a_zeroAccess = mod(real(gamma_vec(t))*(x_zeroAccess(i,t) - alpha*x_zeroAccess(i,t-1)) + curr_dither + delta,2*delta) - delta;
            
            % zero access linear subtractor
            norm = (k_vec(t-1)/alpha)^2 * P_error_predict_zeroAccess_lin(i,t) + (1-(k_vec(t-1)/alpha)^2)*W;
            a_zeroAccess_lin = sqrt(P/norm) * (x_zeroAccess_lin(i,t) - (alpha - k_vec(t-1))*x_zeroAccess_lin(i,t-1));
        end
        
        %% channel
        z = sigma_z*randn;
        b_fullAccess = a_fullAccess + z;
        b_zeroAccess = a_zeroAccess + z;
        b_zeroAccess_lin = a_zeroAccess_lin + z;
        
        %% controller
        
        % estimation
        if t==1
            x_hat_zeroAccess(i,t) = snrLin(i)/(1+snrLin(i))*sqrt(W/P)*b_zeroAccess;
            x_hat_estim_zeroAccess_lin(i,t) = snrLin(i)/(1+snrLin(i))*sqrt(W/P)*b_zeroAccess_lin;
            x_hat_estim_fullAccess(i,t) = snrLin(i)/(1+snrLin(i)) * sqrt(W/P) * b_fullAccess;
        else
            
            % full access : update current state estimation 
            x_hat_estim_fullAccess(i,t) = x_hat_predict_fullAccess(i,t) + ...
                snrLin(i)/(1+snrLin(i)) * sqrt(P_error_predict(i,t)/P) * b_fullAccess;
            P_error_estim(i,t) = P_error_predict(i,t)/(1+snrLin(i));
            
            % zero access : modulo receiver and update W estimation
            if abs(gamma_vec(t)*curr_w + z) > delta
                cnt = cnt + 1;
            end
            y_tilde = mod(b_zeroAccess - real(gamma_vec(t))*u_zeroAccess(i,t-1) - curr_dither + delta,2*delta) - delta;
            x_hat_zeroAccess(i,t) = alpha*x_hat_zeroAccess(i,t-1) + u_zeroAccess(i,t-1) + (real(gamma_vec(t))*W/((real(gamma_vec(t))*sqrt(W))^2 + sigma_z^2))*y_tilde;
            
            % zero access linear: recover a CUBE and update estimate
            zeta = (sqrt(P/norm) / (P + sigma_z^2))*((k_vec(t-1)/alpha) * P_error_predict_zeroAccess_lin(i,t) + (1 - k_vec(t-1)/alpha)*W);
            y_tilde = zeta*b_zeroAccess_lin;
            x_hat_estim_zeroAccess_lin(i,t) = x_hat_predict_zeroAccess_lin(i,t) + y_tilde ;
            P_error_estim_zeroAccess_lin(i,t) = (1 - sqrt(P/norm)*zeta*k_vec(t-1)/alpha)*P_error_predict_zeroAccess_lin(i,t) + ...
            zeta^2 * (sigma_z)^2 + ...
            ((zeta*sqrt(P/norm)*(1-k_vec(t-1)/alpha))^2 - 2*zeta*sqrt(P/norm)*(1-k_vec(t-1)/alpha)*(1-k_vec(t-1)*zeta*sqrt(P/norm)/alpha))*W;
        end
        
        % control generation
        u_zeroAccess(i,t) = -k_vec(t)*x_hat_zeroAccess(i,t);
        u_zeroAccess_lin(i,t) = -k_vec(t)*x_hat_estim_zeroAccess_lin(i,t);
        u_fullAccess(i,t) = -k_vec(t)*x_hat_estim_fullAccess(i,t);
        
        % update estimates
        x_hat_predict_fullAccess(i,t+1) = alpha*x_hat_estim_fullAccess(i,t) + u_fullAccess(i,t);
        P_error_predict(i,t+1) = alpha^2/(1+snrLin(i)) * P_error_predict(i,t) + W;
        
        x_hat_predict_zeroAccess_lin(i,t+1) = alpha*x_hat_estim_zeroAccess_lin(i,t) + u_zeroAccess_lin(i,t);
        P_error_predict_zeroAccess_lin(i,t+1) = alpha^2*P_error_estim_zeroAccess_lin(i,t) + W;    
        
        %% Calculate current cost
        if t==1
            J_fullAccess(i,t) = (Q*(x_fullAccess(i,t))^2 + R*(u_fullAccess(i,t))^2)/T;
            J_zeroAccess(i,t) = (Q*(x_zeroAccess(i,t))^2 + R*(u_zeroAccess(i,t))^2)/T;
            J_zeroAccess_lin(i,t) = (Q*(x_zeroAccess_lin(i,t))^2 + R*(u_zeroAccess_lin(i,t))^2)/T;
        else
            J_fullAccess(i,t) = J_fullAccess(i,t-1) + (Q*(x_fullAccess(i,t))^2 + R*(u_fullAccess(i,t))^2)/T;
            J_zeroAccess(i,t) = J_zeroAccess(i,t-1) + (Q*(x_zeroAccess(i,t))^2 + R*(u_zeroAccess(i,t))^2)/T;
            J_zeroAccess_lin(i,t) = J_zeroAccess_lin(i,t-1) + (Q*(x_zeroAccess_lin(i,t))^2 + R*(u_zeroAccess_lin(i,t))^2)/T;
        end
    end
    
    figure;
    % generate sample path plots, and cost
    subplot(121); hold all
    plot(1:T,x_fullAccess(i,:),'-g','LineWidth',2)
    plot(1:T,x_zeroAccess_lin(i,:),'-k','LineWidth',2)
    grid on; grid minor;
    xlabel('t'); ylabel('state')
    legend('x(t) for full access scenarion','x(t) for zero access, linear scheme')
    title(strcat('state , snr = ',num2str(SNR(i)),'[dB]'))
    
    subplot(122); hold all
    plot(1:T,10*log10(J_fullAccess(i,:)),'-gs','LineWidth',2)
    plot(1:T,10*log10(J_zeroAccess_lin(i,:)),'-kp','LineWidth',2)
    grid on; grid minor;
    xlabel('t'); ylabel('cost [dB]');
    title(strcat('Cost in [dB] over time , snr = ',num2str(SNR(i)),'[dB]'));
    legend('Full Access','Zero Access, linear scheme');
    
end

function [gamma_vec,k_vec] = calcLQG(sigma_z,Q,R,W,T,delta,alpha,p_alias)

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

% calculate gamma
gamma_vec = zeros(1,T);
p_alias_overall = p_alias/T;
gamma_vec(1:end) = sqrt((delta/qfuncinv(p_alias_overall/2))^2 - (sigma_z)^2)/sqrt(W);

end
