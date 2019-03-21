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
W = 1; % input gaussian variance
T = 120; % horizon
delta = 3;
P = delta^2/3;

% cost parameters
Q = 1;
R = 2;


SNR = [25 20 15 10 8 5];
deltaSNR = 5;
snrLin = 10.^(SNR/10);
snrFeedback_Lin = 10.^((SNR + deltaSNR)/10);

initArrays;

P_alias = 5*1e-2;
N_feedback = 8;
% set timesharing parameters - we need to scale the powers to the overall
% avergae power will be P
uncoded_timesharing_factor = 4;

cnt = 0;
for i=1:length(SNR)
    sigma_z = sqrt(P/10^(SNR(i)/10));
    sigma_v = sqrt(P/10^((SNR(i) + deltaSNR)/10));
    % calculate LQG parameters for current run
    k_vec = calcLQG(Q,R,T,alpha);
    P_X_linear = zeros(1,T);
    P_X_hat_linear = zeros(1,T);
    
    P_X_modulo = zeros(1,T);
    P_X_hat_modulo = zeros(1,T);
    
    P_X_feedback = zeros(1,T);
    P_X_hat_feedback = zeros(1,T);
    for t=1:T
        
        %% plant
        if t==1
            % initialization
            curr_w = sqrt(W)*randn;
            x_fullAccess(i,t) = curr_w;
            x_zeroAccess(i,t) = curr_w;
            x_partialAccess(i,t) = curr_w;
            x_zeroAccess_modulo(i,t) = curr_w;
            x_partialAccess_feedback(i,t) = curr_w;
            
        elseif t==T
            curr_w = sqrt(W)*randn;
            % no control in the last stage - evolve state and break
            x_fullAccess(i,t) = alpha*x_fullAccess(i,t-1) + curr_w + u_fullAccess(i,t-1);
            J_fullAccess(i,t) = J_fullAccess(i,t-1) + (Q*(x_fullAccess(i,t))^2)/T;
            
            % no control in the last stage - evolve state and break
            x_zeroAccess(i,t) = alpha*x_zeroAccess(i,t-1) + curr_w + u_zeroAccess(i,t-1);
            J_zeroAccess(i,t) = J_zeroAccess(i,t-1) + (Q*(x_zeroAccess(i,t))^2)/T;
            
            % no control in the last stage - evolve state and break
            x_partialAccess(i,t) = alpha*x_partialAccess(i,t-1) + curr_w + u_partialAccess(i,t-1);
            J_partialAccess(i,t) = J_partialAccess(i,t-1) + (Q*(x_partialAccess(i,t))^2)/T;
            
            % no control in the last stage - evolve state and break
            x_zeroAccess_modulo(i,t) = alpha*x_zeroAccess_modulo(i,t-1) + curr_w + u_zeroAccess_modulo(i,t-1);
            J_zeroAccess_modulo(i,t) = J_zeroAccess_modulo(i,t-1) + (Q*(x_zeroAccess_modulo(i,t))^2)/T;
            
            % no control in the last stage - evolve state and break
            x_partialAccess_feedback(i,t) = alpha*x_partialAccess_feedback(i,t-1) + curr_w + u_partialAccess_feedback(i,t-1);
            J_partialAccess_feedback(i,t) = J_partialAccess_feedback(i,t-1) + (Q*(x_partialAccess_feedback(i,t))^2)/T;
            
            break
            
            
        else
            % evolution
            curr_w = sqrt(W)*randn;
            x_zeroAccess(i,t) = alpha*x_zeroAccess(i,t-1) + curr_w + u_zeroAccess(i,t-1);
            x_fullAccess(i,t) = alpha*x_fullAccess(i,t-1) + curr_w + u_fullAccess(i,t-1);
            x_partialAccess(i,t) = alpha*x_partialAccess(i,t-1) + curr_w + u_partialAccess(i,t-1);
            x_zeroAccess_modulo(i,t) = alpha*x_zeroAccess_modulo(i,t-1) + curr_w + u_zeroAccess_modulo(i,t-1);
            x_partialAccess_feedback(i,t) = alpha*x_partialAccess_feedback(i,t-1) + curr_w + u_partialAccess_feedback(i,t-1);
            
        end
        
        %% observer
        if t==1
            
            % in the first stage, we just normalize the current instant and
            % send to the controller
            a_fullAccess = sqrt(P/W)*x_fullAccess(i,t);
            a_zeroAccess = sqrt(P/W)*x_zeroAccess(i,t);
            a_partialAccess = sqrt(P/W)*x_partialAccess(i,t);
            a_zeroAccess_modulo = sqrt(P/W)*x_zeroAccess_modulo(i,t);
            a_partialAccess_feedback = sqrt(P/W)*x_partialAccess_feedback(i,t);
            
            P_X_linear(t) = W;
            P_X_modulo(t) = W;
            P_X_feedback(t) = W;
        else
            %% full access :
            % subtract the last stage receive side estimation, and normalize
            a_fullAccess = sqrt(P/P_error_predict(i,t))*(x_fullAccess(i,t) - x_hat_predict_fullAccess(i,t));
            
            %% zero access : linear subtractor
            norm = (k_vec(t-1))^2 * P_error_estim_zeroAccess(i,t-1) + W;
            a_zeroAccess = sqrt(P/norm) * (x_zeroAccess(i,t) - (alpha - k_vec(t-1))*x_zeroAccess(i,t-1));
            
            %% zero access with modulo : scale the unknown source part using "zooming" factor, and dither and send through the channel
            
            if mod(t-1,uncoded_timesharing_factor) == 0
                P_X_modulo(t) = W + (k_vec(t-1))^2 * P_error_estim_zeroAccess_modulo(i,t-1);
                P_X_hat_modulo(t) = P_X_modulo(t) + P_error_predict_zeroAccess_modulo(i,t-1);
                curr_zoom = calc_zoom(P_alias,W,alpha,P,P_X_modulo(t),snrLin(i));
            else
                P_X_modulo(t) = W + P_X_modulo(t-1) * (alpha - k_vec(t-1))^2 + (k_vec(t-1))^2 * P_error_estim_zeroAccess_modulo(i,t-1);
                P_X_hat_modulo(t) = P_X_modulo(t) + P_error_predict_zeroAccess_modulo(i,t-1);
                curr_zoom = calc_zoom(P_alias,W,alpha,P,P_X_modulo(t-1),snrLin(i));
            end
            
            curr_dither = delta*(rand - 0.5);
            
            % combine uncoded transmission sometime
            if mod(t,uncoded_timesharing_factor) == 0
                a_zeroAccess_modulo = x_zeroAccess_modulo(i,t);
            else
                a_zeroAccess_modulo = mod(curr_zoom*x_zeroAccess_modulo(i,t) + curr_dither + delta,2*delta) - delta;
            end
            
            %% partial access : estimate the noisy control, and then subtract noisy control estimate
            P_X_linear(t) = W + P_X_linear(t-1) * (alpha - k_vec(t-1))^2 + (k_vec(t-1))^2 * P_error_estim_partialAccess(i,t-1);
            P_X_hat_linear(t) = P_X_linear(t) + P_error_predict_partialAccess(i,t-1);
            
            r = sqrt(P/P_X_hat_linear(t))*x_hat_predict_partialAccess(i,t) + sigma_v*randn;
            estim_coeff = sqrt(P_X_hat_linear(t)/P);
            x_estim = estim_coeff*r;
            
            norm_partialAccess = P_error_predict_partialAccess(i,t) + estim_coeff^2 * sigma_v^2;
            a_partialAccess = sqrt(P/norm_partialAccess) * (x_partialAccess(i,t) - x_estim);
            
            %% partial access with feedback : estimate the noisy control after feedback iterations,then subtract noisy control estimate
            P_X_feedback(t) = W + P_X_feedback(t-1) * (alpha - k_vec(t-1))^2 + (k_vec(t-1))^2 * P_error_estim_partialAccess_feedback(i,t-1);
            P_X_hat_feedback(t) = P_X_feedback(t) + P_error_predict_partialAccess_feedback(i,t-1);
            
            [x_final,MSE_final] = feedback_over_awgn(x_hat_predict_partialAccess_feedback(i,t),P,P_X_hat_linear(t),SNR(i) + deltaSNR,...
                -1*deltaSNR,N_feedback,P_alias);
            
            norm_partialAccess_feedback = P_error_predict_partialAccess_feedback(i,t) + MSE_final;
            a_partialAccess_feedback = sqrt(P/norm_partialAccess_feedback) * (x_partialAccess_feedback(i,t) - x_final);
        end
        
        %% channel
        z = sigma_z*randn;
        b_fullAccess = a_fullAccess + z;
        b_zeroAccess = a_zeroAccess + z;
        b_partialAccess = a_partialAccess + z;
        b_zeroAccess_modulo = a_zeroAccess_modulo + z;
        b_partialAccess_feedback = a_partialAccess_feedback + z;
        
        %% controller
        
        % estimation
        if t==1
            x_hat_estim_zeroAccess(i,t) = snrLin(i)/(1+snrLin(i))*sqrt(W/P)*b_zeroAccess;
            x_hat_estim_fullAccess(i,t) = snrLin(i)/(1+snrLin(i)) * sqrt(W/P) * b_fullAccess;
            x_hat_estim_partialAccess(i,t) = snrLin(i)/(1+snrLin(i)) * sqrt(W/P) * b_partialAccess;
            x_hat_estim_zeroAccess_modulo(i,t) = snrLin(i)/(1+snrLin(i)) * sqrt(W/P) * b_zeroAccess_modulo;
            x_hat_estim_partialAccess_feedback(i,t) = snrLin(i)/(1+snrLin(i)) * sqrt(W/P) * b_partialAccess_feedback;
        else
            
            %% full access : update current state estimation
            x_hat_estim_fullAccess(i,t) = x_hat_predict_fullAccess(i,t) + ...
                snrLin(i)/(1+snrLin(i)) * sqrt(P_error_predict(i,t)/P) * b_fullAccess;
            P_error_estim(i,t) = P_error_predict(i,t)/(1+snrLin(i));
            
            %% zero access linear: recover a CUBE and update estimate
            zeta = (sqrt(P/norm) / (P + sigma_z^2))*((k_vec(t-1)/alpha) * P_error_predict_zeroAccess(i,t) + (1 - k_vec(t-1)/alpha)*W);
            y_tilde = zeta*b_zeroAccess;
            x_hat_estim_zeroAccess(i,t) = x_hat_predict_zeroAccess(i,t) + y_tilde ;
            P_error_estim_zeroAccess(i,t) = (1 - sqrt(P/norm)*zeta*k_vec(t-1)/alpha)^2*P_error_predict_zeroAccess(i,t) + ...
                zeta^2 * (sigma_z)^2 + ...
                ((zeta*sqrt(P/norm)*(1-k_vec(t-1)/alpha))^2 - 2*zeta*sqrt(P/norm)*(1-k_vec(t-1)/alpha)*(1-k_vec(t-1)*zeta*sqrt(P/norm)/alpha))*W;
            
            %% partial access : update current state estimation
            gamma = sqrt(P/norm_partialAccess);
            partial_zeta = gamma*P_error_predict_partialAccess(i,t)/(P + (sigma_z)^2);
            x_hat_estim_partialAccess(i,t) = x_hat_predict_partialAccess(i,t) + partial_zeta * b_partialAccess;
%             P_error_estim_partialAccess(i,t) = (1-gamma*partial_zeta)^2*P_error_predict_partialAccess(i,t) + ...
%                 (partial_zeta*gamma)^2 * (sigma_v)^2 + zeta^2*(sigma_z)^2;
            P_error_estim_partialAccess(i,t) = (1-gamma*partial_zeta)^2*P_error_predict_partialAccess(i,t) + ...
                (partial_zeta*gamma*estim_coeff)^2 * (sigma_v)^2 + zeta^2*(sigma_z)^2;                
            %% zero access : modulo receiver
            
            % handle uncoded transmissions
            if mod(t,uncoded_timesharing_factor) == 0
                b_tilde = b_zeroAccess_modulo;
                
                % update state estimate
                x_hat_estim_zeroAccess_modulo(i,t) = b_tilde;
                P_error_estim_zeroAccess_modulo(i,t) = sigma_z^2;
                
            elseif mod(t-1,uncoded_timesharing_factor) == 0
                alpha_mmse = snrLin(i)/(1 + snrLin(i));
                b_tilde = mod(delta + alpha_mmse*b_zeroAccess_modulo - curr_dither,2*delta) - delta;
                
                % update state estimate
                x_hat_estim_zeroAccess_modulo(i,t) = b_tilde/curr_zoom;
                
                P_error_estim_zeroAccess_modulo(i,t) = (P/(1+snrLin(i)))/curr_zoom^2;
            else
                alpha_mmse = snrLin(i)/(1 + snrLin(i));
                b_tilde = mod(delta + alpha_mmse*b_zeroAccess_modulo - curr_dither - curr_zoom*u_zeroAccess_modulo(i,t-1),2*delta) - delta;
                
                % update state estimate
                x_hat_estim_zeroAccess_modulo(i,t) = b_tilde/curr_zoom + u_zeroAccess_modulo(i,t-1);
                
                P_error_estim_zeroAccess_modulo(i,t) = (P/(1+snrLin(i)))/curr_zoom^2;
            end
            
            %% partial access with feedback : we assume that the feedback output is perfect, so we just treat it
            %% as in the full access scheme
            x_hat_estim_partialAccess_feedback(i,t) = x_hat_predict_partialAccess_feedback(i,t) + ...
                snrLin(i)/(1+snrLin(i)) * sqrt(norm_partialAccess_feedback/P) * b_partialAccess_feedback;
            P_error_estim_partialAccess_feedback(i,t) = P_error_predict_partialAccess_feedback(i,t)/(1+snrLin(i));
            
        end
        
        %% control generation
        u_zeroAccess(i,t) = -k_vec(t)*x_hat_estim_zeroAccess(i,t);
        u_partialAccess(i,t) = -k_vec(t)*x_hat_estim_partialAccess(i,t);
        u_fullAccess(i,t) = -k_vec(t)*x_hat_estim_fullAccess(i,t);
        u_partialAccess_feedback(i,t) = -k_vec(t)*x_hat_estim_partialAccess_feedback(i,t);
        
        if mod(t,uncoded_timesharing_factor) == 0
            u_zeroAccess_modulo(i,t) = -alpha*x_hat_estim_zeroAccess_modulo(i,t);
        else
            u_zeroAccess_modulo(i,t) = -k_vec(t)*x_hat_estim_zeroAccess_modulo(i,t);
        end
        %% update estimates
        x_hat_predict_fullAccess(i,t+1) = alpha*x_hat_estim_fullAccess(i,t) + u_fullAccess(i,t);
        P_error_predict(i,t+1) = alpha^2 * P_error_estim(i,t) + W;
        
        x_hat_predict_zeroAccess(i,t+1) = alpha*x_hat_estim_zeroAccess(i,t) + u_zeroAccess(i,t);
        P_error_predict_zeroAccess(i,t+1) = alpha^2*P_error_estim_zeroAccess(i,t) + W;
        
        x_hat_predict_partialAccess(i,t+1) = alpha*x_hat_estim_partialAccess(i,t) + u_partialAccess(i,t);
        P_error_predict_partialAccess(i,t+1) = alpha^2*P_error_estim_partialAccess(i,t) + W;
        
        x_hat_predict_zeroAccess_modulo(i,t+1) = alpha*x_hat_estim_zeroAccess_modulo(i,t) + u_zeroAccess_modulo(i,t);
        P_error_predict_zeroAccess_modulo(i,t+1) = alpha^2*P_error_estim_zeroAccess_modulo(i,t) + W;
        
        x_hat_predict_partialAccess_feedback(i,t+1) = alpha*x_hat_estim_partialAccess_feedback(i,t) + u_partialAccess_feedback(i,t);
        P_error_predict_partialAccess_feedback(i,t+1) = alpha^2*P_error_estim_partialAccess_feedback(i,t) + W;
        
        %% Calculate current cost
        if t==1
            J_fullAccess(i,t) = (Q*(x_fullAccess(i,t))^2 + R*(u_fullAccess(i,t))^2)/T;
            J_zeroAccess(i,t) = (Q*(x_zeroAccess(i,t))^2 + R*(u_zeroAccess(i,t))^2)/T;
            J_partialAccess(i,t) = (Q*(x_partialAccess(i,t))^2 + R*(u_partialAccess(i,t))^2)/T;
            J_zeroAccess_modulo(i,t) = (Q*(x_zeroAccess_modulo(i,t))^2 + R*(u_zeroAccess_modulo(i,t))^2)/T;
            J_partialAccess_feedback(i,t) = (Q*(x_partialAccess_feedback(i,t))^2 + R*(u_partialAccess_feedback(i,t))^2)/T;
        else
            J_fullAccess(i,t) = J_fullAccess(i,t-1) + (Q*(x_fullAccess(i,t))^2 + R*(u_fullAccess(i,t))^2)/T;
            J_zeroAccess(i,t) = J_zeroAccess(i,t-1) + (Q*(x_zeroAccess(i,t))^2 + R*(u_zeroAccess(i,t))^2)/T;
            J_partialAccess(i,t) = J_partialAccess(i,t-1) + (Q*(x_partialAccess(i,t))^2 + R*(u_partialAccess(i,t))^2)/T;
            J_zeroAccess_modulo(i,t) = J_zeroAccess_modulo(i,t-1) + (Q*(x_zeroAccess_modulo(i,t))^2 + R*(u_zeroAccess_modulo(i,t))^2)/T;
            J_partialAccess_feedback(i,t) = J_partialAccess_feedback(i,t-1) + (Q*(x_partialAccess_feedback(i,t))^2 + R*(u_partialAccess_feedback(i,t))^2)/T;
        end
    end
    
    figure;
    % generate sample path plots, and cost
    hold all
    plot(1:T,x_fullAccess(i,:),'-g','LineWidth',1.5)
    plot(1:T,x_zeroAccess(i,:),'-k','LineWidth',1.5)
    plot(1:T,x_partialAccess(i,:),'-r','LineWidth',1.5)
    plot(1:T,x_zeroAccess_modulo(i,:),'--c','LineWidth',1.5)
    plot(1:T,x_partialAccess_feedback(i,:),'-b','LineWidth',1.5)
    
    xlabel('t'); ylabel('state')
    legend('x(t) for full control access','x(t) for zero control access','x(t) for partial control access',...
        'x(t) for zero control access, with Kochman-Zamir',...
        strcat('x(t) for partial control access, with feedback :',num2str(N_feedback),' iterations'))
    title(strcat('state , snr = ',num2str(SNR(i)),'[dB]',' \Deltasnr = ',num2str(deltaSNR)))
    grid on; grid minor;
    
    figure; hold all
    plot(1:T,10*log10(J_fullAccess(i,:)),'-gs','LineWidth',1.5)
    plot(1:T,10*log10(J_zeroAccess(i,:)),'-kp','LineWidth',1.5)
    plot(1:T,10*log10(J_partialAccess(i,:)),'-r*','LineWidth',1.5)
    plot(1:T,10*log10(J_zeroAccess_modulo(i,:)),'-cd','LineWidth',1.5)
    plot(1:T,10*log10(J_partialAccess_feedback(i,:)),'-b>','LineWidth',1.5)
    grid on; grid minor;
    xlabel('t'); ylabel('cost [dB]');
    title(strcat('Cost in [dB] over time , snr = ',num2str(SNR(i)),'[dB]',' \Deltasnr = ',num2str(deltaSNR)));
    legend('Full Access','Zero Access','Partial Access',...
        'Zero Access - kochman zamir',strcat('Partial Access - feedback : ',num2str(N_feedback),' iterations'));
    
end

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

function [curr_zoom] = calc_zoom(P_alias,W,alpha,P,P_X,snrLin)

a = 3*P / ((qfuncinv(P_alias/2))^2);
b = P/(1+snrLin);
c = alpha^2 * P_X + W;

curr_zoom = sqrt((a - b)/c);

end