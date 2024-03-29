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

close all; clear all; clc;

% Init parameters and arrays structures
Params = initParams_Control();
[Results,Debug,PredictErr,EstimErr] = initStructs_Control(Params);

for i=1:length(Params.SNR)
    sigma_z = sqrt(Params.P/Params.snrLin(i));
    sigma_v = sqrt(Params.P/Params.snrFeedback_Lin(i));
    
    xVecEstim = zeros(Params.NumOfSchemes,1); xVecPredict = zeros(Params.NumOfSchemes,1);
    xVec = zeros(Params.NumOfSchemes,1);
    aVec = zeros(Params.NumOfSchemes,1);bVec = zeros(Params.NumOfSchemes,1);
    uVec = zeros(Params.NumOfSchemes,1);
    
    alphaMMSE = Params.snrLin(i)/(1+Params.snrLin(i));
    
    for n = 1 : Params.N_avg
        for t=1:Params.T
            
            % The input driving the schemes is the same input
            curr_w = sqrt(Params.W)*randn;
            %% plant
            if t==1
                % initialization - x_{0} = w_{0}
                xVec = curr_w*ones(Params.NumOfSchemes,1);
            elseif t==Params.T
                % no control in the last stage - evolve state and break
                xVec = Params.alpha*xVec + curr_w*ones(Params.NumOfSchemes,1) + uVec;
                Results(:,i,t,n) = (Results(:,i,t-1,n)*(Params.T-1) + Params.Q*xVec.^2) / Params.T;
                break
            else
                % evolution
                xVec = Params.alpha*xVec + curr_w*ones(Params.NumOfSchemes,1) + uVec;
            end
            
            %% observer
            if t==1
                % in the first stage, we just normalize the current instant and
                % send to the controller
                aVec = sqrt(Params.P/Params.W) * xVec;
                Debug(:,i,t,n) = Params.W;
                
                P_X_linear = Params.P; P_X_modulo = Params.P; P_X_Tuncel = Params.P;
                P_X_modulo_feedback = Params.P; P_X_modulo_feedback_linear = Params.P;
            else
                
                %% full access :
                % subtract the last stage receive side estimation, and normalize
                aVec(1) = sqrt(Params.P/PredictErr(1,i,t,n))*(xVec(1) - xVecPredict(1));
                
                %% Zero access , linear scheme : normalize according to some strategy, rest is same as Kochman-Zamir
                curr_k_for_power = Params.k_for_power(t);
                P_X_linear = Params.W + P_X_linear * (Params.alpha - curr_k_for_power)^2 + ...
                    (curr_k_for_power^2 + 2*curr_k_for_power*(Params.alpha - curr_k_for_power)) * EstimErr(2,i,t-1,n);
                aVec(2) = sqrt(Params.P/P_X_linear) * xVec(2);
                
                %% zero access with modulo : scale the unknown source part using "zooming" factor, and dither and send through the channel
                
                P_X_modulo = Params.W + P_X_modulo * (Params.alpha - curr_k_for_power)^2 + ...
                    (curr_k_for_power^2 + 2*curr_k_for_power*(Params.alpha - curr_k_for_power)) * EstimErr(3,i,t-1,n);
                curr_zoom = Params.gammaKochman*sqrt(Params.P/P_X_modulo);
                
                % Time share between modulo encoded signals with uncoded
                % signals
                if mod(t,Params.uncoded_timesharing_factor) == 0
                    aVec(3) = sqrt(Params.P/P_X_modulo)*xVec(3);
                else
                    curr_dither = 0;
                    aVec(3) = myModulo(curr_zoom*xVec(3) + curr_dither,Params.deltaKochman);
                end
                
                
                %% Tuncel scheme : quantize and send the quantized value + quantization value after proper scaling
                % Tuncel coding
                P_X_Tuncel_Prev = P_X_Tuncel;
                P_X_Tuncel = Params.W + P_X_Tuncel * (Params.alpha - curr_k_for_power)^2 + ...
                    (curr_k_for_power^2 + 2*curr_k_for_power*(Params.alpha - curr_k_for_power)) * EstimErr(4,i,t-1,n);
                
                normX_Tuncel = sqrt(Params.P/P_X_Tuncel) * xVec(4);
                currDistance = abs(normX_Tuncel - Params.centersTuncel(:));
                [~,minIdx] = min(currDistance,[],1);
                T_Tuncel = Params.codebookTuncel(minIdx);
                S_Tuncel = normX_Tuncel - T_Tuncel;
                
                aVec(4) = Params.alphaTuncel*T_Tuncel + Params.betaTuncel*S_Tuncel;
                
                %% partial access with feedback : estimate the noisy control after feedback iterations,then subtract noisy control estimate
                %  feedback is with modulo
                P_X_modulo_feedback = Params.W + P_X_modulo_feedback * (Params.alpha - Params.k_for_power(t))^2 + ...
                    (Params.k_for_power(t)^2 + 2*Params.k_for_power(t)*(Params.alpha - Params.k_for_power(t))) * EstimErr(5,i,t-1,n);
                
                if Params.sign_cheat
                    [aVec(5),MSE_final] = feedback_over_awgn(abs(xVec(5)),Params.P,(1-2/pi)*P_X_modulo_feedback,...
                        Params.SNR(i),deltaSNR,N_feedback,1e-2,'modulo');
                else
                    [aVec(5),MSE_final] = feedback_over_awgn(xVec(5),Params.P,P_X_modulo_feedback,...
                        Params.SNR(i),Params.deltaSNR,Params.N_feedback,1e-2,'modulo');
                end
                
                %% partial access with feedback : estimate the noisy control after feedback iterations,then subtract noisy control estimate
                %  feedback is with linear scaling
                P_X_modulo_feedback_linear = Params.W + P_X_modulo_feedback_linear * (Params.alpha - Params.k_for_power(t))^2 + ...
                    (Params.k_for_power(t)^2 + 2*Params.k_for_power(t)*(Params.alpha - Params.k_for_power(t))) * EstimErr(6,i,t-1,n);
                
                if Params.sign_cheat
                    [aVec(6),MSE_final] = feedback_over_awgn(abs(xVec(6)),Params.P,(1-2/pi)*P_X_modulo_feedback_linear,...
                        Params.SNR(i),deltaSNR,N_feedback,1e-6,'linear');
                else
                    [aVec(6),MSE_final] = feedback_over_awgn(xVec(6),Params.P,P_X_modulo_feedback_linear,...
                        Params.SNR(i),Params.deltaSNR,Params.N_feedback,1e-6,'linear');
                end
                
                %% reversed feedback : The decoder send its last prediction, and feedbacl iterations improves the estimate
                [aVec(7),MSE_final_reverse] = feedback_over_awgn(xVecPredict(7),Params.P,...
                    P_X_modulo_feedback + PredictErr(7,i,t,n),Params.SNR(i) + Params.deltaSNR,...
                    Params.SNR(i),Params.N_feedback,1e-6,'modulo');
                aVec(7) = sqrt(Params.P/(MSE_final_reverse + PredictErr(7,i,t,n)))*(xVec(7) - aVec(7));
                
                %% reversed feedback : The decoder send its last prediction, and feedbacl iterations improves the estimate
                [aVec(8),MSE_final_reverse_linear] = feedback_over_awgn(xVecPredict(8),Params.P,...
                    P_X_modulo_feedback + PredictErr(8,i,t,n),Params.SNR(i) + Params.deltaSNR,...
                    Params.SNR(i),Params.N_feedback,1e-6,'linear');
                aVec(8) = sqrt(Params.P/(MSE_final_reverse + PredictErr(8,i,t,n)))*(xVec(8) - aVec(8));
                
            end
            Debug(:,i,t,n) = aVec(:).^2;
            
            %% channel
            z = sigma_z*randn;
            bVec = aVec + z*ones(Params.NumOfSchemes,1);
            
            if t>1
                bVec(5:6) = aVec(5:6);
            end
            %% controller
            
            % estimation
            if t==1
                xVecEstim = alphaMMSE*sqrt(Params.W/Params.P)*bVec;
                EstimErr(:,i,t,n) = (Params.W/(1+Params.snrLin(i)))*ones(Params.NumOfSchemes,1);
            else
                
                %% full access : update current state estimation
                xVecEstim(1) = xVecPredict(1) + ...
                    alphaMMSE * sqrt(PredictErr(1,i,t,n)/Params.P) * bVec(1);
                EstimErr(1,i,t,n) = PredictErr(1,i,t,n)/(1+Params.snrLin(i));
                
                %% zero access linear: recover a CUBE and update estimate

                zeta = sqrt(Params.P/P_X_linear);
                xi = zeta*PredictErr(2,i,t,n) / (zeta^2 * PredictErr(2,i,t,n) + sigma_z^2);
                y_tilde = xi*(bVec(2) - zeta * (Params.alpha - Params.k_vec(t-1)) * xVecEstim(2));
                xVecEstim(2) = xVecPredict(2) + y_tilde ;
                EstimErr(2,i,t,n) = (1 - zeta*xi)^2*PredictErr(2,i,t,n) + ...
                    xi^2 * (sigma_z)^2;

                %% zero access : modulo receiver
                
                % handle uncoded transmissions
                if mod(t,Params.uncoded_timesharing_factor) == 0
                    
                    zeta_modulo = sqrt(Params.P/P_X_modulo);
                    xi_modulo = zeta_modulo *PredictErr(3,i,t,n) / (zeta_modulo^2 * PredictErr(3,i,t,n) + sigma_z^2);
                    y_tilde_modulo = xi_modulo*(bVec(3) - zeta_modulo * (Params.alpha - Params.k_vec(t-1)) * xVecEstim(3));
                    xVecEstim(3) = xVecPredict(3) + y_tilde_modulo ;
                    EstimErr(3,i,t,n) = (1 - zeta_modulo*xi_modulo)^2*PredictErr(3,i,t,n) + ...
                        xi_modulo^2 * (sigma_z)^2;
                    
                else
                    zeta_modulo = curr_zoom;
                    b_tilde = bVec(3) - zeta_modulo*(Params.alpha - Params.k_vec(t-1))*xVecEstim(3);
                    b_tilde_mod = myModulo(b_tilde,Params.deltaKochman);
                    
                    r_estim = zeta_modulo*PredictErr(3,i,t,n)/(zeta_modulo^2*PredictErr(3,i,t,n) + sigma_z^2);
                    e_hat = r_estim * b_tilde_mod;
                    
                    % update state estimate
                    xVecEstim(3) = xVecPredict(3) + e_hat;
                    
                    EstimErr(3,i,t,n) = (1 - curr_zoom^2*PredictErr(3,i,t,n)/(curr_zoom^2*PredictErr(3,i,t,n) + sigma_z^2))^2*PredictErr(3,i,t,n) + ...
                        (r_estim)^2 * sigma_z^2;
                end
                
                %% Zero Access : Tuncel Receiver
                
                RxSig_tuncel = sqrt(Params.P/P_X_Tuncel) * bVec(4);
                T_hat = MMSE_decoder_Int(RxSig_tuncel,Params.deltaTuncel,Params.codebookTuncel,Params.alphaTuncel,Params.betaTuncel,...
                    PredictErr(4,i,t,n),xVecPredict(4),Params.P,sigma_z);
                
                S_hat = MMSE_decoder_Frac(RxSig_tuncel,Params.deltaTuncel,Params.alphaTuncel,Params.betaTuncel,...
                    PredictErr(4,i,t,n),xVecPredict(4),Params.P,sigma_z,Params.codebookTuncel);
                
                if t < 16
                    T_hat = T_Tuncel;
                end
                xVecEstim(4) = sqrt(P_X_Tuncel/Params.P) * (T_hat + S_hat);
                
                
                %% partial access with feedback : we assume that the feedback output is perfect, so we just treat it
                %% as in the full access scheme
                
                xi = PredictErr(5,i,t,n) / (PredictErr(5,i,t,n) + MSE_final);
                y_tilde = xi*(bVec(5) -(Params.alpha - Params.k_vec(t-1))*xVecEstim(5));
                xVecEstim(5) = xVecPredict(5) + y_tilde ;
                EstimErr(5,i,t,n) = (1 - xi)^2*PredictErr(5,i,t,n) + ...
                    xi^2 * MSE_final;
                
                %% partial access with feedback : we assume that the feedback output is perfect, so we just treat it
                %% as in the full access scheme
                
                xi = PredictErr(6,i,t,n) / (PredictErr(6,i,t,n) + MSE_final);
                y_tilde = xi*(bVec(6) -(Params.alpha - Params.k_vec(t-1))*xVecEstim(6));
                xVecEstim(6) = xVecPredict(6) + y_tilde ;
                EstimErr(6,i,t,n) = (1 - xi)^2*PredictErr(5,i,t,n) + ...
                    xi^2 * MSE_final;
                
                %% reversed feedback : update current state estimation
                xVecEstim(7) = xVecPredict(7) + ...
                    alphaMMSE * sqrt((MSE_final_reverse + PredictErr(7,i,t,n))/Params.P) * bVec(7);
                EstimErr(7,i,t,n) = (MSE_final_reverse + PredictErr(7,i,t,n))/(1+Params.snrLin(i));

                %% reversed feedback : update current state estimation
                xVecEstim(8) = xVecPredict(8) + ...
                    alphaMMSE * sqrt((MSE_final_reverse_linear + PredictErr(8,i,t,n))/Params.P) * bVec(8);
                EstimErr(8,i,t,n) = (MSE_final_reverse_linear + PredictErr(8,i,t,n))/(1+Params.snrLin(i));
                
            end
            
            %% control generation
            uVec = -Params.k_vec(t) * xVecEstim;

            %% update estimates
            xVecPredict = Params.alpha * xVecEstim + uVec;
            PredictErr(:,i,t+1,n) = Params.alpha^2 * EstimErr(:,i,t,n) + Params.W;

            %% Calculate current cost
            if t == 1
               Results(:,i,t,n) =  Params.Q*xVec.^2 + Params.R*uVec.^2 ;
            else
               Results(:,i,t,n) =  (Results(:,i,t-1,n)*(t-1) + (Params.Q*xVec.^2 + Params.R*uVec.^2))/t;
            end

        end
                
        if mod(n,50) == 0
            display(strcat('n = ',num2str(n)));
        end
        
    end
    
    figure;hold all
    lineStyle = ['--g','-.k','-oc','-*g','-*b','-om'];
    linesToPlot = [1 2 7 8];
    for k = 1:length(linesToPlot)
        currMean = mean(reshape(Results(linesToPlot(k),i,:,:),Params.T,[]),2);
        plot(1:(Params.T-1),10*log10(currMean(1:end-1)),lineStyle(3*(k-1) + 1 : 3*k),'LineWidth',2)
    end
    plot(1:Params.T,10*log10(Params.s_vec(2)*Params.W + ((Params.Q + (Params.alpha^2-1)*Params.s_vec(2))/(1+Params.snrLin(i)-Params.alpha^2))*Params.W)*ones(1,Params.T),'--y','LineWidth',2)
    plot(1:Params.T,10*log10(Params.s_vec(2)*Params.W)*ones(1,Params.T),'--m','LineWidth',2)
    grid on; grid minor;
    xlabel('t'); ylabel('cost [dB]');
    legend('Full Access','Zero Access',...
        strcat('Modulo Reversed feedback : ',num2str(Params.N_feedback),' iterations'),...
        strcat('Linear Reversed feedback : ',num2str(Params.N_feedback),' iterations'),'LQG_{\infty}(SDR_{0})','LQG_{\infty}(SNR = \infty)');
%     legend('Full Access','Zero Access','Zero Access - Modulo Encoder (Kochman)','Zero Access - HDA Encoder (Tuncel)',...
%         strcat('Modulo feedback : ',num2str(Params.N_feedback),' iterations'),...
%         strcat('Linear feedback : ',num2str(Params.N_feedback),' iterations'),...
%         strcat('Reversed feedback : ',num2str(Params.N_feedback),' iterations'),'SDR OPTA');
%         strcat('Linear reversed feedback : ',num2str(Params.N_feedback),' iterations'),'SDR OPTA');
    
    title(strcat('SNR = ',num2str(Params.SNR),'[dB], Feedback SNR = ',num2str(Params.SNR + Params.deltaSNR),'[dB]'))
    powerAnalysis(reshape(Debug(3,i,:,:),Params.T,[]),Params.N_avg,'Kochman',sigma_z)
    powerAnalysis(reshape(Debug(4,i,:,:),Params.T,[]),Params.N_avg,'Tuncel',sigma_z)
    
    
end

function [] = powerAnalysis(P_x,N_avg,run,sigma)

inputPower = sum(P_x,2)/N_avg;

figure;hold all
plot(inputPower)
title([run strcat('measured SNR = ',num2str(10*log10(mean(inputPower(2:end-1))/sigma^2)))])

end

function [T_hat] = MMSE_decoder_Int(RxSig_tuncel,delta,codebook,alpha,beta,PredictErr,PredictState,StatePower,sigmaW)

ds = 0.001;
s = -delta/2:ds:delta/2;

% start with 2D calculations
calc1 = codebook(:) + s;
calc2 = alpha*codebook(:) + beta*s;

expTerm1 = -1*(0.5*(calc1/sqrt(StatePower)).^2);
expTerm2 = - (0.5/(PredictErr^2))*(PredictState - calc1).^2;
expTerm3 = 0.5*(1/sigmaW^2)*(RxSig_tuncel - calc2).^2;

Integrand = exp(expTerm1 + expTerm2 - expTerm3);

num = sum(Integrand,2);
T_hat = sum(codebook(:).*num(:),1)./sum(num(:));

% quantize to the nearset r_k
% [~,minIdx] = min(abs(T_hat - codebook));
% T_hat = codebook(minIdx);
end

function [S_hat] = MMSE_decoder_Frac(RxSig_tuncel,delta,alpha,beta,PredictErr,PredictState,StatePower,sigmaW,codebook)


ds = 0.001;
s = ones(length(codebook),1) * (-delta/2:ds:delta/2);

% start with 2D calculations
calc1 = codebook(:) + s;
calc2 = alpha*codebook(:) + beta*s;

expTerm1 = -1*(0.5*(calc1/sqrt(StatePower)).^2);
expTerm2 = - (0.5/(PredictErr^2))*(PredictState - calc1).^2;
expTerm3 = 0.5*(1/sigmaW^2)*(RxSig_tuncel - calc2).^2;
Integrand = exp(expTerm1 + expTerm2 - expTerm3);


num = sum(sum(Integrand.*s,2),1);
den = sum(sum(Integrand));
S_hat = num/den;

end


