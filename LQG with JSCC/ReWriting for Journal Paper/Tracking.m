%% Simulation of Markov Source tracking over power limited awgn channel
% we simulate 3 cases :
% - sensor and observer has full access to the last state estimation \hat{s}_t|t
% - only sensor has the estimation
% - noisy feedback between both

% The evolution of s_t is given by :
% s_(t+1) = lambda*s_t + w_(t+1)
% where w_(t+1) is wgn with variance W

% The goal is to minimize the estimation error ||s_{t} - \hat{s}_{t|t}||^2
clear; clc; close all

% Init parameters and arrays structures
Params = initParams_Tracking();
[Results,Debug,PredictErr,EstimErr] = initStructs(Params);

% Actual simulation starts here
for i=1:length(Params.SNR)
    sigma_z = sqrt(Params.P/Params.snrLin(i));
    sigma_v = sqrt(Params.P/Params.snrFeedback_Lin(i));
    
    xVec = zeros(Params.NumOfSchemes,1);xVecPredict = zeros(Params.NumOfSchemes,1);
    xVecEstim = zeros(Params.NumOfSchemes,1);
    aVec = zeros(Params.NumOfSchemes,1);bVec = zeros(Params.NumOfSchemes,1);
    
    alphaMMSE = Params.snrLin(i)/(1+Params.snrLin(i));
    for n = 1:Params.N_avg
        for t=1:Params.T
            
            % The input driving the schemes is the same input
            curr_w = sqrt(Params.W)*randn;
            %% plant
            if t==1
                % initialization - x_{0} = w_{0}
                xVec = curr_w*ones(Params.NumOfSchemes,1);
            elseif t==Params.T
                % no control in the last stage - evolve state and break
                xVec = Params.lambda*xVec + curr_w*ones(Params.NumOfSchemes,1);
                Results(:,i,t,n) = (xVec - xVecPredict).^2;
                break
            else
                % evolution
                xVec = Params.lambda*xVec + curr_w*ones(Params.NumOfSchemes,1);
            end
            
            %% observer
            if t==1
                % in the first stage, we just normalize the current instant and
                % send to the controller
                aVec = sqrt(Params.P/Params.W) * xVec;
                Debug(:,i,t,n) = Params.W;
            else
                %% full access :
                % subtract the last stage receive side estimation, and normalize
                aVec(1) = sqrt(Params.P/PredictErr(1,i,t,n))*(xVec(1) - xVecPredict(1));
                
                %% Zero access , linear scheme :
                aVec(2) = sqrt(Params.P/Params.statePower(t)) * xVec(2);
                
                %% zero access with modulo : scale the unknown source part using "zooming" factor, and dither and send through the channel
                
                curr_zoom = Params.gammaKochman * sqrt(Params.P/Params.statePower(t));
                % calc_zoom(Params.P_alias,Params.P,PredictErr(3,i,t,n),Params.snrLin(i));
                curr_dither = 0;
                
                % combine uncoded transmission sometime
                if mod(t,Params.uncoded_timesharing_factor) == 0
                    aVec(3) = sqrt(Params.P_timesharing/Params.statePower(t))*xVec(3);
                else
%                     aVec(3) = myModulo(curr_zoom*xVec(3) + curr_dither,Params.deltaKochman);
                        aVec(3) = myModulo(curr_zoom*xVec(3) + curr_dither,Params.deltaKochman);
                end
                
                %% Tuncel scheme : quantize and send the quantized value + quantization value after proper scaling
                % Tuncel coding
                normX_Tuncel = sqrt(Params.P/Params.statePower(t)) * xVec(4);
                currDistance = abs(normX_Tuncel - Params.centersTuncel(:));
                [~,minIdx] = min(currDistance,[],1);
                T_Tuncel = Params.codebookTuncel(minIdx);
                S_Tuncel = normX_Tuncel - T_Tuncel;
                
                aVec(4) = Params.alpha*T_Tuncel + Params.beta*S_Tuncel;
                
                %% partial access with feedback : estimate the noisy control after feedback iterations,then subtract noisy control estimate
                
                [aVec(5),MSE_final] = feedback_over_awgn(xVec(5),Params.P,Params.statePower(t),Params.SNR(i),...
                    Params.deltaSNR,Params.N_feedback,1e-3,'modulo');
                
                %% partial access with feedback : estimate the noisy control after feedback iterations,then subtract noisy control estimate
                %  feedback is with linear scaling
                [aVec(6),MSE_final_linear] = feedback_over_awgn(xVec(6),Params.P,Params.statePower(t),Params.SNR(i),...
                    Params.deltaSNR,Params.N_feedback,1e-3,'linear');

                %% reversed feedback : The decoder send its last prediction, and feedbacl iterations improves the estimate 
                [aVec(7),MSE_final_reverse] = feedback_over_awgn(xVecPredict(7),Params.P,...
                    Params.statePower(t) + PredictErr(7,i,t,n),Params.SNR(i) + Params.deltaSNR,...
                    Params.SNR(i),Params.N_feedback,1e-5,'modulo');
                aVec(7) = sqrt(Params.P/(MSE_final_reverse + PredictErr(7,i,t,n)))*(xVec(7) - aVec(7));
                
                %% reversed feedback : The decoder send its last prediction, and feedbacl iterations improves the estimate 
                [aVec(8),MSE_final_reverse_linear] = feedback_over_awgn(xVecPredict(8),Params.P,...
                    Params.statePower(t) + PredictErr(8,i,t,n),Params.SNR(i) + Params.deltaSNR,...
                    Params.SNR(i),Params.N_feedback,1e-5,'modulo');
                aVec(8) = sqrt(Params.P/(MSE_final_reverse + PredictErr(8,i,t,n)))*(xVec(8) - aVec(8));
            end
            Debug(:,i,t,n) = aVec(:).^2;
            
            %% channel
            z = sigma_z*randn;
            bVec = aVec + z*ones(Params.NumOfSchemes,1);
            
            if t>1
                bVec(5:6) = aVec(5:6);
            end
            %% sensor
            
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
                
                zeta = sqrt(Params.P/Params.statePower(t));
                xi = zeta*PredictErr(2,i,t,n) / (zeta^2 * PredictErr(2,i,t,n) + sigma_z^2);
                y_tilde = xi*(bVec(2) - zeta*Params.lambda * xVecEstim(2));
                xVecEstim(2) = xVecPredict(2) + y_tilde ;
                EstimErr(2,i,t,n) = (1 - zeta*xi)^2*PredictErr(2,i,t,n) + ...
                    xi^2 * (sigma_z)^2;
                
                %% zero access : modulo receiver
                
                % handle uncoded transmissions
                if mod(t,Params.uncoded_timesharing_factor) == 0
                    
                    zeta_modulo = sqrt(Params.P_timesharing/Params.statePower(t));
                    xi_modulo = zeta_modulo *PredictErr(3,i,t,n) / (zeta_modulo^2 * PredictErr(3,i,t,n) + sigma_z^2);
                    y_tilde_modulo = xi_modulo*(bVec(3) - zeta_modulo * Params.lambda * xVecEstim(3));
                    xVecEstim(3) = xVecPredict(3) + y_tilde_modulo ;
                    EstimErr(3,i,t,n) = (1 - zeta_modulo*xi_modulo)^2*PredictErr(3,i,t,n) + ...
                        xi_modulo^2 * (sigma_z)^2;
                    
                else
                    
                    b_tilde = bVec(3) - curr_dither - curr_zoom*xVecPredict(3);
                    b_tilde_mod = myModulo(b_tilde,Params.deltaKochman);
                    
                    r_estim = curr_zoom*PredictErr(3,i,t,n)/(curr_zoom^2*PredictErr(3,i,t,n) + sigma_z^2);
                    e_hat = r_estim * b_tilde_mod;
                    
                    % update state estimate
                    xVecEstim(3) = xVecPredict(3) + e_hat;
                    
                    EstimErr(3,i,t,n) = (1 - curr_zoom^2*PredictErr(3,i,t,n)/(curr_zoom^2*PredictErr(3,i,t,n) + sigma_z^2))^2*PredictErr(3,i,t,n) + ...
                        (r_estim)^2 * sigma_z^2;
                    
                end
                
                %% Zero Access : Tuncel Receiver
                
                % start with parameter estimation
                % rho = Params.lambda^2*Params.statePower(t) / (Params.lambda^2*Params.statePower(t) + Params.W);
                rho = Params.lambda ^ 2 ;
                y = sqrt(Params.P/Params.statePower(t))*xVecPredict(4);
                
                T_hat = MMSE_decoder_Int(bVec(4),y,...
                    Params.deltaTuncel,Params.codebookTuncel,Params.alpha,Params.beta,rho,sigma_z);
                S_hat = MMSE_decoder_Frac(bVec(4),y,...
                    Params.deltaTuncel,Params.alpha,Params.beta,rho,sigma_z,Params.codebookTuncel);
                
                if t < 16
                    T_hat = T_Tuncel;
                end
                xVecEstim(4) = sqrt(Params.statePower(t)/Params.P) * (T_hat + S_hat);
                
                %% partial access with feedback : we assume that the feedback output is perfect, so we just treat it
                %% as in the full access scheme
                
                xi = PredictErr(5,i,t,n) / (PredictErr(5,i,t,n) + MSE_final);
                y_tilde = xi*(bVec(5) - Params.lambda*xVecEstim(5));
                xVecEstim(5) = xVecPredict(5) + y_tilde ;
                EstimErr(5,i,t,n) = (1 - xi)^2*PredictErr(5,i,t,n) + ...
                    xi^2 * MSE_final;
                
                %% partial access with feedback : we assume that the feedback output is perfect, so we just treat it
                %% as in the full access scheme
                
                xi = PredictErr(6,i,t,n) / (PredictErr(6,i,t,n) + MSE_final_linear);
                y_tilde = xi*(bVec(6) - Params.lambda*xVecEstim(6));
                xVecEstim(6) = xVecPredict(6) + y_tilde ;
                EstimErr(6,i,t,n) = (1 - xi)^2*PredictErr(6,i,t,n) + ...
                    xi^2 * MSE_final_linear;
                
                %% reversed feedback : update current state estimation
                xVecEstim(7) = xVecPredict(7) + ...
                    alphaMMSE * sqrt((MSE_final_reverse + PredictErr(7,i,t,n))/Params.P) * bVec(7);
                EstimErr(7,i,t,n) = (MSE_final_reverse + PredictErr(7,i,t,n))/(1+Params.snrLin(i));

                %% reversed feedback : update current state estimation
                xVecEstim(8) = xVecPredict(8) + ...
                    alphaMMSE * sqrt((MSE_final_reverse_linear + PredictErr(8,i,t,n))/Params.P) * bVec(8);
                EstimErr(8,i,t,n) = (MSE_final_reverse_linear + PredictErr(8,i,t,n))/(1+Params.snrLin(i));
                
            end
            
            %% update estimates
            xVecPredict = Params.lambda*xVecEstim;
            PredictErr(:,i,t+1,n) = Params.lambda^2 * EstimErr(:,i,t,n) + Params.W;
            
            %% Calculate current cost
            Results(:,i,t,n) = (xVec - xVecEstim).^2;
            
        end
        
        if mod(n,50) == 0
            display(strcat('n = ',num2str(n)));
        end
    end
    figure;hold all
    lineStyle = ['-g','--k',':c',':g','b',':m'];
    for k=1:(size(Results,1)-1)
        currMean = mean(reshape(Results(k,i,:,:),Params.T,[]),2);
        plot(1:(Params.T-1),10*log10(currMean(1:end-1)),'LineWidth',2)
    end
    plot(1:(Params.T-1),10*log10(Params.statePower(1:end-1) * (1-Params.lambda^2)/Params.snrLin),'--','LineWidth',2)
    grid on; grid minor;
    xlabel('t'); ylabel('cost [dB]');
    legend('Full Access','Zero Access','Zero Access - Modulo Encoder (Kochman)','Zero Access - HDA Encoder (Tuncel)',...
        strcat('Modulo feedback : ',num2str(Params.N_feedback),' iterations'),...
        strcat('Linear feedback : ',num2str(Params.N_feedback),' iterations'),...
        strcat('Reversed feedback : ',num2str(Params.N_feedback),' iterations'),'SDR OPTA');
%         strcat('Linear reversed feedback : ',num2str(Params.N_feedback),' iterations'),'SDR OPTA');
    
    title(strcat('SNR = ',num2str(Params.SNR),'[dB], Feedback SNR = ',num2str(Params.SNR + Params.deltaSNR),'[dB]'))
    powerAnalysis(reshape(Debug(3,i,:,:),Params.T,[]),Params.N_avg,'Kochman',sigma_z)
    powerAnalysis(reshape(Debug(4,i,:,:),Params.T,[]),Params.N_avg,'Tuncel',sigma_z)
    
end

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
