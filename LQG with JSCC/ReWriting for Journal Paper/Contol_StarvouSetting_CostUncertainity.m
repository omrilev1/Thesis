%% Simulation of LQG over power AWGN channel, with receiver SI
% We assume full knowledge of the control signal u_t in the observer, and
% the existence of receiver SI signal, z_t.
% we simulate 3 cases :
% - Observer has full access to the SI
% - Only controller has the control signal
% - noisy feedback between both

% The sensor has access to the plant, and needs to transmit the current
% state to the controller who controls the plant. The sensor has acces to
% x_t , u_t and limited access to z_t , while the controller has access to u_t and z_t, but no
% access to x_t

% The evolution of x_t is given by :
% x_(t+1) = alpha*x_t + w_(t+1) + u_t
% where w_(t+1) is wgn with variance W

% We assume that only the receiver aware of the current cost coefficents
% Q,R , and the transmitter need to design control signals accordingly
% (with cost uncertainity)

close all; clear all; clc;

% Init parameters and arrays structures
Params = initParams_Control();
Params.NumOfSchemes = 3;
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
                
                P_X_linear = Params.P; P_X_Tuncel = Params.P;
                P_X_modulo_feedback = Params.P; P_X_modulo_feedback_linear = Params.P;
            else
                %% Full Access : Subtract the receiver SI
                SI = Params.rho * sqrt(Params.P/PredictErr(1,i,t,n)) * (xVec(1) - xVecPredict(1)) + eta;
                powerNorm = (1 - Params.rho^2) * PredictErr(1,i,t,n) ;
                aVec(1) = sqrt(Params.P/powerNorm)*(xVec(1) - xVecPredict(1) - Params.rho*sqrt(PredictErr(1,i,t,n)/Params.P)*SI);
                
                %% Linear Precoding :
                % Send X_{t} as is, with proper power normalization
                curr_k_for_power = Params.k_for_power(t);
                P_X_linear = Params.W + P_X_linear * (Params.alpha - curr_k_for_power)^2 + ...
                    (curr_k_for_power^2 + 2*curr_k_for_power*(Params.alpha - curr_k_for_power)) * EstimErr(2,i,t-1,n);
                aVec(2) = sqrt(Params.P/P_X_linear) * xVec(2);
                
                %% Tuncel scheme : quantize and send the quantized value + quantization value after proper scaling
                % Tuncel coding
                P_X_Tuncel = Params.W + P_X_Tuncel * (Params.alpha - curr_k_for_power)^2 + ...
                    (curr_k_for_power^2 + 2*curr_k_for_power*(Params.alpha - curr_k_for_power)) * EstimErr(3,i,t-1,n);

                
                normX_Tuncel = sqrt(Params.P/P_X_Tuncel) * xVec(3);
                currDistance = abs(normX_Tuncel - Params.centersTuncel(:));
                [~,minIdx] = min(currDistance,[],1);
                T_Tuncel = Params.codebookTuncel(minIdx);
                S_Tuncel = normX_Tuncel - T_Tuncel;
                
                aVec(3) = Params.alphaTuncel*T_Tuncel + Params.betaTuncel*S_Tuncel;
                
            end
            Debug(:,i,t,n) = aVec(:).^2;
            
            %% channel
            z = sigma_z*randn;
            eta = sqrt(1 - Params.rho^2) * randn;
            bVec = aVec + z*ones(Params.NumOfSchemes,1);
            
            %% controller
            
            % estimation
            if t==1
                xVecEstim = alphaMMSE*sqrt(Params.W/Params.P)*bVec;
                EstimErr(:,i,t,n) = (Params.W/(1+Params.snrLin(i)))*ones(Params.NumOfSchemes,1);
            else
                %% Full Access, Linear Receiver
                
                % calc LMMSE Params
                Cxy = [sqrt(Params.P*(1 - Params.rho^2)*PredictErr(1,i,t,n)) Params.rho * sqrt(Params.P*PredictErr(1,i,t,n))];
                Cyy = [Params.P*(1 +  1/Params.snrLin(i)),0 ; 0,Params.rho^2 * Params.P + (1 - Params.rho^2)];
                
                xVecEstim(1) = xVecPredict(1) + Cxy * (Cyy \ [bVec(1);SI]);
                % EstimErr(1,i,t,n) = 10^(-3.95/20) * PredictErr(1,i,t,n)/(1+Params.snrLin(i));
                EstimErr(1,i,t,n) = (1 - Params.rho^2) * PredictErr(1,i,t,n)/(1+Params.snrLin(i));
                
                %% zero Access, Linear Receiver
                
                % generate SI
                SI = Params.rho * (xVec(2) - xVecPredict(2)) + eta * sqrt(P_X_linear);
                
                % generate the receiver side signals for estimation
                b_tilde_linear = bVec(2) - sqrt(Params.P/P_X_linear) * xVecPredict(2);
                
                % calc LMMSE Params
%                Cxy = [sqrt(Params.P/P_X_linear)*PredictErr(2,i,t,n) Params.rho*PredictErr(2,i,t,n)];
%                Cyy = [Params.P*(1 + 1/Params.snrLin(i)) Params.rho*sqrt(Params.P/P_X_linear)*PredictErr(2,i,t,n) ;...
%                     Params.rho*sqrt(Params.P/P_X_linear)*PredictErr(2,i,t,n) ((Params.rho)^2 * PredictErr(2,i,t,n) + (1 - Params.rho^2)*P_X_linear)];
               Cxy = PredictErr(2,i,t,n) * [sqrt(Params.P/P_X_linear) Params.rho];
               Cyy = [Params.P*(PredictErr(2,i,t,n)/P_X_linear + 1/Params.snrLin(i)) Params.rho*sqrt(Params.P/P_X_linear)*PredictErr(2,i,t,n) ;...
                    Params.rho*sqrt(Params.P/P_X_linear)*PredictErr(2,i,t,n) ((Params.rho)^2 * PredictErr(2,i,t,n) + P_X_linear*(1 - Params.rho^2))]; 

                
                xVecEstim(2) = xVecPredict(2) + Cxy * (Cyy \ [b_tilde_linear;SI]);
                
                EstimErr(2,i,t,n) = PredictErr(2,i,t,n) - Cxy * (Cyy \ Cxy.');
%                 EstimErr(2,i,t,n) = 10^(-2.2/20) * EstimErr(2,i,t,n); % Randomization, 8[dB]
                EstimErr(2,i,t,n) = 10^(-2.85/20) * EstimErr(2,i,t,n); % Worst Case, 8[dB]

                %% Zero Access, Tuncel Receiver

                RxSig_tuncel = bVec(3);
                
                % generate SI
%                 SI = Params.rho * xVec(3) + eta * sqrt(P_X_Tuncel);
%                 receiverStateSig = (Params.alpha - curr_k_for_power) * xVecEstim(3);
%                x_hat_Tuncel = FullMMSE_decoder(RxSig_tuncel,SI,Params.alphaTuncel,Params.betaTuncel,Params.rho,sigma_z,...
%                     Params.codebookTuncel,Params.centersTuncel(:),ones(size(Params.codebookTuncel)),P_X_Tuncel,Params.P);
%                 
%                 xVecEstim(3) = xVecPredict(3) + sqrt(P_X_Tuncel/Params.P) * (x_hat_Tuncel) - receiverStateSig;
%                EstimErr(3,i,t,n) = 10^(-1/20) * (1 - Params.rho^2) * PredictErr(3,i,t,n)/(1 + (1 - Params.rho^2) * Params.snrLin(i));

                SI = Params.rho * normX_Tuncel + eta;
               x_hat_Tuncel = FullMMSE_decoder(RxSig_tuncel,SI,Params.alphaTuncel,Params.betaTuncel,Params.rho,sigma_z,...
                    Params.codebookTuncel,Params.centersTuncel(:),ones(size(Params.codebookTuncel)));
                
                xVecEstim(3) = sqrt(P_X_Tuncel/Params.P)* (x_hat_Tuncel);   
                
               EstimErr(3,i,t,n) = 10^(-1/20) * (1 - Params.rho^2)/(1 + (1 - Params.rho^2) * Params.snrLin(i));
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
    lineStyle = ['-g','--k',':c',':g','b',':m'];
    for k = [2 3 1]
        currMean = mean(reshape(Results(k,i,:,:),Params.T,[]),2);
        plot(1:(Params.T-1),10*log10(currMean(1:end-1)),lineStyle(k),'LineWidth',2)
    end
    % Bound From Toli Paper
    plot(1:Params.T,10*log10(Params.s_vec(2)*Params.W + ((Params.Q + (Params.alpha^2-1)*Params.s_vec(2))/((1+Params.snrLin(i))/(1 - Params.rho^2)-Params.alpha^2))*Params.W)*ones(1,Params.T),'-m','LineWidth',2)
    plot(1:Params.T,10*log10(Params.s_vec(2)*Params.W + ((Params.Q + (Params.alpha^2-1)*Params.s_vec(2))/((1+10^(0.15/20)*(1 - Params.rho^2)*Params.snrLin(i))/(1 - Params.rho^2)-Params.alpha^2))*Params.W)*ones(1,Params.T),'-','LineWidth',2)
    plot(1:Params.T,10*log10(Params.s_vec(2)*Params.W + ((Params.Q + (Params.alpha^2-1)*Params.s_vec(2))/((1+10^(3.05/20)*((1 - Params.rho^2))*Params.snrLin(i))/(1 - Params.rho^2)-Params.alpha^2))*Params.W)*ones(1,Params.T),'-','LineWidth',2)
    
    grid on; grid minor;
    xlabel('t'); ylabel('cost [dB]');
    legend('Observer is Linear',...
        strcat('Observer With Tuncel : \alpha = ',num2str(Params.alphaTuncel),' \beta = ',num2str(Params.betaTuncel)),...
        'SI in Tx and Rx','LQG_{\infty} With Tx and Rx SI','LQG_{\infty} With Rx SI and Linear Encoder','LQG_{\infty} With Rx SI and Modulo Encoder');
    %         strcat('Linear reversed feedback : ',num2str(Params.N_feedback),' iterations'),'SDR OPTA');
    
    title(strcat('SNR = ',num2str(Params.SNR),'[dB], SI noise is \sigma^{2}_z = ',num2str(1)))
    powerAnalysis(reshape(Debug(1,i,:,:),Params.T,[]),Params.N_avg,'SI Known to Tx and Rx',sigma_z)
    powerAnalysis(reshape(Debug(2,i,:,:),Params.T,[]),Params.N_avg,'Linear Precoding',sigma_z)
    powerAnalysis(reshape(Debug(3,i,:,:),Params.T,[]),Params.N_avg,'Tuncel',sigma_z)
    
    
end

function [] = powerAnalysis(P_x,N_avg,run,sigma)

inputPower = sum(P_x,2)/N_avg;

figure;hold all
plot(inputPower)
title([run strcat('measured SNR = ',num2str(10*log10(mean(inputPower(2:end-1))/sigma^2)))])

end

function [T_hat] = MMSE_decoder_Int(RxSig_tuncel,y,delta,codebook,alpha,beta,rho,sigmaW,P_x,P)

ds = 0.001;
s = -delta/2:ds:delta/2;

% start with 2D calculations
calc1 = codebook(:) + s;
calc2 = alpha*codebook(:) + beta*s;

expTerm1 = -1*(0.5*(calc1/sqrt(P)).^2);
expTerm2 = - (0.5/(P_x * (1-rho^2)))*(y - rho*sqrt(P/P_x)*calc1).^2;
expTerm3 = 0.5*(1/sigmaW^2)*(RxSig_tuncel - calc2).^2;

Integrand = exp(expTerm1 + expTerm2 - expTerm3);

num = sum(Integrand,2);
T_hat = sum(codebook(:).*num(:),1)./sum(num(:));

% quantize to the nearset r_k
% [~,minIdx] = min(abs(T_hat - codebook));
% T_hat = codebook(minIdx);
end

function [S_hat] = MMSE_decoder_Frac(RxSig_tuncel,y,delta,alpha,beta,rho,sigmaW,codebook,P_x,P)


ds = 0.001;
s = ones(length(codebook),1) * (-delta/2:ds:delta/2);

% start with 2D calculations
calc1 = codebook(:) + s;
calc2 = alpha*codebook(:) + beta*s;

expTerm1 = -1*(0.5*(calc1/sqrt(P)).^2);
expTerm2 = - (0.5/(P_x * (1-rho^2)))*(y - rho*sqrt(P/P_x)*calc1).^2;
expTerm3 = 0.5*(1/sigmaW^2)*(RxSig_tuncel - calc2).^2;
Integrand = exp(expTerm1 + expTerm2 - expTerm3);


num = sum(sum(Integrand.*s,2),1);
den = sum(sum(Integrand));
S_hat = num/den;

end

function [X_hat] = FullMMSE_decoder(RxSig,y,alpha,beta,rho,sigmaW,codebook,centers,probs)

dx = 0.001;
x = (-8:dx:8);

% start with 2D calculations
currDistance = abs(x - centers(:));
[~,minIdx] = min(currDistance,[],1);
T = codebook(minIdx);
S = x - T;

% Tuncel 
f_x = alpha*T + beta*probs(minIdx).*S;

expTerm1 = -1*(0.5*(x).^2);
expTerm2 = - (0.5/(1-rho^2))*(y - rho*x).^2;
expTerm3 = - 0.5*(1/sigmaW^2)*(RxSig - f_x).^2;
Integrand = exp(expTerm1 + expTerm2 + expTerm3);

num = sum(Integrand.*x,2);
den = sum(Integrand);
X_hat = num/den;

end



