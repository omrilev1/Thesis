%% Simulation of LQG, with tracking over reference input, over power AWGN channel, with receiver SI
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
% we extend the model to account also for the tracking error, thus the
% second state equation is given by
% e_(t+1) = e_(t) + (x_(t) - r)
% where w_(t+1) is wgn with variance W

% We assume that only the receiver aware of the current cost coefficents
% Q,R , and the transmitter need to design control signals accordingly
% (with cost uncertainity)

close all; clear all; clc;

% Init parameters and arrays structures
Params = initParams_LQGtracking();
Params.NumOfSchemes = 3;
[Results,Debug,PredictErr,EstimErr] = initStructs_Control(Params);

for i=1:length(Params.SNR)
    sigma_z = sqrt(Params.P/Params.snrLin(i));
    sigma_v = sqrt(Params.P/Params.snrFeedback_Lin(i));
    
    xVecEstim = zeros(2,Params.NumOfSchemes); xVecPredict = zeros(2,Params.NumOfSchemes);
    xVec = zeros(2,Params.NumOfSchemes);
    aVec = zeros(2,Params.NumOfSchemes);bVec = zeros(2,Params.NumOfSchemes);
    uVec = zeros(1,Params.NumOfSchemes);
    
    alphaMMSE = Params.snrLin(i)/(1+Params.snrLin(i));
    
    x_path = zeros(Params.N_avg,Params.T,3); err_path = zeros(Params.N_avg,Params.T,3);
    for n = 1 : Params.N_avg
        ref = 10; % Params.maxRefInput*(rand - 0.5);
        
        for t=1:Params.T
            
            % The input driving the schemes is the same input
            curr_w = sqrt(Params.W)*randn;
            %% plant
            if t==1
                % initialization - x_{0} = w_{0}
                xVec = [curr_w*ones(1,Params.NumOfSchemes);ref*ones(1,Params.NumOfSchemes)];
                
            elseif t==Params.T
                % no control in the last stage - evolve state and break
                xVec = [Params.alpha 0;-1 1]*xVec + [curr_w;0]*ones(1,Params.NumOfSchemes) + [1;0]*uVec + [0;1]*ref;
                Results(:,i,t,n) = (Results(:,i,t-1,n)*(Params.T-1) + Params.Q*(xVec(1,:).^2).') / Params.T; % (Results(:,i,t-1,n)*(Params.T-1) + Params.Q*(xVec(1,:).^2).') / Params.T;
                break
            else
                % evolution
                xVec = [Params.alpha 0;-1 1]*xVec + [curr_w;0]*ones(1,Params.NumOfSchemes) + [1;0]*uVec + [0;1]*ref;
                
            end
            x_path(n,t,:) = xVec(1,:);
            err_path(n,t,:) = xVec(2,:);
            %% observer
            if t==1
                % in the first stage, we just normalize the current instant and
                % send to the controller
                aVec = sqrt(Params.P/Params.W) * xVec(1,:);
                Debug(:,i,t,n) = Params.W;
                
                Corr_Mat_Linear = [Params.W 0; 0 ref^2]; Corr_Mat_Tuncel = [Params.W 0; 0 ref^2];
                Expec_Vec_Linear = [0;ref]; Expec_Vec_Tuncel = [0;ref];
            else
                %% Full Access : Subtract the receiver SI
                aVec(1) = sqrt(Params.P/PredictErr(1,i,t,n))*(xVec(1,1) - xVecPredict(1,1));
                
                %% Linear Precoding :
                % Send X_{t} as is, with proper power normalization
                curr_k_for_power = reshape(Params.k_for_power(1,:,t),1,2);
                
                Expec_Vec_Linear = (Params.Amat - Params.Bmat * curr_k_for_power) * Expec_Vec_Linear + [0;1]*ref;
%                 if t > 2
%                     Corr_Mat_Linear = [Params.W 0; 0 0] + ...
%                         (Params.Amat - Params.Bmat * curr_k_for_power) * Corr_Mat_Linear * (Params.Amat - Params.Bmat * curr_k_for_power).' + ...
%                         Params.Bmat * curr_k_for_power * [EstimErr(2,i,t-1,n) 0;0 EstimErr(2,i,t-2,n)] * (Params.Bmat * curr_k_for_power).'; %  + ref^2 * [0 0;0 1];
%                 else
%                     Corr_Mat_Linear = [Params.W 0; 0 0] + ...
%                         (Params.Amat - Params.Bmat * curr_k_for_power) * Corr_Mat_Linear * (Params.Amat - Params.Bmat * curr_k_for_power).' + ...
%                         Params.Bmat * curr_k_for_power * [EstimErr(2,i,t-1,n) 0;0 0] * (Params.Bmat * curr_k_for_power).'; %  + ref^2 * [0 0;0 1];
%                 end
%                 Corr_Mat_Linear = Corr_Mat_Linear - [0 0;0 1]*ref^2 + ref*[0 Expec_Vec_Linear(1);Expec_Vec_Linear(1) 2*Expec_Vec_Linear(2)];
                
                if t > 2
                    Rn = [EstimErr(2,i,t-1,n) 0;0 sum(EstimErr(2,i,1:t-2,n))];
                else
                    Rn = [EstimErr(2,i,t-1,n) 0;0 0];
                end
                Corr_Mat_Linear = [Params.W 0; 0 0] + ...
                        (Params.Amat - Params.Bmat * curr_k_for_power) * Corr_Mat_Linear * (Params.Amat - Params.Bmat * curr_k_for_power).' + ...
                        Params.Bmat * curr_k_for_power * Rn * (Params.Bmat * curr_k_for_power).' + (Params.Amat - Params.Bmat * curr_k_for_power) * Rn * (Params.Bmat * curr_k_for_power).' + ...
                        Params.Bmat * curr_k_for_power * Rn * (Params.Amat - Params.Bmat * curr_k_for_power).';
                    
                Corr_Mat_Linear = Corr_Mat_Linear - [0 0;0 1]*ref^2 + ref*[0 Expec_Vec_Linear(1);Expec_Vec_Linear(1) 2*Expec_Vec_Linear(2)];
                
%                 P_X_linear = 0.925 * Corr_Mat_Linear(1,1); % (Q,R)
                P_X_linear = 0.925 * Corr_Mat_Linear(1,1); % (Q,R)
                
                aVec(2) = sqrt(Params.P/P_X_linear) * xVec(1,2);
                
                r_Linear = sqrt(1 - PredictErr(2,i,t,n) / P_X_linear);
                
                %% Tuncel scheme : quantize and send the quantized value + quantization value after proper scaling
                % Tuncel coding
                Expec_Vec_Tuncel = (Params.Amat - Params.Bmat * curr_k_for_power) * Expec_Vec_Tuncel + [0;1]*ref;
                
%                 if t > 2
%                     Corr_Mat_Tuncel = [Params.W 0; 0 0] + ...
%                         (Params.Amat - Params.Bmat * curr_k_for_power) * Corr_Mat_Tuncel * (Params.Amat - Params.Bmat * curr_k_for_power).' + ...
%                         Params.Bmat * curr_k_for_power * [EstimErr(3,i,t-1,n) 0;0 EstimErr(3,i,t-2,n)] * (Params.Bmat * curr_k_for_power).'; %  + ref^2 * [0 0;0 1];
%                 else
%                     Corr_Mat_Tuncel = [Params.W 0; 0 0] + ...
%                         (Params.Amat - Params.Bmat * curr_k_for_power) * Corr_Mat_Tuncel * (Params.Amat - Params.Bmat * curr_k_for_power).' + ...
%                         Params.Bmat * curr_k_for_power * [EstimErr(3,i,t-1,n) 0;0 0] * (Params.Bmat * curr_k_for_power).'; %  + ref^2 * [0 0;0 1];
%                 end
                if t > 2
                    Rn = [EstimErr(3,i,t-1,n) 0;0 sum(EstimErr(3,i,1:t-2,n))];
                else
                    Rn = [EstimErr(3,i,t-1,n) 0;0 0];
                end
                Corr_Mat_Tuncel = [Params.W 0; 0 0] + ...
                        (Params.Amat - Params.Bmat * curr_k_for_power) * Corr_Mat_Tuncel * (Params.Amat - Params.Bmat * curr_k_for_power).' + ...
                        Params.Bmat * curr_k_for_power * Rn * (Params.Bmat * curr_k_for_power).' + (Params.Amat - Params.Bmat * curr_k_for_power) * Rn * (Params.Bmat * curr_k_for_power).' + ...
                        Params.Bmat * curr_k_for_power * Rn * (Params.Amat - Params.Bmat * curr_k_for_power).';
                    
                Corr_Mat_Tuncel = Corr_Mat_Tuncel - [0 0;0 1]*ref^2 + ref*[0 Expec_Vec_Tuncel(1);Expec_Vec_Tuncel(1) 2*Expec_Vec_Tuncel(2)];
                
                P_X_Tuncel = 1.65 * Corr_Mat_Tuncel(1,1);
                
                normX_Tuncel = sqrt(Params.P/P_X_Tuncel) * xVec(1,3);
                currDistance = abs(normX_Tuncel - Params.centersTuncel(:));
                [~,minIdx] = min(currDistance,[],1);
                T_Tuncel = Params.codebookTuncel(minIdx);
                S_Tuncel = normX_Tuncel - T_Tuncel;
                
                aVec(3) = Params.alphaTuncel*T_Tuncel + Params.betaTuncel*S_Tuncel;
                
                r_Tuncel = sqrt(1 - PredictErr(3,i,t,n) / P_X_Tuncel);
            end
            Debug(:,i,t,n) = aVec(:).^2;
            
            %% channel
            z = sigma_z*randn;
            eta = sqrt(1 - Params.rho^2) * randn;
            bVec = aVec + z*ones(1,Params.NumOfSchemes);
            
            %% controller
            
            % estimation
            if t==1
                xVecEstim(1,:) = alphaMMSE*sqrt(Params.W/Params.P)*bVec;
                xVecEstim(2,:) = ref;
                EstimErr(:,i,t,n) = (Params.W/(1+Params.snrLin(i)))*ones(Params.NumOfSchemes,1);
            else
                %% Full Access, Linear Receiver
                
                xVecEstim(2,1) = xVecEstim(2,1) + (ref - xVecEstim(1,1));
                xVecEstim(1,1) = xVecPredict(1,1) + sqrt(PredictErr(1,i,t,n)/Params.P) * alphaMMSE * bVec(1);
                EstimErr(1,i,t,n) = PredictErr(1,i,t,n)/(1+Params.snrLin(i));
                
                %% zero Access, Linear Receiver
                
                % generate the receiver side signals for estimation
                b_tilde_linear = bVec(2) - sqrt(Params.P/P_X_linear) * xVecPredict(1,2);
                
                % calc LMMSE Params
                Cxy = sqrt(Params.P/P_X_linear) * PredictErr(2,i,t,n);
                Cyy = Params.P*(PredictErr(2,i,t,n)/P_X_linear + 1/Params.snrLin(i));
                
                xVecEstim(2,2) = xVecEstim(2,2) + (ref - xVecEstim(1,2));
                xVecEstim(1,2) = xVecPredict(1,2) + (Cxy/Cyy) * b_tilde_linear;
                
                EstimErr(2,i,t,n) = PredictErr(2,i,t,n)* (Params.P / Params.snrLin(i)) / Cyy;
                
                
                %% Zero Access, Tuncel Receiver
                
                SI = r_Tuncel * normX_Tuncel + sqrt(1 - r_Tuncel^2)*randn;
                x_hat_Tuncel = FullMMSE_decoder(bVec(3),SI,r_Tuncel,Params.alphaTuncel,Params.betaTuncel,sigma_z,...
                    Params.codebookTuncel,Params.centersTuncel(:),ones(size(Params.codebookTuncel)));
                estimErr = xVec(1,3) - sqrt(P_X_Tuncel/Params.P) * x_hat_Tuncel;
                
                xVecEstim(2,3) = xVecEstim(2,3) + (ref - xVecEstim(1,3));
                xVecEstim(1,3) = sqrt(P_X_Tuncel/Params.P)* x_hat_Tuncel;
                
                EstimErr(3,i,t,n) = 10^(-1/20) * PredictErr(3,i,t,n)* (Params.P / Params.snrLin(i)) / Cyy;
            end
            
            %% control generation
            uVec = -Params.k_mat(:,:,t) * xVecEstim;
            
            %% update estimates
            xVecPredict = [Params.alpha 0;-1 1]*xVecEstim + [1;0]*uVec(1,:) + [0;1]*ref;
            PredictErr(:,i,t+1,n) = Params.alpha^2 * EstimErr(:,i,t,n) + Params.W;
            
            %% Calculate current cost
            if t == 1
                Results(:,i,t,n) =  Params.Q*(xVec(1,:).').^2 + Params.R*(uVec(1,:).').^2; % Params.Q*(xVec(1,:).^2) + Params.R*(uVec(1,:).^2)
            else
                Results(:,i,t,n) =  (Results(:,i,t-1,n)*(t-1) + (Params.Q*(xVec(1,:).').^2 + Params.R*(uVec(1,:).').^2))/t; % (Results(:,i,t-1,n)*(t-1) + (Params.Q*(xVec(1,:).^2) + Params.R*(uVec(1,:).^2)).')/t
            end
            
        end
        
        if mod(n,10) == 0
            display(strcat('n = ',num2str(n)));
        end
        
    end
    
    figure;hold all
    lineStyle = ['-g','--k',':c',':g','b',':m'];
    for k = [2 3 1]
        currMean = mean(reshape(Results(k,i,:,:),Params.T,[]),2);
        plot(1:(Params.T-1),10*log10(currMean(1:end-1)),lineStyle(k),'LineWidth',2)
    end
    
    % Bound From Rate-Cost tradeoffs in control

    grid on; grid minor;
    xlabel('t'); ylabel('cost [dB]');
    legend('Observer is Linear',...
        strcat('Non Linear Mapping [17] : \alpha = ',num2str(Params.alphaTuncel),' \beta = ',num2str(Params.betaTuncel)),...
        'SI in Tx and Rx');
    %         strcat('Linear reversed feedback : ',num2str(Params.N_feedback),' iterations'),'SDR OPTA');
    
    title(strcat('SNR = ',num2str(Params.SNR),'[dB], ref = ',num2str(ref)))
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

function [X_hat] = FullMMSE_decoder(RxSig,y,rho,alpha,beta,sigmaW,codebook,centers,probs)

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



