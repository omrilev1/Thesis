%% Simulation of Analog PPM + MLM JSCC in the infinite bandwidth regime
% We simulate the PPM in an artificial manner, by generating the
% equivalence channel noise

close all; clear all; clc;

% Init parameters and arrays structures
Params = initParams();

E_0 = 1;
optBeta = 0;% -1.63;
gamma_PPM = 10^(optBeta/10)* Params.snrLin(1);

D_Linear = zeros(length(Params.SNR),1);
D_PPM = zeros(length(Params.SNR),1);
D_LinearPrev = zeros(length(Params.SNR),1);
D_PPMPrev = zeros(length(Params.SNR),1);

currNumOfLevels = 1;finished = 0;
currE = zeros(1,currNumOfLevels);currE(1) = 1;

etaLinear = 1; etaPPM = 1;

for i=1:length(Params.SNR)
    
    prevNumOfLevels = currNumOfLevels;
    currNumOfLevels = calcNumOfLevels(Params.SNR,i);
    
    if currNumOfLevels==1
        N_avg = 2^16;
    else
        N_avg = Params.N_avg;
    end
    
    if prevNumOfLevels ~= currNumOfLevels
        gamma_PPM = 10^(optBeta/10)*Params.snrLin(i-1);
    end
    
    % calculate semi-analytical (Wozencraft+Jacobs) PPM MSE, to derive mmse
    % estimation coefficient
%     c_opt = ((144/pi^2) * sqrt(2*pi*gamma_PPM)/(2*gamma_PPM))^(1/3);
    c_opt = ((52*6*sqrt(pi))^(1/3))/(gamma_PPM^(5/6));
    beta = c_opt * exp(gamma_PPM/6);
    W = beta/2;
    
    dt = 1/(100*beta); Fs = 1/dt;
    
    if gamma_PPM < 7.5
        t = -2.5:dt:2.5;
    else
        t = -1.125:dt:1.125;
    end
    tIdx = find(abs(t) <= 0.5);
    
%     ppmPulse = sqrt(2*W)*sin(2*pi*W*t)./(2*pi*W*t);
%     ppmPulse = ppmPulse / sqrt(sum(abs(ppmPulse.^2)*dt));
    ppmPulse = zeros(size(t));
    ppmPulse(abs(t) < 1/(2*beta)) = sqrt(beta);
    ppmPulse = ppmPulse / sqrt(sum(abs(ppmPulse.^2)*dt));

    Lx = length(ppmPulse);
    Ly = length(ppmPulse)+length(ppmPulse)-1;
    Ly2 = pow2(nextpow2(Ly));
    PPMfreq = fft(ppmPulse, Ly2);
    
    currDist_Linear = 0;
    currDist_PPM = 0;
    
    DistVec_Linear = zeros(1,currNumOfLevels);
    DistVec_PPM = zeros(1,currNumOfLevels);
    
    energy = 0;energy_PPM = 0;
    
    dx = 1e-4;
    x = -12.5:dx:12.5;
    
    Pe = 1 - (1/sqrt(2*pi))*dx*sum(((1 - qfunc(x + sqrt(2*Params.snrLin(i)))).^((W-1))).*exp(-x.^2/2));
    msePPM = (1 - Pe).*(1./(2*Params.snrLin(i)))*(3/(pi^2))*(1/beta^2) + Pe.* (1/6);
    
    if prevNumOfLevels ~= currNumOfLevels
        
        currEtaLinear = optimizeEta(D_Linear(i-1),Params.IntLin(currNumOfLevels-1),1/Params.snrLin(i));
        currEtaPPM =  optimizeEta(D_PPM(i-1),Params.IntLin(currNumOfLevels-1),msePPM);
        currEnergy = calcModuloEnergy(Params.IntLin(currNumOfLevels-1),currEtaLinear);
        
        etaLinear = [etaLinear currEtaLinear];
        etaPPM = [etaPPM currEtaPPM];
        
        currE = [currE currEnergy];
        %         for k=1:currNumOfLevels
        %             if k>1
        %                 Params.etaLinear(k-1) = optimizeEta(D_Linear(i-1),Params.IntLin(k-1),1/Params.snrLin(i));
        %                 Params.etaPPM(k-1) = optimizeEta(D_PPM(i-1),Params.IntLin(k-1),msePPM);
        %                 currE(k) = calcModuloEnergy(Params.IntLin(k-1),Params.etaLinear(k-1));
        %             else
        %                 currE(k) = 1;
        %             end
    end
    
    for n = 1 : N_avg
        
        % calculate the sources of the current transmission levels
        S = rand - 0.5;
        TxSourceLinear = zeros(1,Params.maxStages);
        TxSourcePPM = zeros(1,Params.maxStages);
        dither = zeros(1,Params.maxStages);
        
        TxSourceLinear(1) = S;
        TxSourcePPM(1) = S;
        
        if currNumOfLevels > 1
            for idx=2:currNumOfLevels
                delta = Params.IntLin(idx-1);
                TxSourceLinear(idx) = mod(etaLinear(idx)*S + delta,2*delta) - delta;
                TxSourcePPM(idx) = mod(etaPPM(idx)*S + delta,2*delta) - delta;
            end
            
        end
        %         energy = energy + TxSourceLinear.^2;
        %         energy_PPM = energy_PPM + TxSourcePPM.^2;
        % simulate the successive transmission and reception
        
        y_linear = zeros(1,currNumOfLevels);
        y_ppm = zeros(1,currNumOfLevels);
        
        sHat_Linear = 0;
        sHat_PPM = 0;
        for k=1:currNumOfLevels
            %
            %             if k==1
            %                 currE = 1;
            %             else
            %                 Params.etaLinear(k-1) = optimizeEta(D_Linear(i-1),Params.IntLin(k-1),1/Params.snrLin(i));
            %                 etaPPM(k-1) = optimizeEta(D_PPM(i-1),Params.IntLin(k-1),msePPM(k));
            %                 currE = calcModuloEnergy(Params.IntLin(k-1),Params.etaLinear(k-1));
            %
            %             end
            
            %% channel part
            if k==1
                % Linear
                y_linear(1) = sqrt(12) * sqrt(currE(k)) * S + sqrt(1/Params.snrLin(i))*randn;
                
                % PPM
%                 TxPulse = sqrt(2*W)*sin(2*pi*W*(t - S))./(2*pi*W*(t - S));
%                 TxPulse(t == S) = max(ppmPulse);
%                 TxPulse = TxPulse/sum(abs(TxPulse.^2)*dt);
                TxPulse = zeros(size(t));
                TxPulse(abs(t - S) < 1/(2*beta)) = sqrt(beta);
                TxPulse = TxPulse/sum(abs(TxPulse.^2)*dt);
        
                r_ppm = sqrt(currE(k)) * TxPulse + sqrt((1/(2*dt))/Params.snrLin(i))*randn(size(t));
                
%                 % ppm correlator receiver (ML)
%                 y_ppm(1) = ppmCorrelator(r_ppm,PPMfreq,tIdx,t,Lx,Ly,Ly2);
%                 
                
                % PPM Correlator receiver
                PPMcorr = fconv(r_ppm,Lx,Ly,Ly2,PPMfreq);
                
                %         PPMcorr2 = conv(r,ppmPulse,'same');
                PPMcorr = PPMcorr(tIdx);
                
                [~,maxIdx] = max(PPMcorr);
                y_ppm(1) = t(min(tIdx) + maxIdx - 1);
                if y_ppm(1) > 0.5
                    y_ppm(1) = 0.5;
                elseif y_ppm(1) <= -0.5
                    y_ppm(1) = -0.5;
                end
                
            else
                % Linear
                y_linear(k) = TxSourceLinear(k) + sqrt(1/Params.snrLin(i))*randn;
                
                % PPM
                TxPulse = sqrt(currE(k))*sqrt(2*W)*sin(2*pi*W*(t - TxSourcePPM(k)/(2*Params.IntLin(k-1))))./(2*pi*W*(t - TxSourcePPM(k)/(2*Params.IntLin(k-1))));
                r_ppm = TxPulse + sqrt((1/(2*dt))/Params.snrLin(i))*randn(size(t));
                
                % ppm correlator receiver (ML)
                % y_ppm(k) = ppmCorrelator(r_ppm,PPMfreq,tIdx,t,Lx,Ly,Ly2);
                                % PPM Correlator receiver
                PPMcorr = fconv(r_ppm,Lx,Ly,Ly2,PPMfreq);
                
                %         PPMcorr2 = conv(r,ppmPulse,'same');
                PPMcorr = PPMcorr(tIdx);
                
                [~,maxIdx] = max(PPMcorr);
                y_ppm(k) = t(min(tIdx) + maxIdx - 1);
                if y_ppm(k) > 0.5
                    y_ppm(k) = 0.5;
                elseif y_ppm(k) <= -0.5
                    y_ppm(k) = -0.5;
                end
                y_ppm(k) = y_ppm(k) * (2*Params.IntLin(k-1));
            end
            
            %% MLM Part
            if k==1
                sHat_Linear = y_linear(1) / (sqrt(currE(k)) * sqrt(12)) ;
                sHat_PPM = y_ppm(1);
            else
                
                alphaMMSE_Linear = currE(k)/(currE(k) + 1/Params.snrLin(i));
                sHat_Linear = sHat_Linear + (1/etaLinear(k)) * (mod(alphaMMSE_Linear * y_linear(k) - etaLinear(k)*sHat_Linear + Params.IntLin(k-1),2*Params.IntLin(k-1)) - Params.IntLin(k-1));
                
                alphaMMSE_PPM = currE(k)/(currE(k) + msePPM);
                sHat_PPM = sHat_PPM + (1/etaPPM(k)) * (mod(alphaMMSE_PPM * y_ppm(k) - etaPPM(k)*sHat_PPM + Params.IntLin(k-1),2*Params.IntLin(k-1)) - Params.IntLin(k-1));
                
            end
        
        DistVec_Linear(k) = (sHat_Linear - S)^2 + DistVec_Linear(k);
        DistVec_PPM(k) = (sHat_PPM - S)^2 + DistVec_PPM(k);   
            
        end
        
        currDist_Linear = (sHat_Linear - S)^2 + currDist_Linear;
        currDist_PPM = (sHat_PPM - S)^2 + currDist_PPM;
    end
    %     energy = energy/N_avg;
    %     energy_PPM = energy_PPM/N_avg;
    D_Linear(i) = currDist_Linear/N_avg;
    D_PPM(i) = currDist_PPM/N_avg;
    
    if k>1
        D_LinearPrev(i) = DistVec_Linear(end-1)/N_avg;
        D_PPMPrev(i) = DistVec_PPM(end-1)/N_avg;
    else
        D_LinearPrev(i) = D_Linear(i);
        D_PPMPrev(i) = D_PPM(i);
    end
    
    if mod(i,5) == 0
        disp(strcat('Finished SNR = ',num2str(Params.SNR(i))));
    end
    
end


figure;hold all
plot(Params.SNR(1:i),10*log10(1/12./((Params.snrLin(1:i)/Params.snrLin(1)).^2)),'LineWidth',2)
plot(Params.SNR(1:i),10*log10(1/12./((Params.snrLin(1:i)/Params.snrLin(1)).^3)),'LineWidth',2)
plot(Params.SNR(1:i),10*log10(D_Linear(1:i)),'-o','LineWidth',2)
plot(Params.SNR(1:i),10*log10(min(D_PPM(1:i),D_Linear(1:i))),'-*','LineWidth',2)
xlabel('ENR [dB]'); ylabel('Distortion [dB]');
legend('Profile (L = 2)','Profile (L = 3)','Linear + MLM','PPM + MLM');


function [sHat] = ppmCorrelator(r,PPMfreq,tIdx,t,Lx,Ly,Ly2)
% PPM Correlator receiver
% PPMcorr = conv(r,ppmPulse,'same');
% PPMcorr = PPMcorr(tIdx);
%
% [~,maxIdx] = max(PPMcorr);
% sHat = t(min(tIdx) + maxIdx - 1);
% if sHat > 0.5
%     sHat = 0.5;
% elseif sHat <= -0.5
%     sHat = -0.5;
% end

% PPM Correlator receiver
PPMcorr = fconv(r,Lx,Ly,Ly2,PPMfreq);

%         PPMcorr2 = conv(r,ppmPulse,'same');
PPMcorr = PPMcorr(tIdx);

[~,maxIdx] = max(PPMcorr);
sHat = t(min(tIdx) + maxIdx - 1);
if sHat > 0.5
    sHat = 0.5;
elseif sHat <= -0.5
    sHat = -0.5;
end
end


function currNumOfLevels = calcNumOfLevels(SNR,i)

snrLin = 10.^(SNR/10);
profile = 1./((snrLin/snrLin(1)).^2);

numOfLevels = 1;
prevSNR = 1;
for l=1:i
    if prevSNR/snrLin(l) >= profile(l)
        numOfLevels = numOfLevels + 1;
        prevSNR = prevSNR * 10^(-1.7/10);
    end
end
currNumOfLevels = numOfLevels;
end

function Energy = calcModuloEnergy(Int,eta)

% calculate modulo indices
dx = Int/1e4;
x = -15*Int:dx:(15*Int - dx);
Int_dx = Int/dx;
mod_Idx = reshape(1:length(x),2*Int_dx,[]);

% calculate equivalent input distribution
f_s = zeros(size(x));
f_s(abs(x) < 0.5*eta) = 1;
f_s = f_s/sum(f_s);

f_mod_s = sum(f_s(mod_Idx),2);
f_mod_s = f_mod_s / sum(f_mod_s);
s_mod = -Int:dx:(Int - dx);

Energy = sum((s_mod(:)).^2.*f_mod_s);

end

function etaOpt = optimizeEta(D,Int,sigmaN)
eta = 1:0.0125:12;

% calculate probability of error
P_e = qfunc(Int./sqrt(eta.^2*D + sigmaN));

% calculate distortion
dist = (1 - P_e).*(sigmaN./(eta.^2)) + P_e.*(2*Int^2);

[~,etaOpt_idx] = min(dist);
etaOpt = eta(etaOpt_idx);
end

function y = fconv(x,Lx,Ly,Ly2,H)
% Convolution in frequency domain using power of 2 fft
% Since the input signal is real we use known fft identities to accelerate
% the fft

% input fft
X = fft(x, Ly2);

% multiply with precalculated freq domain signal
Y = X.*H;

% inverse fft and truncation
y = real(ifft(Y, Ly2));
y = y(1:1:Ly);

if mod(Lx,2) == 0
    y = y(Lx/2 + 1 : end - (Lx/2) + 1);
else
    y = y((Lx+1)/2 : end - ((Lx+1)/2) + 1);
end

end
