%% Analog PPM Performance simulation
clear ; clc; close all

ENR = (-1:1:22);
ENRlin = 10.^(ENR/10);

ENR_opt = 12;
ENR_opt_lin = 10^(ENR_opt/10);

c_opt = ((13*sqrt(pi)/(4*sqrt(15)))^(1/3))/(ENR_opt_lin^(5/6));
beta = 10; 2*c_opt * exp(ENR_opt_lin/6);

W = beta/2;

dt = 1/(200*beta); Fs = 1/dt;

overload = 5.125;
t = -overload:dt:overload;

MSE_MAP = zeros(size(ENR));% Mean = zeros(size(ENR));
LargeErrMSE = zeros(size(ENR));
ProbOfLargeErr = zeros(size(ENR));
Nrun =  1e4;

ppmPulse = zeros(size(t));
ppmPulse(abs(t) < 1/(2*beta)) = sqrt(beta);
ppmPulse = ppmPulse / sqrt(sum(abs(ppmPulse.^2)*dt));
% ppmPulse(abs(t) > 0.5) = 0;

Lx = length(ppmPulse);
Ly = length(ppmPulse)+length(ppmPulse)-1;
Ly2 = pow2(nextpow2(Ly));
PPMfreq = fft(ppmPulse, Ly2);

Hist = zeros(length(ENR),Nrun);
TxEnergy = zeros(1,Nrun);
for i=1:length(ENR)
    
    currMSE_MAP = zeros(1,Nrun); currLargeErrMSE = zeros(1,Nrun);
    largeErrorInd = zeros(1,Nrun);sVec = zeros(1,Nrun);sHatVec = zeros(1,Nrun);
    currENR = ENRlin(i);
    largeError = zeros(1,Nrun);
    parfor n=1:Nrun
        
        % generate source - Gaussian Source with edge truncation
        S = randn;
        if abs(S) > max(t)
            S = max(t)*sign(S);
        end
        
        TxPulse = zeros(size(t));
        TxPulse(abs(t - S) < 1/(2*beta)) = sqrt(beta);
        TxPulse = sqrt(currENR)*TxPulse/sum(abs(TxPulse.^2)*dt);
        
        noise = randn(size(t));
        noise = sqrt(1/(2*dt))*noise;
        r = TxPulse + noise;
        
        % PPM Correlator receiver
        PPMcorr = fconv(r,Lx,Ly,Ly2,PPMfreq);

        [~,maxIdx] = max(sqrt(currENR)*PPMcorr*dt - 0.5*(t.^2));
        sHat_MAP = t(maxIdx);
        
        currMSE_MAP(n) = (S - sHat_MAP)^2;
        
        if abs(S - sHat_MAP) >= 1/beta
            currLargeErrMSE(n) = (S - sHat_MAP)^2;
            largeErrorInd(n) = 1;
            largeError(n) = sHat_MAP;
        else
            currLargeErrMSE(n) = 0;
        end
        sVec(n) = S;
        sHatVec(n) = sHat_MAP;
    end
    
    MSE_MAP(i) = sum(currMSE_MAP)/Nrun;
    
    if sum(largeErrorInd > 0) == 0
        LargeErrMSE(i) = 0;
    else
        LargeErrMSE(i) = sum(currLargeErrMSE(largeErrorInd > 0))/sum(largeErrorInd > 0);
    end
    
    ProbOfLargeErr(i) = sum(largeErrorInd)/length(largeErrorInd);
    
    if mod(i,5) == 0
        disp(strcat('Finished ENR = ',num2str(ENR(i))));
    end
end

% ZZLB
deps = 1e-5;
eps = 0:deps:10;
eps(eps == 0) = deps;
zzlb = zeros(size(ENRlin));
autoCorr = zeros(size(eps));
autoCorr(eps < 1/(beta)) = 1 - beta*eps(eps < 1/(beta));
for i=1:length(ENRlin)
    zzlb(i) = 4 * deps * sum((eps).*qfunc(eps/2).*qfunc(sqrt(ENRlin(i) * (1 - autoCorr))));
end

% WWB
deps = 1e-4;
eps = 0:deps:10;
eps(eps == 0) = deps;
wwb = zeros(size(ENRlin));optH = zeros(size(ENRlin));
for i=1:length(ENRlin)
    J = zeros(size(eps));
    autoCorr = zeros(size(eps)); autoCorr(eps < 1/(beta)) = 1 - beta*eps(eps < 1/(beta));
    autoCorr2 = zeros(size(eps)); autoCorr2(eps < 1/(2*beta)) = 1 - (2*beta)*eps(eps < 1/(2*beta));
    J = (0.5*(eps.^2).*(exp(-0.25*eps.^2)).*exp(-0.5*ENRlin(i)*(1 - autoCorr)))./(exp(-0.125*eps.^2) - exp(-0.125*(2*eps).^2).*exp(-0.25*ENRlin(i)*(1 - autoCorr2)));
    
%     autoCorr = zeros(size(eps)); autoCorr(eps < 1/(beta)) = 1 - beta*eps(eps < 1/(beta));
%     autoCorr2 = zeros(size(eps)); autoCorr2(eps < 1/(beta/2)) = 1 - (beta/2)*eps(eps < 1/(beta/2));
    
%     J = (0.25*(eps.^2).*(exp(-0.25*eps.^2)).*exp(-0.5*ENRlin(i)*(1 - autoCorr)))./(exp(-0.125*eps.^2) - exp(-0.125*(2*eps).^2).*exp(-0.25*ENRlin(i)*(1 - autoCorr2)));

    [wwb(i),minIdx] = max(J);
    optH(i) = eps(minIdx);
end

% Upper Bound
% LargeErrProb_UB = beta*(pi^(3/2))*((ENRlin).^(1/2)).*(1 + 0.75*sqrt(1./ENRlin) + 0.5./ENRlin).*exp(-ENRlin/2) + beta*sqrt(pi)*exp(-2*ENRlin);
% LargeErrProb_UB = min(LargeErrProb_UB,1);
LargeErrProb_UB = beta*exp(-ENRlin/2).*sqrt(pi^2 * ENRlin/8).*(1 + 1.85*sqrt(8*pi./ENRlin) + 1./(2*ENRlin)) + beta*sqrt(pi)*exp(-2*ENRlin);
LargeErrProb_UB = min(LargeErrProb_UB,1);

LargeErrMSE_UB = (1.5*(8*ENRlin/beta).^(2/3) + 1.35*(1/beta + sqrt(8/pi))*(8*ENRlin/beta).^(1/3) + 1 + 1/(4*beta^2) + sqrt(2/pi)/(2*beta));
D_S = (1./(beta*ENRlin).^2) .* ((13/8 + 2*sqrt(2*beta*ENRlin).*(1 - 1./(2*beta*ENRlin)).*exp(-ENRlin))./(1 - 1./(2*beta*ENRlin)).^4);

validIdx = (ProbOfLargeErr > 0);
figure; 
semilogy(ENR(validIdx),LargeErrProb_UB(validIdx),'LineWidth',2.5);hold on;
semilogy(ENR(validIdx),ProbOfLargeErr(validIdx),'-o','LineWidth',2.5);hold on;
semilogy(ENR(validIdx),beta*exp(-ENRlin(validIdx)/2).*sqrt(pi^2 * ENRlin(validIdx)/8),'-*','LineWidth',2.5);hold on;
xlabel('ENR [dB]'); ylabel('Probability Of Error');

figure;
semilogy(ENR(validIdx),LargeErrMSE_UB(validIdx),'LineWidth',2.5); hold on;
semilogy(ENR(validIdx),LargeErrMSE(validIdx),'LineWidth',2.5); hold on;
xlabel('ENR [dB]'); ylabel('Large Error MSE');

figure;hold all;
plot(ENR,10*log10(D_S.*(1 - LargeErrProb_UB) + LargeErrProb_UB.*LargeErrMSE_UB),'-','LineWidth',2.5);
plot(ENR,10*log10(MSE_MAP),'-o','LineWidth',1.5);
plot(ENR,10*log10(zzlb),'-*','LineWidth',1.5);
plot(ENR,10*log10(wwb),'-p','LineWidth',1.5);
xlabel('ENR [dB]'); ylabel('MSE [dB]');
legend('Upper Bound','Empiric','ZLLB','WWLB');


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
