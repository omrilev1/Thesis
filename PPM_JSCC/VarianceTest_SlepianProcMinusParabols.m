%% Analog PPM Performance simulation
close all; clear all; clc;

StartPoint = 0;
gamma = 1; 

dt = 1/(100); Fs = 1/dt;
t = -25:dt:25;
tIdx = find(abs(t) <= 0.5);

MeanMax = zeros(size(gamma));% Mean = zeros(size(ENR));
VarMax = zeros(size(gamma));% Mean = zeros(size(ENR));

Nrun = 1e5;

ppmPulse = zeros(size(t));
ppmPulse(abs(t) < 1/(2)) = 1;
ppmPulse = ppmPulse / sqrt(sum(abs(ppmPulse.^2)*dt));
% ppmPulse(abs(t) > 0.5) = 0;

Lx = length(ppmPulse);
Ly = length(ppmPulse)+length(ppmPulse)-1;
Ly2 = pow2(nextpow2(Ly));
PPMfreq = fft(ppmPulse, Ly2);


for i=1:length(gamma)
    currMax = zeros(1,Nrun);
    currMaxSquared = zeros(1,Nrun);

    parfor n=1:Nrun
        
        tForMax = t;
        noise = randn(size(t));
        noise = sqrt(1/(2*dt))*noise;% sqrt(2*Fs/W)*noise;
        corr = fconv(noise,Lx,Ly,Ly2,PPMfreq);
        
        SlepMinusParabola = corr*dt - 0.5*gamma(i)*(t.^2);
        [~,maxIdx] = max(SlepMinusParabola);
        xHat = tForMax(maxIdx);
        
        currMax(n) = xHat;
        currMaxSquared(n) = xHat^2;
    end
    
    MeanMax(i) = sum(currMax)/Nrun;
    VarMax(i) = sum(currMaxSquared)/Nrun;

end

figure;plot(gamma,VarMax,'-','LineWidth',2);
hold on; plot(gamma,(StartPoint)^2 + (gamma).^(-4/3),'-','LineWidth',2);
grid on; grid minor; 
legend('Variance of Maximum','\gamma^{-4/3}');
xlabel('\gamma'); ylabel('Second Moment of Max Index');



% figure;plot(gamma,VarMax,'-o','LineWidth',2);
% hold on; plot(StartPoint,StartPoint.^2,'-*','LineWidth',2);
% grid on; grid minor; 
% legend('Variance of Maximum','Squared of Start Point');
% xlabel('Start Point'); ylabel('Second Moment of Max Index');
% 
% 
% figure;plot(StartPoint,StartPoint.^2 - VarMax,'LineWidth',2);
% grid on; grid minor; 
% xlabel('Start Point'); ylabel('Diff');

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
