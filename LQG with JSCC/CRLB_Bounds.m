
% Calculate and print some CRLB bounds for different channel transforms
% We concentrate on the case of JSCC with zero delay, with side information
% available only @ the decoder
% According to Busgang theorem we have 
% f(s) = k*s + eta , where eta in not correlated with s. 
% distribution type referes to distribution of eta. 

close all; clear all; clc;

dist = {'Gaussian';'Uniform';'laplace'}; 

dz = 0.005;
z = -5:dz:5;
s = [-3:dz:(0-10*dz) (0+10*dz):dz:3];
k = 0.01:0.005:0.99;
snr = 6:1:15;
snrLin = 10.^(snr/10);

rho = 0.9; 

GaussianBound = zeros(size(snr));
UniformBound = zeros(size(snr));
LaplaceBound = zeros(size(snr));

for i=1:length(snrLin)
    
    GaussianBound(i) = (1-rho^2)/(1 + snrLin(i) * (1-rho^2));
    
    currUniBound = zeros(size(k));
    currLaplaceBound = zeros(size(k));
    for j = 1:length(k)
        
        delta = sqrt(3*(1 - k(j)^2));
        f_uni = (1/(2*delta)) * (qfunc((z-delta)*sqrt(snrLin(i))) - qfunc((z+delta)*sqrt(snrLin(i))));
        
         f_dot_uni = sqrt(snrLin(i)/(2*pi)) * (1/(2*delta))*(exp(-0.5*snrLin(i)*(z+delta).^2) - ...
             exp(-0.5*snrLin(i)*(z-delta).^2));
        
        lambda = sqrt(2/(1-k(j)^2));
        f_laplace = 0.5*lambda * exp(0.5*lambda^2 / snrLin(i)) * (exp(-lambda*z).*qfunc(sqrt(snrLin(i))*(lambda/snrLin(i) - z)) + ...
            exp(lambda*z).*qfunc(sqrt(snrLin(i))*(lambda/snrLin(i) + z)));
        
        f_dot_laplace = 0.5*lambda * exp(0.5*lambda^2 / snrLin(i)) * ...
            (sqrt(snrLin(i)/(2*pi)) * exp(lambda*z).*exp(-0.5*snrLin(i)*(lambda/snrLin(i) + z).^2) - ...
            sqrt(snrLin(i)/(2*pi)) * exp(-1*lambda*z).*exp(-0.5*snrLin(i)*(lambda/snrLin(i) - z).^2) + ...
            lambda*exp(lambda*z).*qfunc(sqrt(snrLin(i))*(lambda/snrLin(i) + z)) - ...
            lambda*exp(-1*lambda*z).*qfunc(sqrt(snrLin(i))*(lambda/snrLin(i) - z)));
        
        validIdx = abs(f_uni) > 0;
        currUniBound(j) = (k(j))^2 * sum(dz * (f_dot_uni(validIdx).^2)./f_uni(validIdx));
        
        validIdx = abs(f_laplace) > 0;
        currLaplaceBound(j) = (k(j))^2 * sum(dz * (f_dot_laplace(validIdx).^2)./f_laplace(validIdx));
    end
    UniformBound(i) = 1/(1/(1-rho^2) + max(currUniBound));
    LaplaceBound(i) = 1/(1/(1-rho^2) + max(currLaplaceBound));
end

figure;hold all
plot(snr , 10*log10((snrLin/(1 - rho^2) + 1)),'LineWidth',3)
plot(snr , 10*log10(1./GaussianBound),'g.-','LineWidth',3)
plot(snr , 10*log10(1./UniformBound),'-kd','LineWidth',3)
plot(snr , 10*log10(1./LaplaceBound),'r--','LineWidth',3)
grid on; grid minor;
xlabel('SNR [dB]'); ylabel('SDR [dB]');
legend('SDR OPTA - SI case : 1 + SNR/(1-\rho^2)',...
    '1/CRLB For Gaussian input noise',...
    '1/CRLB For Uniform input noise',...
    '1/CRLB For Laplace input noise');
title({strcat('Distortion Vs SNR, \rho = ',num2str(rho));'\eta independent with s'});


GaussianBound = zeros(size(snr));
UniformBound = zeros(size(snr));
LaplaceBound = zeros(size(snr));
eps = 1e-4;

f_s = (1/sqrt(2*pi)) * exp(-0.5 * s(:).^2);
z = repmat(z,length(s),1);
s = repmat(s(:),1,size(z,2));

for i=1:length(snrLin)
    
    GaussianBound(i) = (1-rho^2)/(1 + snrLin(i) * (1-rho^2));
    
    currUniBound = zeros(size(k));
    currLaplaceBound = zeros(size(k));
    for j = 1:length(k)
        
        delta = sqrt(3*(1 - k(j)^2)) * s;
        f_uni = (1./(2*delta)) .* (qfunc((z-delta)*sqrt(snrLin(i))) - qfunc((z+delta)*sqrt(snrLin(i))));
        
         f_dot_uni = sqrt(snrLin(i)/(2*pi)) * (1./(2*delta)).*(exp(-0.5*snrLin(i)*(z+delta).^2) - ...
             exp(-0.5*snrLin(i)*(z-delta).^2));
        
%         lambda = sqrt(2/(1-k(j)^2)) * (1./s);
%         f_laplace = 0.5*lambda .* exp(0.5*lambda.^2 / snrLin(i)) .* (exp(-lambda.*z).*qfunc(sqrt(snrLin(i))*(lambda./snrLin(i) - z)) + ...
%             exp(lambda.*z).*qfunc(sqrt(snrLin(i))*(lambda./snrLin(i) + z)));
%         
%         f_dot_laplace = 0.5*lambda .* exp(0.5*lambda.^2 / snrLin(i)) .* ...
%             (sqrt(snrLin(i)/(2*pi)) * exp(lambda.*z).*exp(-0.5*snrLin(i)*(lambda./snrLin(i) + z).^2) - ...
%             sqrt(snrLin(i)/(2*pi)) * exp(-1*lambda.*z).*exp(-0.5*snrLin(i)*(lambda./snrLin(i) - z).^2) + ...
%             lambda.*exp(lambda.*z).*qfunc(sqrt(snrLin(i))*(lambda./snrLin(i) + z)) - ...
%             lambda.*exp(-1*lambda.*z).*qfunc(sqrt(snrLin(i))*(lambda./snrLin(i) - z)));
        
        div = f_dot_uni.^2./f_uni;
        div(abs(f_uni) < eps) = 0;
        tempUniBound = (k(j))^2 * sum(dz * div,2);
        currUniBound(j) = sum(f_s .* tempUniBound) * dz;
        
        
%         div = f_dot_laplace.^2./f_laplace;
%         div(abs(f_laplace) < eps) = 0;
%         tempLaplaceBound = (k(j))^2 * sum(dz * div,2);
%         currLaplaceBound(j) = sum(f_s .* tempLaplaceBound) * dz;
    end
    UniformBound(i) = 1/(1/(1-rho^2) + max(currUniBound));
%     LaplaceBound(i) = 1/(1/(1-rho^2) + max(currLaplaceBound));
end

figure;hold all
plot(snr , 10*log10((snrLin/(1 - rho^2) + 1)),'LineWidth',3)
plot(snr , 10*log10(1./GaussianBound),'g.-','LineWidth',3)
plot(snr , 10*log10(1./UniformBound),'-kd','LineWidth',3)
plot(snr , 10*log10(1./LaplaceBound),'r--','LineWidth',3)
grid on; grid minor;
xlabel('SNR [dB]'); ylabel('SDR [dB]');
legend('SDR OPTA - SI case : 1 + SNR/(1-\rho^2)',...
    '1/CRLB For Gaussian input noise',...
    '1/CRLB For Uniform input noise',...
    '1/CRLB For Laplace input noise');
title({strcat('Distortion Vs SNR, \rho = ',num2str(rho));'\eta depends with s'});



        
        
        