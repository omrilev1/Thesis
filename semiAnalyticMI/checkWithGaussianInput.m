dX = 0.01;
x = -50:dX:50;

Px = 1;
SNR = -5:1:15;
SNRlin = 10.^(SNR/10);

GaussianInputPDF = dX*(1/sqrt(2*pi))*exp(-x.^2/(2));
MI = zeros(1,length(SNR));
for i=1:length(SNRlin)
    sigma = 1/SNRlin(i);
    GaussianNoisePDF = dX*(1/sqrt(2*pi*sigma^2))*exp(-x.^2/(2*sigma^2));
    
    % calculate pdf(x,y) and pdf(Y)
    [ pdfY,intervalY ] = pdfOfSum(GaussianInputPDF,x,GaussianInputPDF,x);
    [ pdfXY,intervalXY ] = pdfAfterConstAddition( GaussianNoisePDF,x,x);
    
    % calculate H(y) and H(Y|X)
    Hy = numeric_Entropy_integral(pdfY.',intervalY.');
    Hy_given_x = numeric_Entropy_integral(pdfXY.',intervalXY);
    Hxy = sum(Hy_given_x.*GaussianInputPDF)*dX;
    
    MI(i) = Hy - Hxy;
end
awgnCapTheory = 0.5*log2(1 + SNRlin);
figure;
plot(SNR,awgnCapTheory,'--','LineWidth',1.5);hold on;
plot(SNR,MI,'-*','LineWidth',1.5);hold on;
grid on; grid minor;
xlabel('SNR [dB]'); ylabel('I(X;Y)');
title('Capacity')
legend('Theoretic AWGN','Simulation AWGN');
    
    