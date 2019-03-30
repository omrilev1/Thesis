
snr = [25 20 15];
dx = 1e-3;
x = -12:dx:12;
P_x = 1;
delta = sqrt(3*P_x);

origVec = zeros(size(x));
curr_uni = zeros(size(x));
curr_uni(abs(x) <= delta) = 1/(2*delta);

figure; hold all
for i=1:length(snr)
    sigma_z = 10^(-snr(i)/20);
    z = 1/(sqrt(2*pi)*sigma_z) * exp(-0.5*x.^2/(sigma_z^2));
    
    f = conv(curr_uni,z)*dx;
    f_cut = f(length(curr_uni)/2:end-length(curr_uni)/2);
    
    CDF = cumsum(f_cut)*dx;
    plot(x,log10(1 - CDF),'LineWidth',2)
    plot(x,log10(qfunc(x/sqrt(P_x/(1+10^(snr(i)/10))))),'LineWidth',2)
end
grid on; grid minor;
legend('SNR = 10[dB]','Gaussian Approx','SNR = 8[dB]','Gaussian Approx','SNR = 5[dB]','Gaussian Approx')