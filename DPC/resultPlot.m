function [] =   resultPlot(snr,capacity_AWGN,capacity_DPC,capacity_FiniteDPC,inputType,SIR)

snrLin = 10.^(snr/10);
% AWGN Channel capacity
GaussianCapacity = 0.5*log2(1+ 10.^(snr/10));

% Plot Finite DPC + DPC + AWGN + Lower Bound , for uniform input , EsN0
uniform_idx = find(strcmp(inputType,'uni'));
% Bounds and awgn Capacity
figure;hold all
plot(snr,GaussianCapacity,'-*','LineWidth',2.5);
plot(snr,capacity_AWGN(uniform_idx,:),'-o','LineWidth',2.5);
plot(snr,capacity_DPC(uniform_idx,:),'-p','LineWidth',2.5);
plot(snr,capacity_FiniteDPC(uniform_idx,:),'-^','LineWidth',2.5);

set(gca,'FontSize', 12, 'FontName', 'Times New Roman');
grid on; grid minor;
box on;
xlabel('SNR (Es/N0) [dB]','FontSize', 14, 'FontName', 'Times New Roman');
ylabel('Capacity (bits/Tx)','FontSize',14,'FontName', 'Times New Roman');

legend('Gaussian Capacity', 'Gaussian Capacity Uniform Input', ...
    'DPC uniform interferer and input I(v;Y)','Finite Power DPC - Uniform Input');
title({'Capacity for AWGN , With uniform interferer and alpha_{MMSE}',strcat('SIR = ',num2str(SIR),' [dB]')});

% Plot Vs EbN0
% Bounds and awgn Capacity
figure;hold all
plot(10*log10(snrLin./(2*GaussianCapacity)),GaussianCapacity,'-*','LineWidth',2.5);
plot(10*log10(snrLin./(2*capacity_AWGN(uniform_idx,:))),capacity_AWGN(uniform_idx,:),'-o','LineWidth',2.5);
plot(10*log10(snrLin./(2*capacity_DPC(uniform_idx,:))),capacity_DPC(uniform_idx,:),'-p','LineWidth',2.5);
plot(10*log10(snrLin./(2*capacity_FiniteDPC(uniform_idx,:))),capacity_FiniteDPC(uniform_idx,:),'-^','LineWidth',2.5);

set(gca,'FontSize', 12, 'FontName', 'Times New Roman');
grid on; grid minor;
box on;
xlabel('Eb / N0 [dB]','FontSize', 14, 'FontName', 'Times New Roman');
ylabel('Capacity (bits/Tx)','FontSize',14,'FontName', 'Times New Roman');

legend('Gaussian Capacity', 'Gaussian Capacity Uniform Input', ...
    'DPC uniform interferer and input I(v;Y)','Finite Power DPC - Uniform Input');
title({'Capacity for AWGN , With uniform interferer and alpha_{MMSE}',strcat('SIR = ',num2str(SIR),' [dB]')});
xlim([0 6]); ylim([0 2]);
end

