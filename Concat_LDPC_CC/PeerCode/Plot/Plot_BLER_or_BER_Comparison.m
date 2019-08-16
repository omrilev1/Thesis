function [] = Plot_BLER_or_BER_Comparison(Is_BLER,SNR_Type,SNR_Vec_dB,BLER_or_BER_Vec,Num_Of_Colors,Min_BLER_or_BER,title_name,legend_name,plot_path,plot_name)
    
    Max_BLER_or_BER = 10^0;
    Max_Min_BLER_or_BER = 10^-1;
    SNR_Axis_Buffer = 0.5;
    
    %color = ['blue','green','red','magenta','cyan','yellow','black'];
    color = ['b','g','r','m','c','y'];
    shape = ['o','s','^','*','d','+','x','v','>','<','p','h','.'];
    
    Num_Of_Colors = min(Num_Of_Colors,length(color));
    
    BLER_or_BER_Plot_Handel = figure();
    
    set(gcf,'name',title_name,'numbertitle','off');
    set(gca,'FontSize',15);

    BLER_or_BER_Vec(1,:) = max(BLER_or_BER_Vec(1,:),Min_BLER_or_BER*ones(1,length(BLER_or_BER_Vec(1,:))));
    
    semilogy(SNR_Vec_dB,BLER_or_BER_Vec(1,:),sprintf('-%s%s',color(1),shape(1)),'MarkerSize',6,'LineWidth',1.5);
    hold on;
    for i=2:1:length(BLER_or_BER_Vec(:,1))
        BLER_or_BER_Vec(i,:) = max(BLER_or_BER_Vec(i,:),Min_BLER_or_BER*ones(1,length(BLER_or_BER_Vec(i,:))));
        semilogy(SNR_Vec_dB,BLER_or_BER_Vec(i,:),sprintf('-%s%s',color(mod(i-1,Num_Of_Colors)+1),shape(floor((i-1)/Num_Of_Colors)+1)),'MarkerSize',6,'LineWidth',1.5);
    end

    grid on;
    set(legend(legend_name,'Location','best'),'FontSize',15);
    set(title(title_name),'FontSize',12);
    set(xlabel(sprintf('SNR(%s)[dB]',SNR_Type)),'FontSize',12);
    xlim([min(SNR_Vec_dB)-SNR_Axis_Buffer,max(SNR_Vec_dB)+SNR_Axis_Buffer]);
    if(Is_BLER)
        set(ylabel('BLER'),'FontSize',12);
    else
        set(ylabel('BER'),'FontSize',12);
    end
    ylim([10^floor(log10(min(min(min(BLER_or_BER_Vec)),Max_Min_BLER_or_BER))),Max_BLER_or_BER]);

%     saveas(BLER_or_BER_Plot_Handel,fullfile(plot_path,plot_name),'bmp');
%     saveas(BLER_or_BER_Plot_Handel,fullfile(plot_path,plot_name),'fig');
    
end

