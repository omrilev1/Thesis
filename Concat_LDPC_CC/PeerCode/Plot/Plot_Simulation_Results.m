function [] = Plot_Simulation_Results(m,SNR_Type,SNR_Vec_dB,BLER,Min_BLER,BER,First_Symbol,LLR_Avg,LLR_Std,Is_Frozen_Bit_Index_Vec,title_name,legend_name,plot_path,plot_name,is_visible)

    N = size(First_Symbol,2);

    r = size(Is_Frozen_Bit_Index_Vec,2)/size(First_Symbol,2);
    
    Is_Frozen_Bit_Index_Vec = Is_Frozen_Bit_Index_Vec(1:r:end);
    
    Max_BLER = 10^0;
    Max_Min_BLER = 10^-1;
    SNR_Axis_Buffer = 0.5;
    
    %color = ['blue','green','red','magenta','cyan','yellow','black'];
    color = ['b','g','r','m','c','y'];
    
    %BLER_Plot_Handel = figure('Position', [100, 100, 1000, 800],'visible',sprintf('%s',is_visible));
    BLER_Plot_Handel = figure('Position', [100, 100, 1000, 800]);
    
    BLER = max(BLER,Min_BLER*ones(1,length(BLER)));
    
    set(gcf,'name',title_name,'numbertitle','off');
    set(gca,'FontSize',15);
    
    semilogy(SNR_Vec_dB,BLER,'-bo','MarkerSize',6,'LineWidth',1.5);
    grid on;
    hold on;
    set(legend(legend_name,'Location','SouthOutside'),'FontSize',15);
    set(title(title_name),'FontSize',20);
    set(xlabel(sprintf('SNR(%s)[dB]',SNR_Type)),'FontSize',20);
    xlim([min(SNR_Vec_dB)-SNR_Axis_Buffer,max(SNR_Vec_dB)+SNR_Axis_Buffer]);
    set(ylabel('BLER'),'FontSize',20);
    
    try
        ylim([10^floor(log10(min(min(BLER),Max_Min_BLER))),Max_BLER]);
    end
    
    try
        saveas(BLER_Plot_Handel,fullfile(plot_path,sprintf('BLER_%s',plot_name)),'bmp');
        saveas(BLER_Plot_Handel,fullfile(plot_path,sprintf('BLER_%s',plot_name)),'fig');
    end
    
    Max_BER = 10^0;
    Max_Min_BER = 10^-1;
    SNR_Axis_Buffer = 0.5;

    BER_Plot_Handel = figure('Position', [100, 100, 1000, 800],'visible',sprintf('%s',is_visible));
        
    set(gcf,'name',title_name,'numbertitle','off');
    set(gca,'FontSize',15);
    
    semilogy(SNR_Vec_dB,BER,'-bo','MarkerSize',6,'LineWidth',1.5);
    grid on;
    hold on;
    set(legend(legend_name,'Location','SouthOutside'),'FontSize',15);
    set(title(title_name),'FontSize',20);
    set(xlabel(sprintf('SNR(%s)[dB]',SNR_Type)),'FontSize',20);
    xlim([min(SNR_Vec_dB)-SNR_Axis_Buffer,max(SNR_Vec_dB)+SNR_Axis_Buffer]);
    set(ylabel('BER'),'FontSize',20);
    
    try
        ylim([10^floor(log10(min(min(BER),Max_Min_BER))),Max_BER]);
    end
    
    try
        saveas(BER_Plot_Handel,fullfile(plot_path,sprintf('BER_%s',plot_name)),'bmp');
        saveas(BER_Plot_Handel,fullfile(plot_path,sprintf('BER_%s',plot_name)),'fig');
    end
        
    if(strcmp(is_visible,'off'))
        close(BER_Plot_Handel)
    end
    
    Indexes = [1:1:length(First_Symbol)];
    
    First_Symbol_Plot_Handel = figure('Position', [100, 100, 1000, 800],'visible',sprintf('%s',is_visible));
        
    set(gcf,'name',title_name,'numbertitle','off');
    set(gca,'FontSize',15);
    
    temp_Indexes = Indexes(1:N/m);
    temp_First_Symbol = First_Symbol(1:N/m)./sum(First_Symbol);
    temp_Is_Frozen_Bit_Index_Vec = Is_Frozen_Bit_Index_Vec(1:N/m);
    plot(temp_Indexes(~temp_Is_Frozen_Bit_Index_Vec),temp_First_Symbol(~temp_Is_Frozen_Bit_Index_Vec),sprintf('.%s',color(1)));
    hold on;
    for i=2:1:m
        temp_Indexes = Indexes((i-1)*(N/m)+1:i*(N/m));
        temp_First_Symbol = First_Symbol((i-1)*(N/m)+1:i*(N/m))./sum(First_Symbol);
        temp_Is_Frozen_Bit_Index_Vec = Is_Frozen_Bit_Index_Vec((i-1)*(N/m)+1:i*(N/m));
        plot(temp_Indexes(~temp_Is_Frozen_Bit_Index_Vec),temp_First_Symbol(~temp_Is_Frozen_Bit_Index_Vec),sprintf('.%s',color(mod(i-1,length(color))+1)));
    end
    
    grid on;
    hold on;        
    set(legend(legend_name,'Location','SouthOutside'),'FontSize',15);
    set(title(title_name),'FontSize',20);
    set(xlabel('Index'),'FontSize',20);
    set(ylabel('pdf(first symbol error)'),'FontSize',20);
    xlim([0,N]);
    
    try
        ylim([0,max(First_Symbol(~Is_Frozen_Bit_Index_Vec)./sum(First_Symbol))]);
    end
    
    try
        saveas(First_Symbol_Plot_Handel,fullfile(plot_path,sprintf('First_Symbol_%s',plot_name)),'bmp');
        saveas(First_Symbol_Plot_Handel,fullfile(plot_path,sprintf('First_Symbol_%s',plot_name)),'fig');
    end
        
    if(strcmp(is_visible,'off'))
        close(First_Symbol_Plot_Handel)
    end
        
    Indexes = [1:1:length(LLR_Avg)];

    LLR_Plot_Handel = figure('Position', [100, 100, 1000, 800],'visible',sprintf('%s',is_visible));

    set(gcf,'name',title_name,'numbertitle','off');
    set(gca,'FontSize',15);
    
    plot(Indexes(~Is_Frozen_Bit_Index_Vec),LLR_Avg(~Is_Frozen_Bit_Index_Vec),'.b');
    hold on;
    plot(Indexes(~Is_Frozen_Bit_Index_Vec),LLR_Avg(~Is_Frozen_Bit_Index_Vec)-4*LLR_Std(~Is_Frozen_Bit_Index_Vec),'.r');
    plot(Indexes(First_Symbol~=zeros(1,length(First_Symbol))),LLR_Avg(First_Symbol~=zeros(1,length(First_Symbol))),'oblack');
    grid on;
    set(legend(legend_name,'Location','SouthOutside'),'FontSize',15);
    set(title(title_name),'FontSize',20);
    set(xlabel('Index'),'FontSize',20);
    set(ylabel('avg(LLR)'),'FontSize',20);
    xlim([0,N]);
    try
        %ylim([min(0,min(LLR_Avg(~Is_Frozen_Bit_Index_Vec)-4*LLR_Std(~Is_Frozen_Bit_Index_Vec))),max(LLR_Avg(and(~Is_Frozen_Bit_Index_Vec,~isnan(LLR_Avg))))]);
        ylim([0,max(LLR_Avg(and(~Is_Frozen_Bit_Index_Vec,~isnan(LLR_Avg))))]);
    end

    try
        saveas(LLR_Plot_Handel,fullfile(plot_path,sprintf('LLR_%s',plot_name)),'bmp');
        saveas(LLR_Plot_Handel,fullfile(plot_path,sprintf('LLR_%s',plot_name)),'fig');
    end
        
    if(strcmp(is_visible,'off'))
        close(LLR_Plot_Handel)
    end
    
end

