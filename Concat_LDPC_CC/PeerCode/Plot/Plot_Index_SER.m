function [] = Plot_Index_SER(m,Estimated_Bhattacharyya,SER_Vec,Estimated_SER_Vec,LLR_Avg,LLR_Std,Is_Frozen_Bit_Index_Vec,title_name,legend_name,plot_path,plot_name,is_visible)
            
    N = size(SER_Vec,2);

    r = size(Is_Frozen_Bit_Index_Vec,2)/size(SER_Vec,2);
    
    Is_Frozen_Bit_Index_Vec = Is_Frozen_Bit_Index_Vec(1:r:end);
    
    %color = ['blue','green','red','magenta','cyan','yellow','black'];
    color = ['b','g','r','m','c','y'];
    
    Indexes = [1:1:N];
    SER_Plot_Handel = figure('Position', [100, 100, 1000, 800],'visible',sprintf('%s',is_visible));
    set(gca,'FontSize',15);
    set(gcf,'name',title_name,'numbertitle','off');
    
    plot(Indexes(1:N/m),SER_Vec(1:N/m),sprintf('.%s',color(1)));
    hold on;
    plot(Indexes(~Is_Frozen_Bit_Index_Vec),SER_Vec(~Is_Frozen_Bit_Index_Vec),'oblack');
    for i=2:1:m
        plot(Indexes((i-1)*(N/m)+1:i*(N/m)),SER_Vec((i-1)*(N/m)+1:i*(N/m)),sprintf('.%s',color(mod(i-1,length(color))+1)));
    end
    if(m==1)
        plot(Indexes,Estimated_SER_Vec,sprintf('.%s',color(2)));
        set(legend({legend_name,'',strcat(legend_name,' estimated')},'Location','SouthOutside'),'FontSize',15);
    else
        set(legend(legend_name,'Location','SouthOutside'),'FontSize',15);
    end
    grid on;
    set(title(title_name),'FontSize',20);
    set(xlabel('Index'),'FontSize',20);
    xlim([0,N]);
    set(ylabel('SER'),'FontSize',20);
    ylim([0,1]);
    
    try
    	saveas(SER_Plot_Handel,fullfile(plot_path,sprintf('SER_%s',plot_name)),'bmp');
        saveas(SER_Plot_Handel,fullfile(plot_path,sprintf('SER_%s',plot_name)),'fig');
    end
       
    if(strcmp(is_visible,'off'))
        close(SER_Plot_Handel)
    end
    
    Indexes = [1:1:length(LLR_Avg)];
    
    LLR_Plot_Handel = figure('Position', [100, 100, 1000, 800],'visible',sprintf('%s',is_visible));
    set(gca,'FontSize',15);
    set(gcf,'name',title_name,'numbertitle','off');
    
    plot(Indexes,LLR_Avg,'.b');
    hold on;
    plot(Indexes,LLR_Avg-4*LLR_Std,'.r');
    %plot(Indexes(~Is_Frozen_Bit_Index_Vec),LLR_Avg(~Is_Frozen_Bit_Index_Vec),'oblack');
    grid on;
    set(legend(legend_name,'Location','SouthOutside'),'FontSize',15);
    set(title(title_name),'FontSize',20);
    set(xlabel('Index'),'FontSize',20);
    set(ylabel('avg(LLR)'),'FontSize',20);
    xlim([0,N]);
    try
        ylim([0,max(LLR_Avg)]);
    end
    
    try
    	saveas(LLR_Plot_Handel,fullfile(plot_path,sprintf('Ginie_Aided_LLR_%s',plot_name)),'bmp');
        saveas(LLR_Plot_Handel,fullfile(plot_path,sprintf('Ginie_Aided_LLR_%s',plot_name)),'fig');
    end
        
    if(strcmp(is_visible,'off'))
        close(LLR_Plot_Handel)
    end
    
    Indexes = [1:1:length(Estimated_Bhattacharyya)];
    
    Bhattacharyya_Plot_Handel = figure('Position', [100, 100, 1000, 800],'visible',sprintf('%s',is_visible));
    set(gcf,'name',title_name,'numbertitle','off');
    set(gca,'FontSize',15);
    
    plot(Indexes,Estimated_Bhattacharyya,'.b');
    hold on;
    grid on;
    set(legend(legend_name,'Location','SouthOutside'),'FontSize',15);
    set(title(title_name),'FontSize',20);
    set(xlabel('Index'),'FontSize',20);
    set(ylabel('Bhattacharyya'),'FontSize',20);
    xlim([0,N]);
    try
        ylim([min(Estimated_Bhattacharyya),max(Estimated_Bhattacharyya)]);
    end
    
    try
    	saveas(Bhattacharyya_Plot_Handel,fullfile(plot_path,sprintf('Bhattacharyya_%s',plot_name)),'bmp');
        saveas(Bhattacharyya_Plot_Handel,fullfile(plot_path,sprintf('Bhattacharyya_%s',plot_name)),'fig');
    end
        
    if(strcmp(is_visible,'off'))
        close(Bhattacharyya_Plot_Handel)
    end
    
    Indexes = [1:1:length(Estimated_Bhattacharyya)];
    
    Symmetric_Capacity_Plot_Handel = figure('Position', [100, 100, 1000, 800],'visible',sprintf('%s',is_visible));
    set(gcf,'name',title_name,'numbertitle','off');
    set(gca,'FontSize',15);
    
    upper_bound = r*sqrt(1-(Estimated_Bhattacharyya/r).^2);
    plot(Indexes,upper_bound,'.b');
    hold on;
    %lower_bound = log2((2^r)./(1+Estimated_Bhattacharyya*(2^r-1)/r));
    lower_bound = 1-Estimated_Bhattacharyya;
    plot(Indexes,lower_bound,'.r');
    grid on;
    set(legend('Upper bound','lower bound','Location','SouthOutside'),'FontSize',15);
    set(title(title_name),'FontSize',20);
    set(xlabel('Index'),'FontSize',20);
    set(ylabel('Symmetric Capacity'),'FontSize',20);
    xlim([0,N]);
    try
%         ylim([min(lower_bound),max(upper_bound)]);
        ylim([0,r]);
    end
    
    try
    	saveas(Symmetric_Capacity_Plot_Handel,fullfile(plot_path,sprintf('Symmetric_Capacity_%s',plot_name)),'bmp');
        saveas(Symmetric_Capacity_Plot_Handel,fullfile(plot_path,sprintf('Symmetric_Capacity_%s',plot_name)),'fig');
    end
        
    if(strcmp(is_visible,'off'))
        close(Symmetric_Capacity_Plot_Handel)
    end
    
    Indexes = [1:1:length(Estimated_Bhattacharyya)];
    
    Bounds_Gap_Plot_Handel = figure('Position', [100, 100, 1000, 800],'visible',sprintf('%s',is_visible));
    set(gcf,'name',title_name,'numbertitle','off');
    set(gca,'FontSize',15);
    
    bound_gap = upper_bound-lower_bound;
    plot(Indexes,bound_gap,'.b');
    hold on;
    plot(Indexes(~Is_Frozen_Bit_Index_Vec),bound_gap(~Is_Frozen_Bit_Index_Vec),'oblack');
    plot(Indexes,ones(1,size(Indexes,2))*max(bound_gap(~Is_Frozen_Bit_Index_Vec)),'r');
    grid on;
    set(legend({legend_name,strcat(legend_name,' unfrozen bits'),strcat(legend_name,' max bound gap')},'Location','SouthOutside'),'FontSize',15);
    set(title(title_name),'FontSize',20);
    set(xlabel('Index'),'FontSize',20);
    set(ylabel('Upper bound - lower bound'),'FontSize',20);
    xlim([0,N]);
    try
%         ylim([0,max(bound_gap)]);
        ylim([0,r]);
    end
    
    try
    	saveas(Bounds_Gap_Plot_Handel,fullfile(plot_path,sprintf('Bounds_Gap_%s',plot_name)),'bmp');
        saveas(Bounds_Gap_Plot_Handel,fullfile(plot_path,sprintf('Bounds_Gap_%s',plot_name)),'fig');
    end
        
    if(strcmp(is_visible,'off'))
        close(Bounds_Gap_Plot_Handel)
    end
    
    Indexes = [1:1:length(Estimated_SER_Vec)];
    
    BER_Plot_Handel = figure('Position', [100, 100, 1000, 800],'visible',sprintf('%s',is_visible));
    set(gcf,'name',title_name,'numbertitle','off');
    set(gca,'FontSize',15);
    
    temp_Indexes = Indexes(1:N/m);
    temp_Estimated_SER_Vec = Estimated_SER_Vec(1:N/m)./sum(Estimated_SER_Vec(~Is_Frozen_Bit_Index_Vec));
    temp_Is_Frozen_Bit_Index_Vec = Is_Frozen_Bit_Index_Vec(1:N/m);
    plot(temp_Indexes(~temp_Is_Frozen_Bit_Index_Vec),temp_Estimated_SER_Vec(~temp_Is_Frozen_Bit_Index_Vec),sprintf('.%s',color(1)));
    hold on;
    for i=2:1:m
        temp_Indexes = Indexes((i-1)*(N/m)+1:i*(N/m));
        temp_Estimated_SER_Vec = Estimated_SER_Vec((i-1)*(N/m)+1:i*(N/m))./sum(Estimated_SER_Vec(~Is_Frozen_Bit_Index_Vec));
        temp_Is_Frozen_Bit_Index_Vec = Is_Frozen_Bit_Index_Vec((i-1)*(N/m)+1:i*(N/m));
        plot(temp_Indexes(~temp_Is_Frozen_Bit_Index_Vec),temp_Estimated_SER_Vec(~temp_Is_Frozen_Bit_Index_Vec),sprintf('.%s',color(mod(i-1,length(color))+1)));
    end
    grid on;
    set(legend(legend_name,'Location','SouthOutside'),'FontSize',15);
    set(title(title_name),'FontSize',20);
    set(xlabel('Index'),'FontSize',20);
    set(ylabel('BER'),'FontSize',20);
    xlim([0,N]);
    try
        ylim([0,max(Estimated_SER_Vec(~Is_Frozen_Bit_Index_Vec)./sum(Estimated_SER_Vec(~Is_Frozen_Bit_Index_Vec)))]);
    end
    
    try
    	saveas(BER_Plot_Handel,fullfile(plot_path,sprintf('Ginie_Aided_BER_%s',plot_name)),'bmp');
        saveas(BER_Plot_Handel,fullfile(plot_path,sprintf('Ginie_Aided_BER_%s',plot_name)),'fig');
    end
        
    if(strcmp(is_visible,'off'))
        close(BER_Plot_Handel)
    end
    
end
