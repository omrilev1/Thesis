function [] = Plot_LLR_PDF(m,LLR_Bin_Border_Vec,LLR_PDF_Vec,title_name,legend_name,plot_path,plot_name,is_visible)

    color = ['b','g','r','m','c','y'];
    
    LLR_PDF_Plot_Handel = figure('units','normalized','outerposition',[0 0 1 1],'visible',sprintf('%s',is_visible));
        
    set(gcf,'name',title_name,'numbertitle','off');
    set(gca,'FontSize',15);
    
    temp_LLR_Bin_Border_Vec = zeros(1,size(LLR_Bin_Border_Vec,3));
    temp_LLR_PDF_Vec = zeros(1,size(LLR_PDF_Vec,3));
    
    min_Non_Zero = inf*ones(1,m);
    max_Non_Zero = -inf*ones(1,m);
    
    for i=1:1:size(LLR_Bin_Border_Vec,1)
        
        for j=1:1:m
                        
            temp_Non_Zero_Index = find(LLR_PDF_Vec(i,j,:));
            
            min_Non_Zero(1,j) = min(min_Non_Zero(1,j),LLR_Bin_Border_Vec(i,j,temp_Non_Zero_Index(1)));
            max_Non_Zero(1,j) = max(max_Non_Zero(1,j),LLR_Bin_Border_Vec(i,j,temp_Non_Zero_Index(end)));
            
        end
        
    end
    
    for i=1:1:size(LLR_Bin_Border_Vec,1)
        
        for j=1:1:m
            
            subplot(size(LLR_Bin_Border_Vec,1),m,(i-1)*m+j);
            temp_LLR_Bin_Border_Vec(:) = LLR_Bin_Border_Vec(i,j,:);
            temp_LLR_PDF_Vec(:) = LLR_PDF_Vec(i,j,:);
            bar(temp_LLR_Bin_Border_Vec,temp_LLR_PDF_Vec,sprintf('%s',color(j)));
            set(title(sprintf('%s %dbit BER=%.3f',legend_name{i},j,sum(temp_LLR_PDF_Vec(temp_LLR_Bin_Border_Vec<0)))),'FontSize',16);
            set(xlabel('LLR'),'FontSize',20);
            xlim([min_Non_Zero(1,j),max_Non_Zero(1,j)]);
            set(ylabel('PDF'),'FontSize',20);
            
        end
        
    end
    
    try
    	saveas(LLR_PDF_Plot_Handel,fullfile(plot_path,sprintf('LLR_PDF_%s',plot_name)),'bmp');
        saveas(LLR_PDF_Plot_Handel,fullfile(plot_path,sprintf('LLR_PDF_%s',plot_name)),'fig');
    end
        
    if(strcmp(is_visible,'off'))
        close(LLR_PDF_Plot_Handel)
    end
    
end
