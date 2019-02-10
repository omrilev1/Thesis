function [indices] = resultPlot(deltaToPlot,deltaVec,pdfToPlot,pdfVec,H,D,ditherLength)

% Calculate the Avg over the dither
avgH = reshape(mean(H,2),size(H,1),[]);
avgDist = reshape(mean(D,2),size(D,1),[]);

% find Delta's indices in the delta vector
indices = zeros(1,length(deltaToPlot));
for i=1:length(deltaToPlot)
    [~,indices(i)] = min(abs(deltaVec - deltaToPlot(i)));
end

% Figure titles
pdfIdx = zeros(1,length(pdfToPlot));
currTitle = [];
for i=1:length(pdfToPlot)
    for j=1:length(pdfVec)
        count = 1;
        if strcmp(pdfVec(j),pdfToPlot(i))
            pdfIdx(i) = count;
            break;
        else
            count = count + 1;
        end
    end
    currTitle = [currTitle pdfToPlot(i)];
end

% iterate over the distributions and plot the resulting curves
for i=1:length(pdfToPlot)
    
    if strcmp(pdfToPlot(i),'Gaussian')
        % R(D) Shannon curve
        shannonD = 0.001:0.001:1;
        Rd = max(0.5*log(0.25./shannonD),0);
    else
        % R(D) Shannon curve
        shannonD = 0.001:0.001:1;
        lambda = 1/0.25;
        Rd = max(-0.5*log(lambda*shannonD),0);
    end
    % convex hull of R(D|U)
    currH = reshape(H(i,:,:),1,[]);
    currD = reshape(D(i,:,:),1,[]);
    convHull_Idx = convhull(currD,currH);
    
    % take only the points on the relevant region
    [~,relevant_ConvHull_Idx] = find(currD(convHull_Idx) < 4);
    convHull_Idx_ToPlot = convHull_Idx(relevant_ConvHull_Idx);
    
    % Plot
    figure; hold all
    plot(shannonD,max(Rd,0),'-','LineWidth',1.5)
%     plot(currD(convHull_Idx_ToPlot),currH(convHull_Idx_ToPlot)/2,'-cp')
    plot(avgDist(i,:),avgH(i,:)/2,'-','LineWidth',2)
%     currLegend = {'Shannon : R(D) = 0.5*log(1/D)','Optimum Dither Convex Hull','ECDQ'};
    currLegend = {'R(D)','ECDQ'};
    
    curveStyle = ['ro';'ks';'cp';'mo';'gs'];
    for j=1:length(deltaToPlot)
        plot(D(i,:,indices(j)),H(i,:,indices(j))/2,curveStyle(j,:),'LineWidth',1.7)
        currLegend = [currLegend strcat('\Delta =',num2str(deltaToPlot(j)))];
    end
    grid on; grid minor;
    xlabel('D'); ylabel('R [bits]')
    legend(currLegend);
    title(currTitle(i));
    
%     xlim([0 2]); ylim([0 3.5])
    
    % Plot the optimal Delta of the points on the convex hull
    [I,J] = ind2sub([size(H,2) size(H,3)],convHull_Idx_ToPlot);
    optimalDither = zeros(2,length(I));
    for k = 1 : length(I)
        G = [0 sqrt(3); 2 1];
        basic_latticePoints = latticeGen(G,50);
        currDither = generate2D_dither(deltaVec(J(k)),ditherLength,2,deltaVec(J(k))*basic_latticePoints);
        optimalDither(:,k) =  currDither(:,I(k))/deltaVec(J(k));
    end
    figure;
    plot3(optimalDither(1,:),optimalDither(2,:),deltaVec(J),'gp','LineWidth',1.5);
    grid on; grid minor;
    xlabel('\Delta'); ylabel('dither/\Delta');
    title(strcat('optimal Dither Value - ',currTitle(i)));
    
    
end
end