addpath(genpath('..\Utils'))

%% This Script Calculate the R(D) Curve for ECDQ scheme
% We calculate the R(D) For constant dither , and then average and plot
% for the random dither case

% Run Parameters : 
% Signal Variance is 1 , but the distribution can be set as input

clear ; close all; clc;

pdfType = {'Gaussian'}; % 'Exp' , 'Uniform' , 'PAM'
n = 1; % Lattice dimension
nPoints = 5;
DeltaVec = [0.05:0.05:nPoints];
variance = 1;
ditherLength = 40;

avgH = zeros(length(pdfType),length(DeltaVec));
avgDist = zeros(length(pdfType),length(DeltaVec));

currH = zeros(length(pdfType),ditherLength,length(DeltaVec));
currDist = zeros(length(pdfType),ditherLength,length(DeltaVec));

for i=1:length(pdfType)
    for j=1:length(DeltaVec)
    
        Delta = DeltaVec(j);
        
        dither = linspace(-Delta/2,Delta/2,ditherLength);

        effectiveLength = ceil(5*sqrt(variance) / Delta);
        % make number of Delta segments odd 
        if mod(effectiveLength,2) == 0
            effectiveLength = effectiveLength + 1;
        end
        
        codebook = -((effectiveLength+1)/2)*Delta : Delta : ((effectiveLength+1)/2)*Delta;
        partition = -((effectiveLength-1)/2)*Delta - Delta/2 : Delta : ((effectiveLength-1)/2)*Delta + Delta/2;
        
        [pdf,cordX,cordY,dx] = pdfGenerator(pdfType,variance,n);
        
        for k=1 : length(dither)
            % quantize distribution and calculate the probabilites : Pr{s+u in lambda}
            quantIndex = quantiz(cordX + dither(k),partition,codebook);
            segmentsVal = unique(quantIndex);
            
            currH(i,k,j) = 0;
            for valIdx = 1 : length(segmentsVal)
                quantSegment = (quantIndex == segmentsVal(valIdx));
                prob = sum(pdf(quantSegment))*dx;
                currDist(i,k,j) = currDist(i,k,j) + dx*sum(pdf(quantSegment) .* (cordX(quantSegment) + dither(k) - codebook(valIdx)).^2);
                
                if prob > 0 
                    currH(i,k,j) = currH(i,k,j) - prob*log2(prob);
                end 
            end
        end
        
        entropyForAvg = currH(i,:,j);
        distForAvg = currDist(i,:,j);
        
        avgH(i,j) = sum(entropyForAvg)/ditherLength;
        avgDist(i,j) = sum(distForAvg)/ditherLength;
    end
end

% relevant Delta Indices 
idx_1 = find(DeltaVec == 1);
idx_2 = find(DeltaVec == 2);
idx_25 = find(DeltaVec == 2.5);
idx_3 = find(DeltaVec == 3);
idx_35 = find(DeltaVec == 3.5);
idx_4 = find(DeltaVec == 4);
idx_5 = find(DeltaVec == 5);



D = 0.001:0.001:1;
Rd = 0.5*log(1./D);
figure; hold all
plot(D,Rd,'-','LineWidth',1.5)
plot(avgDist,avgH,'-','LineWidth',1.5)
plot(currDist(1,:,idx_1),currH(1,:,idx_1),'--','LineWidth',1.7)
plot(currDist(1,:,idx_2),currH(1,:,idx_2),'--','LineWidth',1.7)
plot(currDist(1,:,idx_25),currH(1,:,idx_25),'--','LineWidth',1.7)
plot(currDist(1,:,idx_3),currH(1,:,idx_3),'--','LineWidth',1.7)
plot(currDist(1,:,idx_35),currH(1,:,idx_35),'--','LineWidth',1.7)
plot(currDist(1,:,idx_4),currH(1,:,idx_4),'--','LineWidth',1.7)
plot(currDist(1,:,idx_5),currH(1,:,idx_5),'--','LineWidth',1.7)
grid on; grid minor;
xlabel('D'); ylabel('R [bits]')
title('Distortion Vs Rate : Shannon and ECDQ')
legend('Shannon','ECDQ','\Delta = 1','\Delta = 2','\Delta = 2.5','\Delta = 3','\Delta = 3.5','\Delta = 4','\Delta = 5')

xlim([0 2]); ylim([0 3.5])



        
        
            