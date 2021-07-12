% median check

N = 500;
iters = 1000;

medVec = zeros(1,N);
medVecBound = zeros(1,N);
medConst = zeros(1,N);
for i=1:N
    currMed = 0; currMedBound = 0;

    xVec = randn(iters,i);
    yVec = ((1:i)/i).^2;
    yVec = repmat(yVec,iters,1);
    
    maxVec1 = max(xVec + yVec,[],2);
    maxVec2 = max(xVec,[],2) ;
    medVec(i) = median(maxVec1);
    medVecBound(i) = median(maxVec2) + median(((1:i)/i).^2);
    
    medConst(i) = median(((1:i)/i).^2);
end

figure;hold all;
plot(1:N,medVec,'LineWidth',2);
plot(1:N,medVecBound,'--','LineWidth',2);
plot(1:N,sqrt(2*log(1:N)),'-.','LineWidth',2);
grid on; grid minor; xlabel('N'); ylabel('Median');
legend('Median','Sum Of Medians','sqrt(2*log(N))');

figure;plot(1:N,medConst,'LineWidth',2);
grid on; grid minor; xlabel('N'); ylabel('Median');
