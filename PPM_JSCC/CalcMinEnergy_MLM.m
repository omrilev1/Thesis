% minimum energy for Infinite BW JSCC with MLM

alpha = 1.1:0.1:10;
Delta = 0.01:0.01:1;

i = 1:1:1e5;
[DeltaGrid,AlphaGrid] = meshgrid(Delta,alpha);
0.5*min(min(exp(AlphaGrid)./DeltaGrid + DeltaGrid.*(exp(AlphaGrid) - exp(-1*AlphaGrid)).*exp(-2*AlphaGrid)./(1 - exp(-1*AlphaGrid))))

% L = 2 with PPM 
k_G = 0.263; k_U = 13.76;

N1 = 0.0901:0.0901:3;
minEnergy = 0.5./N1 + N1.*Lambert_W(((1./k_G)*((N1.^2)./(1 + N1.^2))).^3);
figure;semilogy(N1,minEnergy)
min(minEnergy)

N1 = 0.0901:0.0901:80;
minEnergy = 0.5./N1 + N1.*Lambert_W(((1./k_U)*((N1.^2)./(1 + N1.^2))).^3);
figure;semilogy(N1,minEnergy)
min(minEnergy)