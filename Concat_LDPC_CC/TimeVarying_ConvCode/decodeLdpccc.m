function [vHat,finalL] = decodeLdpccc(rx, baseHT, T, N0, iteration)
% Decode LDPC-CC 
%
%  rx        : Received symbol vector
%  baseHT    : Base of parity check matrix H transpose
%  T         : Convolutional code period
%  N0        : Noise variance
%  iteration : Number of iteration
%
%  vHat      : Decoded vector (0/1)
%  finalL    : Final LLR
% Bagawan S. Nugroho 2007

% Convolutional code memory
Ms = T + 1;

% Received vector length
dataLen = length(rx);

% number of repetitions of basic Conv Code
numOfChuncks = dataLen/(2*(Ms + 1));

% ----------- Initialization ----------

% Log likelihood.
L = (-4/N0).*reshape(rx, 2*(Ms + 1), dataLen/(2*(Ms + 1)));

finalL = zeros(2*(Ms + 1),dataLen/(2*(Ms + 1)));
newL = zeros(2*(Ms + 1),2*(Ms + 1) - 1,dataLen/(2*(Ms + 1)));

% H transpose starting point
offsetV = 0;

for k = 1:numOfChuncks - 1

   % Create H transpose
   [upperH,lowerH,nextOffsetV] = parityCheckMatrixHT(T, baseHT, offsetV);
   offsetV = nextOffsetV;
   [M,N] = size(lowerH);
   
   % Compose the H transpose
   HT(:, :, k) = [upperH; lowerH];
   
   % Log-likelihood matrix
   Lqij(:, :, k) = lowerH.*repmat([L(:, k); L(1:end - 2, k + 1)], 1, N);
        
end % for k

% The last sub statistics 
k = k + 1;

% Create H transpose
[upperH,lowerH,nextOffsetV] = parityCheckMatrixHT(T, baseHT, offsetV);

% Compose the H transpose
HT(:, :, k) = [upperH; lowerH];

% Log-likelihood matrix
Lqij(:, :, k) = lowerH.*repmat([L(:, k); zeros(2*(Ms + 1) - 2, 1)], 1, N);


% ---------- Min-sum decoding ----------

% Iteration
for i = 1:iteration
   
   % Initial Lrji value (likelihood for zeros)
   prevL = 1e10*upperH(:, :, 1);

   for j = 1:numOfChuncks
        
      [currL,newL(:, :, j),finalL(:, j)] = minSum_mex(L(:, j), prevL, Lqij(:, :, j), HT(:, :, j), Ms); 
      prevL = currL;
      
   end % for j
   
   % Update the current and future statistics
   [Lqij,~] = statUpdate(numOfChuncks,Ms,newL,Lqij);
   
   % The last sub statistics
   k = numOfChuncks;
   
   % Current statistics
   Lqij(:, :, k) = [newL(:, :, k); zeros(2*(Ms + 1) - 2, 2*(Ms + 1) - 1)];
   
end % for i


% ------------- Decision ---------------
vHat = zeros(2*(Ms + 1),numOfChuncks - 1);
for i = 1:numOfChuncks - 1
   
   for j = 1:2*(Ms + 1)
   
      % Hard decoding finalL
      if finalL(j, i) < 0
         vHat(j, i) = 1;
      else
         vHat(j, i) = 0;
      end
      
   end % for j

end % for i
