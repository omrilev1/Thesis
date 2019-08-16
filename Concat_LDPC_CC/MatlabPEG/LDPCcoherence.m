% calculate the coherence of LDPC
% matrix. 
function c = LDPCcoherence(A);

% normalize the columns to have unit norm
d = nnz(A(:,1));

An = A ./ sqrt(d);

G = An'*An;
H = G - eye(size(G));
c = max(max(abs(H)));


