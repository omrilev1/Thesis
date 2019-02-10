function [U,T,V] = ggmd(A,N)
% [U,T,V] = ggmd(A,N)
% Matlab implementation of the "K-user Geometric Mean Decomposition".
% Version of Idan Livni, February 3, 2012.
% Copyright 2012, Tel Aviv University, Tel Aviv, Israel.
%
% The function calculates the joint K-GMD of the given matrices.
% such that U{i}'*A{i}*V = T{i},
% where:
%	U{i} and V have orthonormal.
% 	T{i} are triangular with constant diagonals.
% @param A - cell array of K square matrices of the same dimension nxn to be jointly decomposed.
% @param N - number of duplications performed for each matrix, such that N > n^(K-1).
% @output V - right orthonormal-column matrix of the decomposition.
% @output U - cell array of size K with the left orthonormal-column matrices of the decomposition.
% @output T - cell array of size K with the triangular matrices of the decomposition.

K = length(A);
[m,n] = size(A{1});

if (N<n^(K-1))
    ERROR('ggmd: N should be at least n^(K-1)')
end

% 1-gmd of the first matrix
[u_t, t_t, v] =  one_gmd(A{1});

V = kron(eye(N), v);
for i = 1:1:K
    [u{i}, t{i}] = my_qr(A{i}*v);
    U{i} = kron(eye(N), u{i});
    T{i} = kron(eye(N), t{i});
end

% init paramaters
% nn - size of effective extended matrix.
[m,nn]=size(T{1});

% Treat remaining matrices sequentially.
for i = 1:1:K-1
    % distance between consecutive elements of the same extracted matrix.
    delta = n^(K-i)-1;
    
    % selection of columns for which current stage is applied.
    places=[];
    j=0;
    while (j*n+n+(n-1)*delta)< nn
        places=[places, [j*n+n : delta : j*n+n+(n-1)*delta]];
        j=j+1;
    end
    
    % choose selected columns in v
    v= eye(nn);
    v=v(:,places);
    
    % calculate T{i}, U{i} and V accumulated (multiplies) until current stage (not including),
    % after selection of relevant columns
    for l = 1:1:K
        T{l}=v'*T{l}*v;
        U{l}=U{l}*v;
    end
    V=V*v;
    
    % Perform 1-GMD on selected blocks of current stage
    [m,nn]=size(T{1});
    [u_t, t_t, v] =  one_gmd(T{i+1}(1:n,1:n));
    V = V*kron(eye(nn/n), v);
    for l = 1:1:K
        [u{l}, t{l}]=my_qr(T{l}(1:n,1:n)*v);
        U{l}=(U{l}*kron(eye(nn/n), u{l}));
        T{l} = kron(eye(nn/n), t{l});
    end
end
for i = 1:K
    T{i} = U{i}' * kron(eye(N), A{i}) * V;
end
end

function [Q, R, P] = one_gmd(A)
% Performs 1-GMD
% [Q, R, P] = one_gmd (A)

[U,S,V] = svd(A);
[Q,R,P] = gmd(U,S,V);
end

function [Q, R] = my_qr(A)
% [Q, R] = my_qr(A)

[Q,R] = qr(A);

n=length(R);

for i = 1:1:length(diag(R))
    T=eye(n);
    T(i,i) = conj (R(i,i))/abs (R(i,i));
    Q=Q*T;
    R=T*R;
end
end


% Matlab implementation of the "Geometric Mean Decomposition"
% version of Hager, December 3, 2003
% slightly modified by Yi, April 19, 2004
% Copyright 2003, University of Florida, Gainesville, Florida
%
%A = U*S*V' is the singular value decomposition of A
%           U, V unitary, S diagonal matrix with nonnegative
%           diagonal entries in decreasing order
%  = Q*R*P' is the geometric mean decomposition of A
%           P, Q unitary, R real upper triangular with r_ii =
%           geometric mean of the positive singular values of A,
%           1 <= i <= p, p = number of positive singular values
% All singular values smaller than tol treated as zero

function [Q, R, P] = gmd(U, S , V, tol)
if ( nargin < 4 )
    tol = eps ;
end
[m n] = size (S) ;
R = zeros (m, n) ;
P = V ;
Q = U ;
d = diag (S) ;
l = min (m, n) ;
for p = l : -1 : 1
    if ( d (p) >= tol )
        break ;
    end
end
if ( p < 1 )
    return ;
end
if ( p < 2 )
    R (1, 1) = d (1) ;
    return ;
end
z = zeros (p-1, 1) ;
large = 2 ;           % largest diagonal element
small = p ;           % smallest diagonal element
perm = [1 : p] ;      % perm (i) = location in d of i-th largest entry
invperm = [ 1 : p ] ; % maps diagonal entries to perm
sigma_bar = (prod (d (1:p)))^(1/p) ;
for k = 1 : p-1
    flag = 0 ;
    if ( d (k) >= sigma_bar )
        i = perm (small) ;
        small = small - 1 ;
        if ( d (i) >= sigma_bar )
            flag = 1 ;
        end
    else
        i = perm (large) ;
        large = large + 1 ;
        if ( d (i) <= sigma_bar )
            flag = 1 ;
        end
    end
    
    k1 = k + 1 ;
    if ( i ~= k1 )            % Apply permutation Pi of paper
        t = d (k1) ;          % Interchange d (i) and d (k1)
        d (k1) = d (i) ;
        d (i) = t ;
        j = invperm (k1) ;    % Update perm arrays
        perm (j) = i ;
        invperm (i) = j ;
        I = [ k1 i ] ;
        J = [ i k1 ] ;
        Q (:, I) = Q (:, J) ; % interchange columns i and k+1
        P (:, I) = P (:, J) ;
    end
    
    delta1 = d (k) ;
    delta2 = d (k1) ;
    t = delta1 + delta2 ;
    if ( flag )
        c = 1 ;
        s = 0 ;
    else
        f = (delta1 - sigma_bar)/(delta1 - delta2) ;
        s = sqrt (f*(delta1+sigma_bar)/t) ;
        c = sqrt(1-s^2) ;
    end
    d (k1) = delta1*delta2/sigma_bar ;          % = y in paper
    z (k) = s*c*(delta2 - delta1)*t/sigma_bar ; % = x in paper
    R (k, k) = sigma_bar ;
    if ( k > 1 )
        R (1:k-1, k) = z (1:k-1)*c ; % new column of R
        z (1:k-1) = -z (1:k-1)*s ;   % new column of Z
    end
    
    G1 = [ c -s
        s  c ] ;
    J = [ k k1 ] ;
    P (:, J) = P (:, J)*G1 ;        % apply G1 to P
    
    G2 = (1/sigma_bar)*[ c*delta1 -s*delta2
        s*delta2  c*delta1 ] ;
    Q (:, J) = Q (:, J)*G2 ;        % apply G2 to Q
end

R (p, p) = sigma_bar ;
R (1:p-1, p) = z ;
end