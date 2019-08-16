function [finalH] = Gallager_construction_LDPC(n,k,w_v,w_c)
%This code is to generate parity check matrix of LDPC code using Gallager's construction.
% Ref:
% [1] Y. Xiao, M -H Lee,Low complexity MIMO-LDPC CDMA systems over multipath channels,
% IEICE Transactions on Communications, v E89-B, n 5, May, 2006, p 1713-1717
% [2] Gallager classical LDPC paper

flag = 0;
cnt = 0; 
while (~flag && cnt < 1e6)
    if (n*w_v ~= (n-k)*w_c)
        warning('No valid regular LDPC code exists!');
        return
    end
    
    H_sub = zeros(k/w_c,n); % First sub-matrix; there are w_c such sub-matrices.
    H_pre = zeros(n-k,n);
    %% Generation of Basic Sub-matrix
    for i = 1:k/w_v
        for j = (i-1)*w_c+1:i*w_c
            H_sub(i,j) = 1;
        end
    end
    %% Permutation of Basic Sub-matrix
    H_pre(1:k/w_v,1:n) = H_sub;
    for t = 2:w_v
        x = randperm(n);
        H_sub_perm = H_sub(:,x);
        H_pre((t-1)*k/w_v+1:t*k/w_v,1:n) = H_sub_perm;
    end
    
    H_pre = H_pre(:,randperm(n));
    H_pre = H_pre(randperm(n-k),:);
    %% check if H has length 4 cycles
    girth = zeros(1,n-k);
    % We can test girth 4 by using the matrix O, see [1].
    O=H_pre*H_pre';
    for i=1:(n-k)
        O(i,i)=0;
    end
    for i=1:(n-k)
        girth(i)=max(O(i,:));
    end
    girth4=max(girth);
    
    if girth4<2
        flag = 1;
        finalH = H_pre;
    else
        cnt = cnt + 1;
        if (mod(cnt,500) == 0)
            display(strcat('cnt = ',num2str(cnt),' iterations'));
        end
    end
end
end

