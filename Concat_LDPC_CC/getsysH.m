% use row operations to derive the systematic form of a parity check matrix H
% H: parity check matrix whose right hand square sub-matrix has full rank

function Hsys=getsysH(H)
[m,n]=size(H);
k=n-m;
for i=k+1:n
    ind=find(H(:,i),1,'last');
    % exchange (ind)th row and (i-k)th row
    if ind<i-k
        continue;
    end
    if ind~=i-k
        temp=H(ind,:);
        H(ind,:)=H(i-k,:);
        H(i-k,:)=temp;
    end
    I=find(H(:,i));
    % Guassian elimination
    for j=1:length(I)
        if I(j)~=i-k
            H(I(j),:)=mod(H(I(j),:)+H(i-k,:),2);
        end
    end
end
Hsys=H;
return