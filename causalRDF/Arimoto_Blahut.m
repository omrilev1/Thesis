%
%This is MATLAB Code to find channel capacity Cc
%using the Arimoto-Blahut Algorithm
%
%code by Reeno Joseph Baby [reenojoseph@gmail.com]
%4th Year Undergraduate student College of Engineering Trivandrum
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Let C be source symbol matrix and Pc be its probability matrix (column matrix)
%Let Y be channel output matrix and Qy be its probability matrix (column matrix)
%Let Pyc be transition probability matrix
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%Algorithm
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Let |C|=M , |Y|=N , and let F=[f1,f2,f3,....fM].
%Let e be some small positive number (in the range of 10^-7)
%Let j and k be indices having ranges j=1,2,3...M
%                                     k=1,2,3...N
%initialize Pc with element values pj=1/M
%initialize Qy with Qy=Pyc*Pc
%
%REPEAT UNTIL stopping point is reached:
%
%  fj=exp{sum all k[pkj*ln(pkj/qk)]} for j=1,2,3...M
%  x=F*Pc
%  Il=log2(x)
%  Iu=lof2(max(F))
%  IF(Iu-Il)<e THEN
%      Cc=Il    Cc-->channel capacity
%      STOP
%  ELSE
%      pj=fj*pj/x  for j=1,2,3...M
%      Qy=Pyc*Pc
%  END IF
%  END REPEAT
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

Iu=1;Il=0;
e=input('enter error e ');
Pyc=input('enter transition probability matrix Py/c');
[N,M]=size(Pyc);
Pc=(1/M).*ones(M,1);
Qy=Pyc*Pc;
F=zeros(1,M);
while (Iu-Il)>e
    
    for j=1:M
        temp=0;
        for k=1:N
            temp=temp+((Pyc(k,j))*log((Pyc(k,j))/Qy(k)));
        end
        F(j)=exp(temp);
    end
    x=F*Pc;
    Il=log2(x);
    Iu=log2(max(F));
    
    if (Iu-Il)<e
        Cc= Il;
        Cc
        
        Pc
        
        Qy
        
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % to clculate mutual information
        temp2=0
        for j=1:M
            temp1=0;
            for k=1:N
                temp1=temp1+((Pyc(k,j))*log2((Pyc(k,j))/Qy(k)));
            end
            temp2=temp2+Pc(j)*temp1;
        end
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
        I=temp2 %mutual information
        
        break;
    else
        Pc=(1/x).*F'.*Pc;
        Qy=Pyc*Pc;
    end
end


