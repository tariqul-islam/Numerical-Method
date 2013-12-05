function [w,b,Ein,Eout]=rbf_bias(X,Y,Xout,Yout,K,XS,gamma)
%
%[w b Ein Eout] = rbf_bias(X,Y,Xout,Yout,K,SPACE,gamma)
%
%X,Y
%   input training data X, each x in row and corresponding y in each row of Y.
%
%Xout,Yout
%   similar as X,Y
%
%K
%   how many clusters are required
%SPACE
%   the input space. [lower_limit upper_limit]
%
%gamma
%   value of gamma
%

%written by
%tariqul islam ponir
%ponir(dot)bd(at)hotmail(dot)com
%    
    RBFK = @(x1,x2,g) (exp(-g*norm(x1-x2)^2));
    N=length(Y);
    Nout=length(Yout);

    u=kmeans_clustering(X,K,XS);
    PHI=zeros(N,K+1);
    PHI(:,1)=1;
    for i=1:N
        for k=1:K
            PHI(i,k+1)=RBFK(X(i,:),u(k,:),gamma);
        end
    end
    w = inv(PHI'*PHI)*PHI'*Y;
    b=w(1);
    w=w(2:length(w))';
    
    Ein=0;
    for i=1:N
        sum=b;
        for k=1:K
            sum=sum+w(k)*RBFK(X(i,:),u(k,:),gamma);
        end
        if(sign(sum)~=Y(i,:))
            Ein=Ein+1;
        end
    end
    
    Eout=0;
    for i=1:Nout
        sum=b;
        for k=1:K
            sum=sum+w(k)*RBFK(Xout(i,:),u(k,:),gamma);
        end
        if(sign(sum)~=Yout(i,:))
            Eout=Eout+1;
        end
    end
end