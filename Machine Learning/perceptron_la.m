function [per,iter] = perceptron_la(X,Y, itermax)
%[per,iter] = perceptron_la(X,Y,itermax)
%
%    Implements perceprton learning algorithm
%
%INPUT:
%
%X = is the matrix with dimension Nxd
%    where N is the number of vectors
%    d is the dimension of vector
%    each vector is taken in a row
%
%Y = is the matrix with dimension 1xN or Nx1
%    containing the N classification parameters
%    .classification is done as described in getSign
%
%itermax = maximum number of iteration the algorithm will
%          do before stopping. defaults to 100000.
%
%OUTPUT:
%
%per (or w)
%    the w vector in the form [w0, w1, w2,...]
%
%iter
%    number of iterations required
%
%

%%TODO
% Xin, Xout
% Ein calculation
% Xout, Yout
% Eout Calculation

%Written by
%Tariqul Islam
%ponir(dot)bd(at)hotmail(dot)com
%
    if nargin==2
        itermax=100000;
    end

    szx=size(X);
    szy=size(Y);
    
    if (szx(1) ~= length(Y)) || (szy(1)~=1 && szy(2)~=1) || isempty(Y)
        error('Error in data');
    end
    
    N=szx(1);
    Xd = zeros(szx(1),szx(2)+1);
    for i=1:N
        Xd(i,:) = [1 X(i,:)];
    end
    X=Xd;
    Xd=[];
    
    
    per=zeros(1,length(X(1,:)));
    iter=0;
    while true
        mcarr=[];
        for i=1:N
            if getSign(X(i,:),per) ~= Y(i)
                mcarr=[mcarr i];
            end
        end
        
        len=length(mcarr);
        if len==0
            break;
        end
        ch = mcarr(ceil(getRand([1,len],1)));
        per=per+Y(ch).*X(ch,:);
        
        iter=iter+1;
        
        if iter>=itermax
            warning('PLA didn''t finish executing'); 
            break;
        end
    end
end

function y = getSign(x,f)
    y = sum(f.*x);
    
    if y>=0
        y=1;
    else
        y=-1;
    end
end

function y = getRand(x,n)
% Y = getRand(X,N)
% X = [X1,X2] boundary
% N = No of points needed
    if length(x) ~=2
        y=0;
    else
        y = x(1) + (x(2)-x(1))*rand(1,n);
    end
end