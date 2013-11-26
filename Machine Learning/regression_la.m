function [w,Ein,Eout] = regression_la(Xin,Yin,Xout,Yout,lambda,dimension)
%[w,Ein,Eout] = regression_la(Xin,Yin,Xout,Yout,lambda,dimension)
%
%
%   Implements Regression analysis as a Machine Learning Algorithm
%
%TERMS:
%
%   x = [row vector] = [x1, x2, x3,...]
%   x = [row vector] = [x0, x1, x2,...]
%   in regular cases x0=1;
%   dx = [number of members in x vector including x0]
%
%INPUTS:
%   Xin
%       Input matrix. Includes row vectors of x. Each row contains
%       one vector x. Number of rows defines number of x vectors.
%       Number of colums defines number of dimensions
%   Yin
%       A comlumn vector or row vector, each ith member of which
%       corresponds to the group of ith vector in Xin.
%   NOTE: Number of rows in Xin and number of members in Yin is equal
%
%   Xout,Yout
%       Similar to Xin and Xout for estimating out of sample performence
%       empty vectors should be passed if Eout is not needed
%   lambda
%       The regularization parameter.
%       it corresponds to w'w<=C, C is the regularization parameter.
%       default value: 0
%   dimension
%       should be passed the number actual parameter of the x vector,
%       excluding x0, i.e dx-1.
%       default value: dx-1 
%
%OUTPUT:
%   w
%       regression parameters. in the form [w0, w1, w2,...]
%   Ein
%       in sample error
%   Eout
%       out of sample error
%
%

%writtten by
%Tariqul Islam
%ponir(dot)bd(at)hotmail(dot)com
%

    %at least 2 arguments needed
    if nargin <2
        error('not enough input arguments');
    end
    
    %checks the data
    [Xin,Yin,d] = data_sanity_check(Xin,Yin,dimension,'in traindata: ');

    [rowy,~] = size(Yin); %number of rows in Yin
    [~,dx] = size(Xin); %number of columns in Xin, dimension of x

    %takes default value of missing arguments

    if nargin<6
        dimension=dx-1;
        if nargin<5
            lambda=0;
            if nargin<3
                Xout=[];
                Yout=[]; 
            end
        end
    end
    
    %calculates w.
    w = ((Xin'*Xin+lambda*eye(d))\Xin')*Yin;
    w = w';
    
    %calculates Ein
    if nargout>=2
        Ein=0;
        for i=1:rowy
            if sign(sum(w.*Xin(i,:))) ~= Yin(i)
                Ein=Ein+1;
            end
        end
        Ein=Ein/rowy;
    end
    
    %calculates Eout
    if nargout>=3
        if nargin>=4
            [Xout,Yout] = data_sanity_check(Xout,Yout,dimension, 'in testdata: ');
            len = length(Yout);
            Eout=0;
            for i=1:len
                if sign(sum(w.*Xout(i,:))) ~= Yout(i)
                    Eout=Eout+1;
                end
            end
            Eout = Eout/len;
        else
            Eout=[];
        end
    end
end


function [Xin,Yin,d] = data_sanity_check(Xin,Yin,d,erstr)
    [rowy,dy] = size(Yin);
    [rowx,dx] = size(Xin);
    
    if dy~=1 && rowy==1
        Yin=Yin';
        rowy=dy;
    elseif dy==1
        
    else
        error([erstr 'classification data is not supported']);
    end
    
    if (rowy==0 || rowx==0) || (rowy~=rowx)
        error([erstr 'Data missing from training data :: rowx:' int2str(rowx) ' rowy:' int2str(rowy)]);
    end
    
    if rowy>=1
        if dx == d
            Xin = [ones(rowy,1) Xin];
            d=d+1;
        elseif dx == d+1
            d=dx;
        else
            error([erstr 'the dimension of the vectors doesn''t match']);
        end
    end
end