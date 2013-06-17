
%
%IMPR_EULER_DIFF_EQN
%
%Calculates solution to first order linear differential equaltion
%with initial value problem. It uses Improved Euler Method.
%
%Calling Method:
%[x,y] = impr_euler_diff_eqn(df,x0,y0,xfrom,xto,h)
%
%Input Arguments:
%df    : It is the differential equation. It should be a function handle.
%        And df should be of the following form,
%        df = dy/dx = f1(x,y).
%        Or, df = dx/dt = f1(t,x).. etc...
%
%x0    : For initial value problem f(x0)=y0. This is the x0 value, i.e initial
%        value of x.
%y0    : y0, initial value of y.
%        NOTE: (x0,y0) forms a initial (x,y) coordinate for which the
%        solution is known.
%xfrom : we want the values from this point
%xto   : to this point.
%h     : the increment of h. for Improved Euler method h should be
%        sufficiently low in order to have a good result. if a value of
%        h is not passed it defaults to 0.01.

%{
-Mohammad Tariqul Islam
ponir.bd @ hotmail.com
%}
function [x,y] = impr_euler_diff_eqn(df,x0,y0,xfrom,xto,h)
    %if h is not passed h is set
    if nargin < 6
        h=0.01;
    end
    
    %sorting xto and xfrom and remenbering
    m=0;
    if xfrom > xto
        temp = xfrom;
        xfrom = xto;
        xto = temp;
        m=1;
    end
    
    %calculating different ranges stored for three cases
    if x0 <= xfrom && x0 <= xto
        x1 = [];
        x2 = x0:h:xto;
        n=0;
    elseif x0 >= xfrom && x0 <= xto
        x1 = xfrom:h:x0;
        x2 = x0:h:xto;
        n=1;
    elseif x0 >= xfrom && x0 >= xto
        x1 = xfrom:h:x0;
        x2 = [];
        n=2;
    else
        x1 = [];
        x2 = [];
    end
    
    %dealing x1
    len1 = length(x1);
    y1 = zeros(1,len1);
    
    %n= 1 or 2 means x1 has been set
    if n==1 || n==2
        %initializing (x0,y0) pair
        y1(len1) = y0;
        
        hm = - h/2;
        %applying improved euler's method for each x in x1
        for i=len1-1:-1:1
            y1(i) = y1(i+1) - h*df(x1(i+1),y1(i+1));
            y1(i) = hm * (df(x1(i+1),y1(i+1)) + df(x1(i),y1(i))) + y1(i+1);
        end
        
        %dropping (x0,y0) pair
        y1 = y1(1:len1-1);
        x1 = x1(1:len1-1);
    end
    
    %dealing with x2
    len2 = length(x2);
    y2 = zeros(1,len2);
    
    %n=0 or 1 means x2 is set and have to dealt
    if n==1 || n==0
        %initializing x0,y0
        y2(1) = y0;
        
        h=h/2;
        %for each x in xi improved euler's method is applied
        for i=2:len2
            y2(i) = y2(i-1) + h*df(x2(i-1),y2(i-1));
            y2(i) = h*(df(x2(i-1),y2(i-1)) + df(x2(i),y2(i)))+y2(i-1);
        end
    end
    
    %generating x and y pairs
    y = [y1 y2];
    x = [x1 x2];
    
    %getting length of x
    len = length(x);
    
    %obtaining the index of xform
    for i=1:len
        if x(i) >= xfrom
            j1 = i;
            break;
        end
    end
    
    %obtaining the index of xto
    for i=len:-1:1
        if x(i) <= xto
            j2 = i;
            break;
        end
    end
    
    %cutting from xfrom-xto
    %if m=1, i.e-xfrom and xto was swapped, the x,y pairs are reversed
    %otherwise they are the same
    if m
        x=x(j2:-1:j1);
        y=y(j2:-1:j1);
    else
        x=x(j1:j2);
        y=y(j1:j2);
    end
end

