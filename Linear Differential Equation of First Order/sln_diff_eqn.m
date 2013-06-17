%SLN_DIFF_EQN
%
%Calculates solution to first order linear differential equaltion
%with initial value problem. It uses Solution Method.
%Secant Method to solve the non-linear equation.
%
%Calling Method:
%[x,y] = sln_diff_eqn(df,x0,y0,xfrom,xto,h)
%
%Input Arguments:
%df    : It is the differential equation. It should be a function handle.
%        And df should be of the following form,
%        df = dy/dx = @(x,y) (f1(x,y))...
%        Or, df = dx/dt = @(t,x) (f1(t,x)).. etc...
%
%x0    : For initial value problem f(x0)=y0. This is the x0 value, i.e initial
%        value of x.
%y0    : y0, initial value of y.
%        NOTE: (x0,y0) forms a initial (x,y) coordinate for which the
%        solution is known.
%xfrom : we want the values from this point
%xto   : to this point for variable x or t
%h     : the increment of h. for this method method h=0.01 gives
%        remarkebly good reasult for most functions. if a value of h is not
%        passed h=0.01 is assumed.
%

%{
-Mohammad Tariqul Islam
ponir.bd @ hotmail.com
%}


function [x,y] = sln_diff_eqn(df,x0,y0,xfrom,xto,h)
    %if h is not passed h is set
    if nargin<6
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
    
    %calculating different range for three cases
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
        
        %value of h is negetive
        hm = -h;
        
        %solution method for each of x in x1
        for i=len1-1:-1:1
           %function/equation to be solved for y1
           f= @(z) (hm/2*(df(x1(i+1),y1(i+1)) + df(x1(i),z))+y1(i+1) - z);
           %folving y1
           y1(i) = secant_solution(f,1,2); 
        end
        
        %dropping x0,y0 pair
        y1 = y1(1:len1-1);
        x1 = x1(1:len1-1);
    end
    
    %dealing with x2
    len2 = length(x2);
    y2 = zeros(1,len2);
    
    %n-0 or 1 means x2 has been set
    if n==0 || n==1
        %setting x0,y0 pair
        y2(1) = y0;
        
        %solution method for each of x in x2
        for i=2:len2
            %equation to be solved for z=y1
            f= @(z) (h/2*(df(x2(i-1),y2(i-1)) + df(x2(i),z))+y2(i-1) - z);
            %solving y1
            y2(i) = secant_solution(f,1,2);
        end
    end
    
    %getting x,y
    y = [y1 y2];
    x = [x1 x2];
    
    %length of x
    len = length(x);
    
    %obtaining index of xfrom
    for i=1:len
        if x(i) >= xfrom
            j1 = i;
            break;
        end
    end
    
    %obtaining index for xto
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

%{
SECANT_SOLUTION

Calling Method:
y = secant_solution(f,xk,xk0);

Input Arguments:
f        : Function Handle of the equation to solve.
xk,xk0   : initial guesses of x.
%}

function y = secant_solution(f, xk, xk0)
    %number of iterations
    iter=0;
    %current value
    cval = f(xk);
    %previous value
    pval = f(xk0);
    
    %10^-8 tolerance level
    while abs(cval) > 10^-8 && iter < 50
        %getting new assumption for solution
        yk1 = xk - cval*(xk0-xk)/(pval-cval);
        
        %current value becomes previous value
        pval = cval;
        %new current value is valvulated
        cval = f(yk1);
        
        %same as above for x
        xk0 = xk;
        xk = yk1;
        
        %number of iteration
        iter = iter+1;
    end
    
    %setting y
    y = xk;
    
    %gives a warning
    if iter == 50
        warning('maximum iterations reached. solution may not be correct.');
    end
end

