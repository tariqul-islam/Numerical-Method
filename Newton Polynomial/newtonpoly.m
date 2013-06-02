%{
NEWTONPOLY
calculates newton polynomial
Calling Method: newtonpoly(x,y)
Input: (x,y) - one dimentional vectors of equal length N
Output: Polynomial of degree N-1
%}
function P = newtonpoly(x,y)
    %checks if two inputs vectors are obtained
    if nargin ~= 2
        error('The function accepts only two arguments of equal length'); 
    end
   
    sz = size(x); %size of x
    sz2 = size(y); % size of y
   
    %checks if size of x and size of y matches and they are row vectors
    if (sz(1) ~= sz2(1)) || (sz(2) ~= sz2(2)) || (sz(1) ~= 1)
        error('Mismatch in length or unsupported arguments.');
    end
    
    len=sz(2); %length of each vector
    degree = len-1; %degree of the output polynimial
    
    %obtaining newton constants sinc it's the fastes method of the three
    cons = newton_coeff_backsub(x,y);
    
    %initializing polynomial
    P = polynomialmap(cons(1), degree);

    %another container for calculation purpose
    SubP=1;

    %calculating the polynomial
    for i=2:1:len
       SubP = conv(SubP, [1 -x(i-1)]);
       P = P + cons(i) * polynomialmap(SubP,degree);
    end
end