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


%{
NEWTON_COEFF_BACKSUB
Uses back substitution method to calculate coefficients a0,a1,...,an-1
of newton polynomial
Calling Method: newton_coeff_backsub(x,y)
Input: two row vectors x and y that contains the (x,y) pairs.
Ouptut: a row vector containing the values of a0,a1,a2,...,aN-1

-Ponir
ponir.bd @ hotmail.com
%}

function a = newton_coeff_backsub(x,y)
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
   
   %takes the length of thevectors
   len = sz(2);
   
   %creating a vector of 1xlen to store the coefficient
   a = zeros(1,len);
   a(1) = y(1); %sets a0 = y0
   
   %computes other coefficients
   for i=2:1:len
       numerator = y(i) - a(1); %setting y(i) - a0
       denominator = 1; %initializing denominator
       
       for j = 1:1:i-2
           %coeff of a_j*(x_i - x_j)*a_(j+1)
           denominator = denominator*(x(i) - x(j));
           numerator = numerator - denominator*a(j+1);
       end
       denominator = denominator*(x(i) - x(i-1));
       a(i) = numerator/denominator; %storing a(i)
   end
end


%{
POLYNOMIALMAP
Maps polynomials to higher degree polynomials.
Calling Method: polynomialmap(source,degree)
Input: source - NxM matrix, suggesting N polynomials each of M-1 order.
       degree - An integer defining the degree of the output polynomials
                should be mapped. If degree is negetive and/or double it
                will be floored to nearest integer and converted to
                positive.
Output: Nx(M+degree-1) matrix; i.e. all the polynimials are zero padded in
        the front if degree>M-1 or truncated from the front (i.e Higher
        degree members are droppd) if degree<M-1 or unchanged if degree=M-1.

-Ponir
ponir.bd @ hotmail.com
%}

function y = polynomialmap(source, degree)
    %size of source
    sz = size(source);
    
    %checks if length of degree is 1
    if length(degree) ~= 1
       error('2nd argument must be a 1x1 vector.'); 
    end
    
    %if degree is not an integer, it it converts
    degree = abs(floor(degree));
    %length of padding, sz(2) == M
    len = degree-sz(2)+1;
    
    if(len>=0)
        %creates a zero matrix
        front= zeros(1, len);
        %allocates memory for y
        y = zeros(sz(1), len+sz(2));
        %does operation on each row
        for i=1:sz(1)
            y(i,:) = [front source(i,:)];
        end
    elseif len==0
        %unchanged if degree=M-1
        y=source;
    else
        %drops the higher degree values
        %because degree<M-1
        len = 1-len;
        y = source(:,len:sz(2));
    end
end

