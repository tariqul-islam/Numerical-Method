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