%{
NEWTON_COEFF_DIVDIFF
Uses divided difference method to calculate coefficients a0,a1,...,an-1
of newton polynomial
Calling Method: newton_coeff_divdiff(x,y)
Input: two row vectors x and y that contains the (x,y) pairs.
Ouptut: a row vector containing the values of a0,a1,a2,...,aN-1

-Ponir
ponir.bd @ hotmail.com
%}
function a = newton_coeff_divdiff(x,y)
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
    
   %creating a lenxlen order matrix
   A = zeros(len);
   
   %getting f(x) values in the first column
   A(:,1)=y';
   %creating a 1xlen vector for a 
   a = zeros(1,len);
   a(1) = A(1,1); %setting a0 = f(x0)
   
   %running the recursive algorithm
   for i=2:1:len
       for j=2:1:i
          %operation discussed in the algorithm statement
          A(i,j) = ( A(i,j-1) - A(i-1,j-1) ) / (x(i) - x(i-j+1)); 
       end
       a(i) = A(i,i); %setting a(i) = A(i,i), diagonal element of A
   end
end