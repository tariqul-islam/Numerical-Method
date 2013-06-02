%{
LAGRANGEPOLY
The function calculates lagrange polynomial.
Calling Method: lagrangepoly(x,y)
Input: It takes input two 1xlen vectors
Output: Outputs polynomial of degree len-1
%}
function P = lagrangepoly(x, y)
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
   
   %creating a empty 1xlen container
   P = zeros(1, len);
   for i=1:1:len
       %poly evaluates all the convolution of [1 -x(j)] except at x(i)
       %prod evaluates all the product of (x(i) - x(j)) except at x(i)
       P = P + (poly(x((1:len)~=i)) ./ prod(x(i)-x((1:len)~=i))) .* y(i);
   end
end