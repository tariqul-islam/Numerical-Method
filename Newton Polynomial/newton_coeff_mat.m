%{
NEWTON_COEFF_MAT
Uses matrix inversion method/linear equation solving method
to calculate coefficients a0,a1,...,an-1 of newton polynomial
Calling Method: newton_coeff_mat(x,y)
Input: two row vectors x and y that contains the (x,y) pairs.
Ouptut: a row vector containing the values of a0,a1,a2,...,aN-1

-Ponir
ponir.bd @ hotmail.com
%}
function a = newton_coeff_mat(x,y)
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
   
   %creates a matrix container of (len x len) order having all zero
   A=zeros(len);
   A(:,1)=ones(len,1); %setting first colum to 1
   
   %the matrix obtained is a lower triangular matrix
   %(x_i - x_i-n)*previous member of the matrix
   for i=2:1:len
      for j=2:1:i
         A(i,j)= A(i,j-1)*(x(i)-x(j-1));
      end
   end
   
   %computing (inv(A)*(y)T)T
   a = (A\y')';
end