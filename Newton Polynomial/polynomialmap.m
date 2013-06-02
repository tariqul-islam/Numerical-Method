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