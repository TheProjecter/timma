function [rows, cols,G_dec] = graycode(a)

rows = 2^floor(a/2);
cols = 2^(a-floor(a/2));
G_r_dec = grays_2(a);
G_r_bin = dec2bin(G_r_dec)-'0';
Index_G = 1:size(G_r_bin,1);
Index_G = reshape(Index_G,[cols,rows])';
Index_G(2:2:end,:) = fliplr(Index_G(2:2:end,:));
%  G = reshape(num2cell(G_r_bin(Index_G(:),:),2),rows,cols);
G_dec = reshape(G_r_dec(Index_G(:),:),rows,cols);

function  G = grays_2(n)
%  G = grays(n)  is a column of  2^n  distinct  n-bit
%  integers that step through the  Gray Codes  running
%  from  G(1) = 000...000  to  G(2^n) = 100...000 .
%  Each  G(j+1)  is obtained from  G(j)  by changing
%  just one bit of  G(j) .  To see which bit,  compute
%  whichbit = bitxor(G(j), G(j+1))  to get  2^m  for a
%  nonnegative integer  m < n .  Keep  n < 27  because
%  grays(n)  costs time and memory proportional to  2^n .
%  See also  graystep.m,  int2gray.m,  gray2int.m .
%  Display  G + 2^52  in  hex  to see how its last  n
%  bits change.                 W. Kahan,  8 July 2007

n = n(:) ;  T = length(n) ;
% if ( (T~=1)|(n~=round(n))|(n<1)|(n>26) ),  N = n ,
%   error(' grays(N)  needs a small positive integer  N .'),  
% end
G = zeros(2^n,1) ;  G(2) = 1 ;  T = 2 ;
for  k = 2:n
    T2 = T+T ;  ... = 2^k
    G(T+1:T2) = T + flipud(G(1:T)) ;
    T = T2 ;  
end
