function [N] = niceDigits(N,n)
%converts integer N into a string with preceding zeros to a total of n digits

N=num2str(round(N));Nl=length(N);
for ii=1:n-Nl
   N=['0' N];
end

end