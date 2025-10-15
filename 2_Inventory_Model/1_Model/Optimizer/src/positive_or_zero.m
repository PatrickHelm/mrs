%positive_or_zero ensures that a value is set to 0 if it is <=0,
%y=pPositiveOrZero
%return x, if x>0, and 0 if x<=0; if x is a vector positive_or_zero is executed
%to every value of the vector

function res = positive_or_zero(x)

res=zeros(1,length(x));

for i=1:length(x)
    if x(i) <= 0
        res(i)=0;
    else
        res(i)=x(i);
    end
end