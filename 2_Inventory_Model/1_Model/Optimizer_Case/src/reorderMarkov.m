function res = reorderMarkov(nextPrices, b_t, I_max, t)
%REORDER creates a large matrix with copies of 'nextPrices' to fit to the
%observed states in periods t

I_max=I_max+1; % labelimg in Matlab
fac1=numel(b_t);                %Number of different prices  ->5
fac2=I_max*fac1;                %Number of inventory states for each price-> 55

res=zeros(fac1^t*I_max,fac1);
start=1;
a=1;
b=fac2*fac1;
while a<=length(res)
    ende=start+(fac1-1); %start:ende = rows that will be copied
    res(a:b,:)=repmat(nextPrices(start:ende,:),fac2,1); %Copies the rows times the number of possible inventory states
    start=ende+1;
    a=b+1;
    b=b+fac1*fac2;
end








