%aux_matrix creates a large matrix with copies of nextPrices to fit to the
%observed states in periods t

function res = aux_matrix(next_prices,prices,I_max, t)

I_max=I_max+1;                  %labeling in Matlab
fac1=numel(prices);             %number of different prices  
fac2=I_max*fac1;                %number of inventory-price combinations

res=zeros(fac1^(t-1)*I_max,fac1);
start=1;
a=1;
b=fac2;
while a<=length(res)
    ende=start+(fac1-1); %start:ende = rows that will be copied
    res(a:b,:)=repmat(next_prices(start:ende,:),I_max,1); %copies the rows times the number of possible inventory states
    start=ende+1;
    a=b+1;
    b=b+fac2;
end








