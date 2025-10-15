%nextPrice calculates the probabilities for all possible prices in t+1

function res = next_price(PI, P_pr_reg)

Pi_size=size(PI);
res=zeros(Pi_size(1), length(P_pr_reg));


for i=1:Pi_size(1)
    for j=1:length(P_pr_reg)
        res(i,j)=PI(i,1)*P_pr_reg(j,1)+PI(i,2)*P_pr_reg(j,2);
    end
end

if Pi_size(1)==1 % case for t=1, for consistency. 
    res=res';
end
