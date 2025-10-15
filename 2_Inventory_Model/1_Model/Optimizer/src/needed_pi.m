%neededPI returns the reachable beliefs of being in regime 1 dependend on
%the observed price in period 1

function res = needed_pi(PI, observed_price, prices)

t=length(PI);
res=cell(1,t);

actual_price=prices==observed_price; 

for i=1:t
    if i==1         %Reachable belief in period t=1
        res{1}=PI{1}(1);
    elseif i==2     %Reachable belief in period t=2
        res{2}=PI{2}(actual_price,1);
    else
        counter=prices;
        for j=1:i-2 %Compute the reachable belief in period t>=3
            var=combvec(counter,prices); 
            counter=var;
            var=fliplr(var');
            ind=var(:,1)==observed_price;
            res{i}=PI{i}(ind,1);
        end
    end
end
