%update_pi: Calculation of PI_{t}

function res = update_pi(k,m,pi,P_pr_reg,t,solution_method)

Pi_size=size(pi); %size of PI{t-1}
int_res=zeros(length(P_pr_reg)^(t-1),m); %size of the result, hence of PI{t}
index=1; 
    for l=1:Pi_size(1)
        for i=1:length(P_pr_reg)
            for j=1:m   % regimes
                if strcmp(solution_method, 'CEC')  
                    int_res(index,j)=pi(l,j);
                elseif strcmp(solution_method, 'R1')
                    int_res(index,1)=1;
                elseif strcmp(solution_method, 'R2')
                    int_res(index,2)=1;
                else  
                    pi_zaehler=0;
                     pi_nenner=0;
                        for n=1:2
                            pi_zaehler=pi_zaehler+pi(l,n)*k(n,j)*P_pr_reg(i,n); 
                            pi_nenner=pi_nenner+pi(l,n)*P_pr_reg(i,n);
                        end
                    int_res(index,j)=pi_zaehler/pi_nenner;    
                end 
            end
            index=index+1;
        end
    end

res=int_res;
    