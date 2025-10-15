%GETDISTRIBUTION allocates an distribution of the demand

function demand_dist = get_demand_dist(dist,demand,variables)
[ch, cp, d_max, I_max, k, m, pi, P_pr_reg, t_max, prices, normal_mu, normal_sigma,poisson_lambda,nbin_r,nbin_p] = allocate(variables);
real_demand=demand+1; 
x=0:demand;

if strcmp(dist, 'uniform')  
    x=1:real_demand;
    demand_dist=unidpdf(x,real_demand);
elseif strcmp(dist, 'normal') 
    pd=makedist('Normal',normal_mu,normal_sigma);
    t=truncate(pd,0,demand);
    demand_dist=pdf(t,x);
elseif strcmp(dist,'poisson')
    pd=makedist('Poisson',poisson_lambda);
    t=truncate(pd,0,demand);
    demand_dist=pdf(t,x);
elseif strcmp(dist,'nbin')
    pd=makedist('NegativeBinomial',nbin_r,nbin_p);
    t=truncate(pd,0,demand);
    demand_dist=pdf(t,x);
elseif strcmp(dist,'manual')
    x=1:real_demand;
    demand_dist=[0 1];
end

 plot(0:demand,demand_dist);