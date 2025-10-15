def single_period_cost(y, cp, ch, p, d_max, I, demand_dist, dist):
    def positive_or_zero(x):
        return max(x, 0)

    def demand_prob(dist, d):
        # Placeholder for demand_prob function, needs actual implementation
        pass

    costs = p * y

    if dist == 'uniform':
        for d in range(d_max + 1):
            costs += ch / (d_max + 1) * positive_or_zero(I + y - d) + cp / (d_max + 1) * positive_or_zero(d - I - y)
    else:
        for d in range(d_max + 1):
            dProb = demand_prob(demand_dist, d)
            costs += ch * dProb * positive_or_zero(I + y - d) + cp * dProb * positive_or_zero(d - I - y)

    return costs