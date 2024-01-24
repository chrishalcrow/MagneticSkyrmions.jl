
"""
    energy(skyrmion)
    
Returns the energy of `skyrmion`.

"""

function energy(sk; density=false)

    energy_density = zeros(sk.lp[1], sk.lp[2])

    get_energy_density!(energy_density,sk)

    if density == false
        return sum(energy_density)*sk.ls[1]*sk.ls[2]
    else
        return energy_density
    end

end

function get_energy_density!(density, sk )

    for j in sk.sum_grid[2], i in sk.sum_grid[1]

        p = getP(sk, i, j)
        dp = getDP(sk ,i, j)

        density[i,j] = engpt(p,dp,sk.d,sk.K,sk.K2,sk.a,sk.J,sk.D)
        #density[i,j] += sk.d*(p[1]*sk.Bfield[i,j,1] + p[2]*sk.Bfield[i,j,2] + p[3]*sk.Bfield[i,j,3])/sk.a^3

    end

end

function engpt(p,dp,d,K,K2,a,J,D)

        #dEdp[i,j,3] +=  K*d*p[3]/a^3
        #dEdp[i,j,3] +=  K*d/a^3

        

    return (d/a^3)*( 1.0*K2*(1.0-p[3]^2)/2.0 +  K*(1.0 - p[3]) + J*a^2/2.0*(dp[1,1]^2 + dp[2,1]^2 + dp[1,2]^2 + dp[2,2]^2 + dp[1,3]^2 + dp[2,3]^2)  ) + 1.0*D*d*(dp[2,3]*p[1] - dp[1,3]*p[2] + dp[1,2]*p[3] - dp[2,1]*p[3])/a^2

end