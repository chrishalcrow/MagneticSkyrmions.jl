
function gradient_flow!(ϕ; steps = 1, dt=(ϕ.ls[1]*ϕ.ls[2])/50.0, tolerance = 0.0, checks = max(100,steps), print_stuff = true, dEdp = zeros(ϕ.lp[1], ϕ.lp[2], 3) )

    println("Time step is ", dt )

    if tolerance == 0 && checks > steps
        checks = steps
    end
    
    if print_stuff == true
        println("initial: energy: ", energy(ϕ) )

    end

    counter = 0
    prev_error = 1.0e9
    
    while counter < steps
        
        gradient_flow_for_n_steps!(ϕ,dEdp,checks,dt)
        
        err = max_abs_err(dEdp)
        if err > 3*prev_error
            error("Suspected numerical blowup. Please use a smaller dt. Currently, dt = ", dt)
        end
        prev_error = err

        counter += checks
        
        if print_stuff == true
            println("after ", counter, " steps, error = ", round(err, sigdigits=4))
            #println( round(err, sigdigits=8), "," )
        end

        if tolerance != 0.0    # => we are in tol mode    
            if err < tolerance
                counter = steps + 1    # => end the while loop
            else
                steps += checks    # => continue the while loop
            end
        end

    end

    if print_stuff == true
        println("final energy: ", energy(ϕ) )
    end

    return

end

function gradient_flow_for_n_steps!(ϕ,dEdp,n,dt)
    for _ in 1:n
        gradient_flow_1_step!(ϕ,dEdp,dt)
    end
end

function gradient_flow_1_step!(sk, dEdp, dt)

    getdEdp!(sk, dEdp)
    sk.field .+= dt.*dEdp;
    normer!(sk)
   
end 

function getdEdp!(sk, dEdp)

    (d, a, J, K, K2, D) = (sk.d, sk.a, sk.J, sk.K, sk.K2, sk.D)

    for j in sk.sum_grid[2], i in sk.sum_grid[1]
        
        p = getP(sk,i,j)
        dp = getDP(sk ,i, j)
        ddp = getDDP(sk, i, j)

        dEdp[i,j,1] = -2.0*D*d/a^2*(dp[2,3])
        dEdp[i,j,2] = 2.0*D*d/a^2*(dp[1,3])
        dEdp[i,j,3] = 2.0*D*d/a^2*(dp[2,1] - dp[1,2])
        dEdp[i,j,3] +=  K2*d*p[3]/a^3
        dEdp[i,j,3] +=  K*d/a^3

        for b in 1:3
            dEdp[i,j,b] -= d*sk.Bfield[i,j,b]/a^3
            for I in 1:2
                dEdp[i,j,b] += J*d*ddp[I,b]/a 
            end
        end

        DE_dot_field = dEdp[i,j,1]*p[1] + dEdp[i,j,2]*p[2] + dEdp[i,j,3]*p[3]

        for a in 1:3
            dEdp[i,j,a] -= p[a]*DE_dot_field
        end

    end

end

function max_abs_err(A)

    return maximum(abs.(A)) 

end



function time_flow!(ϕ, t0; RHS = zeros(ϕ.lp[1], ϕ.lp[2], 3), B,  steps = 1, dt=(ϕ.ls[1]*ϕ.ls[2])/50.0, checks = min(100,steps), print_stuff = true, dEdp = zeros(ϕ.lp[1], ϕ.lp[2], 3) )

    println(dt)
    
    if print_stuff == true
        println("initial static energy: ", energy(ϕ) )

    end

    counter = 0
    
    while counter < steps
        
        time_flow_for_n_steps!(ϕ,RHS,dEdp,B,checks,t0, dt)
        counter += checks
        
        if print_stuff == true
            println("after ", counter, " static energy is  ", energy(ϕ) )
        end

    end

    if print_stuff == true
        println("final energy: ", energy(ϕ) )
    end

    return

end

function time_flow_for_n_steps!(ϕ,RHS,dEdp,B,n,t, dt)
    for _ in 1:n
        time_flow_1_step!(ϕ,RHS,dEdp,B, t, dt)
    end
end

function time_flow_1_step!(ϕ, RHS, dEdp, B, t, dt)

    γ = 1.0;

    for a in 1:3
        ϕ.Bfield[:,:,a] .= B[a](t[1])
    end

    getmdotRHS!(ϕ, RHS, dEdp)

    ϕ.field .+= dt.*RHS

    normer!(ϕ)

    t[1] += dt;
   
end 


function getmdotRHS!(sk, RHS, dEdp)

    (d, a, J, K, K2,  D, α) = (sk.d, sk.a, sk.J, sk.K, sk.K2, sk.D, sk.α)

    ϵ = getϵ()

    for j in sk.sum_grid[2], i in sk.sum_grid[1]
        
        p = getP(sk,i,j)
        dp = getDP(sk ,i, j)
        ddp = getDDP(sk, i, j)

        dEdp[i,j,1] = -2.0*D*d/a^2*(dp[2,3])
        dEdp[i,j,2] = 2.0*D*d/a^2*(dp[1,3])
        dEdp[i,j,3] = 2.0*D*d/a^2*(dp[2,1] - dp[1,2])
        dEdp[i,j,3] +=  K2*d*p[3]/a^3
        dEdp[i,j,3] +=  K*d/a^3

        for b in 1:3
            dEdp[i,j,b] -= d*sk.Bfield[i,j,b]/a^3
            for I in 1:2
                dEdp[i,j,b] += J*d*ddp[I,b]/a 
            end
        end

        Mi = getMinv(p,α)


        for b in 1:3
            RHS[i,j,b] = 0.0
            for c in 1:3, d in 1:3, e in 1:3
                RHS[i,j,b] += Mi[b,c]*ϵ[c,d,e]*p[d]*dEdp[i,j,e]
            end
        end

        #DE_dot_field = RHS[i,j,1]*p[1] + RHS[i,j,2]*p[2] + RHS[i,j,3]*p[3]

        #for b in 1:3
        #    RHS[i,j,b] -= p[b]*DE_dot_field
        #end

    end

end


function getMinv(n,α)

    Mi = zeros(3,3);

    for i in 1:3
        Mi[i,i] = 1.0
        for j in 1:3
            Mi[i,j] += α^2*n[i]*n[j]
        end
    end

    Mi[1,2] += α*n[3]
    Mi[2,1] -= α*n[3]

    Mi[1,3] -= α*n[2]
    Mi[3,1] += α*n[2]

    Mi[2,3] += α*n[1]
    Mi[3,2] -= α*n[1]

    return Mi

end

function getM(n,α)

    Mi = zeros(3,3);

    for i in 1:3
        Mi[i,i] = 1.0
    end

    Mi[1,2] -= α*n[3]
    Mi[2,1] += α*n[3]

    Mi[1,3] += α*n[2]
    Mi[3,1] -= α*n[2]

    Mi[2,3] -= α*n[1]
    Mi[3,2] += α*n[1]

    return Mi

end

function getϵ()

    epsilon = zeros(Int64,3,3,3)

    epsilon[1,2,3] = 1
    epsilon[2,1,3] = -1

    epsilon[1,3,2] = -1
    epsilon[3,1,2] = 1

    epsilon[2,3,1] = 1
    epsilon[3,2,1] = -1

    return epsilon

end


