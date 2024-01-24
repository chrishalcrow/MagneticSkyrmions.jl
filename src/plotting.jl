function plot_skyrmion(sk)

    fig=Figure()

    ax1 = Axis3(fig[1,1])
    ax2 = Axis3(fig[1,2])
    ax3 = Axis3(fig[1,3])

    surface!(ax1, sk.x[1], sk.x[2], sk.field[:,:,1] )
    surface!(ax2, sk.x[1], sk.x[2], sk.field[:,:,2] )
    surface!(ax3, sk.x[1], sk.x[2], sk.field[:,:,3] )

    return fig

end

function plot_energy(sk)

    ED = energy(sk,density=true);

    fig=Figure()

    ax = Axis3(fig[1,1])
    surface!(ax,sk.x[1], sk.x[2],ED)

    return fig

end


function sk_arrows(P0::Skyrmion; SFreq=1, cut=0)

    Nx = P0.lp[1]
    Ny = P0.lp[2]

    xlr = P0.x[1][1+cut:SFreq:Nx-cut]
    ylr = P0.x[2][1+cut:SFreq:Ny-cut]

    Plr = P0.field[1+cut:SFreq:Nx-cut,1+cut:SFreq:Ny-cut,1:3]

    strength = zeros( size(Plr)[1], size(Plr)[2] )
    angles = zeros(HSL{Float64}, size(Plr)[1], size(Plr)[2] )

    for i in 1:size(Plr)[1], j in 1:size(Plr)[2]
        strength[i,j] = exp(-0.01*(Plr[i,j,1]^2 + Plr[i,j,2]^2))
        angles[i,j] = HSL(atand(Plr[i,j,2], Plr[i,j,1]), 1, 0.5 + 0.5*Plr[i,j,3] )
        #Plr[i,j,1:2] = normalize(Plr[i,j,1:2])
    end

    strength = vec(strength)
    angles = vec(angles)

    ps = [Point3f(xlr[i], ylr[j],0.0) for i in 1:size(xlr)[1] for j in 1:size(ylr)[1] ];
    ns = [Point3f(Plr[i,j,1],Plr[i,j,2],Plr[i,j,3]) for i in 1:size(xlr)[1] for j in 1:size(ylr)[1] ];

    fig = Figure()

    axs = [Axis3(fig[1, i];
        #xgridvisible = false,
        #ygridvisible = false    
        elevation = pi/2.0,
        azimuth = 0.0,
        perspectiveness=0.0
        ) for i=1:1]

    hidedecorations!(axs[1])

    arrows!(axs[1], ps, ns, arrowcolor = angles, linecolor = angles, lengthscale=0.01, arrowsize=0.1*SFreq)

    zlims!(-1,1)

    return fig

end
