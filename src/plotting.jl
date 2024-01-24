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

