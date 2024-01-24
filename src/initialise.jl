
"""
    make_hedgehog!(skyrmion, theta)
    
Create a hedgehog in the field `skyrmion`, with orientation given by `theta`.

"""
function make_hedgehog!(sk, α, prof::Function)

    for i in 1:sk.lp[1], j in 1:sk.lp[2]

        x = sk.x[1][i];
        y = sk.x[2][j];

        r = sqrt( x^2 + y^2 )

        sk.field[i,j,1] = sin(prof(r))*(  cos(α)*x + sin(α)*y )/r
        sk.field[i,j,2] = sin(prof(r))*( -sin(α)*x + cos(α)*y )/r
        sk.field[i,j,3] = cos(prof(r))

    end

end

function make_hedgehog!(sk, α)

    prof(r) = pi*(1.0 - tanh(r))
    make_hedgehog!(sk, α, prof)

end





