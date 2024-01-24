function getP(sk,i,j)

    return SVector{3,Float64}(
        sk.field[i,j,1],
        sk.field[i,j,2],
        sk.field[i,j,3]
    )
    
end


function getDP(ϕ,i,j)
    if ϕ.BC == "periodic"
        return getDPp(ϕ,i,j)
    else
       return getDPnp(ϕ,i,j)
    end
end

function getDPnp(ϕ,i,j)

    return SMatrix{2,3,Float64, 6}(
        dxD(ϕ.field,1,i,j,ϕ.ls[1]),
        dyD(ϕ.field,1,i,j,ϕ.ls[2]),

        dxD(ϕ.field,2,i,j,ϕ.ls[1]),
        dyD(ϕ.field,2,i,j,ϕ.ls[2]),

        dxD(ϕ.field,3,i,j,ϕ.ls[1]),
        dyD(ϕ.field,3,i,j,ϕ.ls[2]),
    )
    
end



function getDPp(ϕ,i,j)

    return SMatrix{2,3,Float64, 6}(
        dxDp(ϕ.field,1,i,j,ϕ.ls[1],ϕ.index_grid_x),
        dyDp(ϕ.field,1,i,j,ϕ.ls[2],ϕ.index_grid_y),

        dxDp(ϕ.field,2,i,j,ϕ.ls[1],ϕ.index_grid_x),
        dyDp(ϕ.field,2,i,j,ϕ.ls[2],ϕ.index_grid_y),

        dxDp(ϕ.field,3,i,j,ϕ.ls[1],ϕ.index_grid_x),
        dyDp(ϕ.field,3,i,j,ϕ.ls[2],ϕ.index_grid_y),
    )
    
end

function getDDP(ϕ,i,j)
    if ϕ.BC == "periodic"
        return getDDPp(ϕ,i,j)
    else
       return getDDPnp(ϕ,i,j)
    end
end

function getDDPnp(ϕ,i,j)

    return SMatrix{3,3,Float64, 9}(
        d2xD(ϕ.field,1,i,j,ϕ.ls[1]),
        d2yD(ϕ.field,1,i,j,ϕ.ls[2]),
        dxdyD(ϕ.field,1,i,j,ϕ.ls[1], ϕ.ls[2]),

        d2xD(ϕ.field,2,i,j,ϕ.ls[1]),
        d2yD(ϕ.field,2,i,j,ϕ.ls[2]),
        dxdyD(ϕ.field,2,i,j,ϕ.ls[1], ϕ.ls[2]),

        d2xD(ϕ.field,3,i,j,ϕ.ls[1]),
        d2yD(ϕ.field,3,i,j,ϕ.ls[2]),
        dxdyD(ϕ.field,3,i,j,ϕ.ls[1], ϕ.ls[2]),
    )
    
end


function getDDPp(ϕ,i,j)

    return SMatrix{3,3,Float64, 9}(
        d2xDp(ϕ.field,1,i,j,ϕ.ls[1],ϕ.index_grid_x),
        d2yDp(ϕ.field,1,i,j,ϕ.ls[2],ϕ.index_grid_y),
        dxdyDp(ϕ.field,1,i,j,ϕ.ls[1], ϕ.ls[2],ϕ.index_grid_x,ϕ.index_grid_y),

        d2xDp(ϕ.field,2,i,j,ϕ.ls[1],ϕ.index_grid_x),
        d2yDp(ϕ.field,2,i,j,ϕ.ls[2],ϕ.index_grid_y),
        dxdyDp(ϕ.field,2,i,j,ϕ.ls[1], ϕ.ls[2],ϕ.index_grid_x,ϕ.index_grid_y),

        d2xDp(ϕ.field,3,i,j,ϕ.ls[1],ϕ.index_grid_x),
        d2yDp(ϕ.field,3,i,j,ϕ.ls[2],ϕ.index_grid_y),
        dxdyDp(ϕ.field,3,i,j,ϕ.ls[1], ϕ.ls[2],ϕ.index_grid_x,ϕ.index_grid_y),
    )
    
end


function dxD(pion_field, a, i, j, lsx)
    @fastmath @inbounds (-pion_field[i+2,j,a] + 8.0*pion_field[i+1,j,a] - 8.0*pion_field[i-1,j,a] + pion_field[i-2,j,a])/(12.0*lsx)
end
function dyD(pion_field, a, i, j, lsy)
    @fastmath @inbounds (-pion_field[i,j+2,a] + 8.0*pion_field[i,j+1,a] - 8.0*pion_field[i,j-1,a] + pion_field[i,j-2,a])/(12.0*lsy)
end

function d2xD(pion_field, a, i, j, lsx)
    @fastmath @inbounds (-pion_field[i+2,j,a] + 16.0*pion_field[i+1,j,a] - 30.0*pion_field[i,j,a] + 16.0*pion_field[i-1,j,a] - pion_field[i-2,j,a])/(12.0*lsx^2)
end
function d2yD(pion_field, a, i, j, lsx)
    @fastmath @inbounds (-pion_field[i,j+2,a] + 16.0*pion_field[i,j+1,a] - 30.0*pion_field[i,j,a] + 16.0*pion_field[i,j-1,a] - pion_field[i,j-2,a])/(12.0*lsx^2)
end

function dxdyD(pion_field, a, i, j, lsx, lsy)
    @fastmath @inbounds 0.5*( dxdydiffD(pion_field, a, i, j,  lsx, lsy) - d2xD(pion_field, a, i, j,  lsx) - d2yD(pion_field, a, i, j,  lsx) )
end
function dxdydiffD(pion_field, a, i, j, lsx, lsy)
    @fastmath @inbounds (-pion_field[i+2,j+2,a] + 16.0*pion_field[i+1,j+1,a] - 30.0*pion_field[i,j,a] + 16.0*pion_field[i-1,j-1,a] - pion_field[i-2,j-2,a])/(12.0*lsx*lsy)
end

# periodic stuff

function dxDp(pion_field, a, i, j, lsx, ig)
    @fastmath @inbounds (-pion_field[ig[i+4],j,a] + 8.0*pion_field[ig[i+3],j,a] - 8.0*pion_field[ig[i+1],j,a] + pion_field[ig[i],j,a])/(12.0*lsx)
end
function dyDp(pion_field, a, i, j, lsy, ig)
    @fastmath @inbounds (-pion_field[i,ig[j+4],a] + 8.0*pion_field[i,ig[j+3],a] - 8.0*pion_field[i,ig[j+1],a] + pion_field[i,ig[j],a])/(12.0*lsy)
end


function d2xDp(pion_field, a, i, j,  lsx, ig)
    @fastmath  @inbounds (-pion_field[ig[i+4],j,a] + 16.0*pion_field[ig[i+3],j,a] - 30.0*pion_field[ig[i+2],j,a] + 16.0*pion_field[ig[i+1],j,a] - pion_field[ig[i],j,a])/(12.0*lsx^2)
end
function d2yDp(pion_field, a, i, j,  lsx, ig)
    @fastmath @inbounds (-pion_field[i,ig[j+4],a] + 16.0*pion_field[i,ig[j+3],a] - 30.0*pion_field[i,ig[j+2],a] + 16.0*pion_field[i,ig[j+1],a] - pion_field[i,ig[j],a])/(12.0*lsx^2)
end
function dxdydiffDp(pion_field, a, i, j, lsx, lsy, igx, igy)
    @fastmath  @inbounds (-pion_field[igx[i+4],igy[j+4],a] + 16.0*pion_field[igx[i+3],igy[j+3],a] - 30.0*pion_field[igx[i+2],igy[j+2],a] + 16.0*pion_field[igx[i+1],igy[j+1],a] - pion_field[igx[i],igy[j],a])/(12.0*lsx*lsy)
end

function dxdyDp(pion_field, a, i, j, lsx, lsy, igx, igy)
    @fastmath @inbounds 0.5*( dxdydiffDp(pion_field, a, i, j,  lsx, lsy, igx, igy) - d2xDp(pion_field, a, i, j,  lsx, igx) - d2yDp(pion_field, a, i, j,  lsx, igy) )
end
