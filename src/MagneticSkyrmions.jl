module MagneticSkyrmions

using StaticArrays, GLMakie

export Skyrmion, overview, normer!

include("initialise.jl")
export make_hedgehog!

include("properties.jl")
export energy


include("derivatives.jl")

include("diff.jl")
export gradient_flow!, time_flow!

include("plotting.jl")
export plot_skyrmion, plot_energy, plot_arrows

"""
    Skyrmion(lp::Int64, ls::Float64)
	Skyrmion([lpx,lpy], [lsx,lsy])
    
Create a skyrme field with `lp` lattice points and `ls` lattice spacing. 

"""
mutable struct Skyrmion
	field::Array{Float64, 3}
    Bfield::Array{Float64, 3}
	lp::Vector{Int64}
	ls::Vector{Float64}
	x::Vector{StepRangeLen{Float64, Base.TwicePrecision{Float64}, Base.TwicePrecision{Float64}, Int64}}
	d::Float64
    hbar::Float64
    a::Float64
    J::Float64
    K::Float64 
    K2::Float64
    D::Float64
    α::Float64
	sum_grid::Vector{UnitRange{Int64}}
    index_grid_x::Vector{Int64}
	index_grid_y::Vector{Int64}
    BC::String
end

Skyrmion(lp::Vector{Int64}, ls::Vector{Float64}; d=1.0, hbar=1.0, a=1.0, J=1.0, K=1.0, K2=1.0, D=1.0, α=0.1,   BC="periodic" ) = Skyrmion( zeros(lp[1],lp[2], 3) ,zeros(lp[1],lp[2],3) ,lp, ls, [ -ls[a]*(lp[a] - 1)/2.0 : ls[a] : ls[a]*(lp[a] - 1)./2.0 for a in 1:2 ], d, hbar, a, J, K, K2, D, α, sum_grid(lp,BC), index_grid(lp[1]), index_grid(lp[2]), BC )

function sum_grid(lp::Vector{Int64}, BC::String)

	if BC == "periodic"
		return [ 1:lp[1], 1:lp[2] ]
	else
		return [ 3:lp[1]-2, 3:lp[2]-2 ]
	end

end


function index_grid(lp)

	index_grid_array = zeros(Int64, lp+4)

	for i in 1:lp+4
		index_grid_array[i] = mod1(i-2,lp)
	end

	return index_grid_array

end


"""
    overview(skyrmion)

Displays an overview of `skyrmion`'s properties.

"""
function overview(sk)

    println("This skyrmion is on a ", sk.lp[1],"x",sk.lp[2]," grid, with lattice spacing [", sk.ls[1],", ", sk.ls[2], "]. ")

    println()
    println("   d = ", round(sk.d, sigdigits=5) ) 
    println("hbar = ", round(sk.hbar, sigdigits=5) ) 
    println("   a = ", round(sk.a, sigdigits=5) ) 
    println("   J = ", round(sk.J, sigdigits=5) ) 
    println("   K = ", round(sk.K, sigdigits=5) ) 
    println("   D = ", round(sk.D, sigdigits=5) ) 


end



"""
    normer!(skyrmion)

Normalises `skyrmion`.

See also [`normer`]
"""
function normer!(sk::Skyrmion)

    for j in 1:sk.lp[2], i in 1:sk.lp[1]
        
        @inbounds normer = 1.0/sqrt( sk.field[i,j,1]^2 + sk.field[i,j,2]^2 + sk.field[i,j,3]^2  )
        for a in 1:3
            @inbounds sk.field[i,j,a] *= normer
        end

    end

end

end