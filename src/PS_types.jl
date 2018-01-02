abstract type AbstractShapeParams end

const R_multipole = 1.1

type ShapeParams <: AbstractShapeParams
	t::Array{Float64,1}
	ft::Array{Float64,2}
	dft::Array{Float64,2}
	R::Float64 #radius of multipole disk

	ShapeParams(t,ft,dft) = new(t,ft,dft,R_multipole*maximum(hypot.(ft[:,1],ft[:,2])))
end

type CircleParams <: AbstractShapeParams
	R::Float64 #radius of multipole disk = radius of circle
end

#### FMM #####

mutable struct FMMgroup
    #TODO: make immutable and change divideSpace accordingly.
    point_ids::Vector{Int64}
    center::Array{Float64,2}
	size::Int64

	FMMgroup(point_ids,center) = new(point_ids,center,length(point_ids))
end

import Base.push!
function push!(g::FMMgroup,p)
	push!(g.point_ids,p)
	g.size = length(g.point_ids)
end

mutable struct FMMbuffer
    rhs::Vector{Complex{Float64}}
    pre_agg::Array{Complex{Float64},2}
    trans::Vector{Complex{Float64}}

    FMMbuffer(Ns::Integer,P::Integer,Q::Integer,Ngroups::Integer) =
                                new(zeros(Complex{Float64},Ns*(2*P+1)),
                                    zeros(Complex{Float64},Q,Ngroups),
                                    zeros(Complex{Float64},Q))
end

mutable struct FMMmatrices
	Agg::Array{Array{Complex{Float64},2},1}
	Trans::Array{Array{Complex{Float64},1},1}
	Disagg::Array{Array{Complex{Float64},2},1}
	Znear::SparseMatrixCSC{Complex{Float64},Int64}
	groups::Vector{ParticleScattering.FMMgroup}
	P2::Int64
	Q::Int64
end

mutable struct FMMoptions
    FMM::Bool       #Is FMM used?
    nx::Integer     #number of groups in x direction (for division)
    dx::Real        #group height/width (alternative division)
    acc::Integer    #accuracy digits for translation truncation, and also for gmres if tol is not given
    tol::Real       #gmres tolerance
    method::String  #method used: can be "pre" or "pre2"
	# symmetric::Bool #are agg = disagg points, and thus Disagg[k] = Agg^*[k] ∀k?

    #empty contructor - for not using FMM
    FMMoptions() = new(false, 0, 0.0, 0, 0.0, "")

    #Full FMMoptions constructor with value checking
    function FMMoptions(FMM; nx = 0, dx = 0.0, acc = 0, tol = 0.0, method = "pre")
        FMM == false && return FMMoptions()
        (nx == 0 && dx == 0.0) &&
            error("FMMoptions: either nx or dx must be specified")
        (nx > 0 && dx > 0.0) &&
            error("FMMoptions: either nx or dx must be specified")
        nx < 0 &&
            error("FMMoptions: nx must be greater than 0")
        dx < 0.0 &&
            error("FMMoptions: dx must be greater than 0")
        !in(acc,1:16) &&
            error("FMMoptions: accuracy digits must be in [1,16]")
        tol < 0.0 &&
            error("FMMoptions: gmres tolerance must be greater than 0.0")
        !in(method,("pre","pre2")) &&
            error("FMMoptions: method must be \"pre\" or \"pre2\"")
        tol == 0.0 && (tol = 10^(-Float64(acc)))
        return new(true, nx, dx, acc, tol, method)
    end
end

### Optimization

mutable struct OptimBuffer
    β::Vector{Complex{Float64}}
    f::Vector{Complex{Float64}}
    ∂β::Array{Complex{Float64},2}
    rhs_grad::Vector{Complex{Float64}}

    OptimBuffer(Ns::Integer,P::Integer,Npoints::Integer) =
                                new(Array{Complex{Float64}}(Ns*(2*P+1)),
                                    Array{Complex{Float64}}(Npoints),
                                    Array{Complex{Float64}}(Ns*(2*P+1),Ns),
                                    Array{Complex{Float64}}(Ns*(2*P+1)))
	#if J ̸= Ns
	OptimBuffer(Ns::Integer,P::Integer,Npoints::Integer,J::Integer) =
                                new(Array{Complex{Float64}}(Ns*(2*P+1)),
                                    Array{Complex{Float64}}(Npoints),
                                    Array{Complex{Float64}}(Ns*(2*P+1),J),
                                    Array{Complex{Float64}}(Ns*(2*P+1)))
end

### Scattering

struct ScatteringProblem
	shapes::Vector{T} where T <: AbstractShapeParams
	ids::Vector{Int64}
	centers::Array{Float64,2}
	φs::Vector{Float64}

	ScatteringProblem(shapes,ids,centers,φs) =
		((length(ids) == size(centers,1)) &&
		(length(ids) == length(φs)) &&
		(maximum(ids) <= length(shapes))) ? new(shapes,ids,centers,φs) :
		error("ScatteringProblem: size mismatch")
end

import Base.size
size(q::ParticleScattering.ScatteringProblem) = length(q.ids)
