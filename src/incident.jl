#plane wave
function u2α(k, u::PlaneWave, centers::Array{T,2}, P) where T <: Real
	#build incoming coefficients for plane wave incident field
	α = Array{Complex{Float64}}(size(centers,1)*(2P + 1))
	for ic = 1:size(centers,1)
		phase = exp(1.0im*k*(cos(u.θi)*centers[ic,1] + sin(u.θi)*centers[ic,2])) #phase shift added to move cylinder coords
		for p = -P:P
			α[(ic-1)*(2P+1) + p + P + 1] = phase*exp(1.0im*p*(π/2-u.θi))
		end
	end
	α
end

function uinc(k, points::Array{T,2}, u::PlaneWave) where T <: Real
	exp.(1.0im*k*(cos(u.θi)*points[:,1] + sin(u.θi)*points[:,2]))
end

function uinc(k, points::Array{T,1}, u::PlaneWave) where T <: Real
	exp(1.0im*k*(cos(u.θi)*points[1] + sin(u.θi)*points[2]))
end

#line source
function u2α(k, u::LineSource, centers::Array{T,2}, P) where T <: Real
	#build incoming coefficients for plane wave incident field
	α = Array{Complex{Float64}}(size(centers,1)*(2P + 1))
	for ic = 1:size(centers,1)
		R = hypot(centers[ic,1] - u.x, centers[ic,2] - u.y)
		θ = atan2(centers[ic,2] - u.y, centers[ic,1] - u.x)
		for p = -P:P
			α[(ic-1)*(2P+1) + p + P + 1] = exp(-1im*p*θ)*besselh(-p, 1, k*R)
		end
	end
	α
end

function uinc(k, points::Array{T,2}, u::LineSource) where T <: Real
	r = hypot.(u.x - points[:,1], u.y - points[:,2])
	if any(r .== 0)
		warn("uinc: encountered singularity in incident field, returned NaN")
		r[r.==0] = NaN
	end
	besselh.(0, 1, k*r)
end

function uinc(k, points::Array{T,1}, u::LineSource) where T <: Real
	r = hypot(u.x - points[1], u.y - points[2])
	if r == 0
		warn("uinc: encountered singularity in incident field, returned NaN")
		return NaN
	end
	besselh(0, 1, k*r)
end
