module ThinPlateSplines

using LinearAlgebra
using Tullio

export ThinPlateSpline
export tps_solve, tps_energy, tps_deform

"""
	ThinPlateSpline{Λ,X,M}

the type (structure) holding the deformation information. This is needed to apply the deformation to other points using `tps_deform`.
(see http://en.wikipedia.org/wiki/Thin_plate_spline for more information)
	
# Members
`λ::Λ`  # Stiffness.
`x1::X` # control points
`Y::M`  # Homogeneous control point coordinates
`Φ::M`  # TPS kernel
`d::M`  # Affine component
`c::M`  # Non-affine component
"""
struct ThinPlateSpline{Λ,X <: AbstractMatrix,M}
	λ::Λ  # Stiffness.
	x1::X # control points
	Y::M  # Homogeneous control point coordinates
	Φ::M  # TPS kernel
	d::M  # Affine component
	c::M  # Non-affine component
end

# Thin-plate splines.
# based on http://en.wikipedia.org/wiki/Thin_plate_spline

# tps_basis is applied to the result of a norm which is either positive or zero
# using ifelse is faster than the (? x:y) notation
tps_basis(r::T) where {T}  = ifelse(r < eps(r), zero(T), r*r*log(r))

my_norm(a) = sqrt(sum(a.^2))

# x: matrix of size KxD
tps_kernel(x) = [tps_basis(my_norm(x[i,:] - x[j,:])) for i in axes(x)[1], j=axes(x)[1]]

"""
	tps_solve(x,y,λ,compute_affine=true)

find solution of tps transformation 
(required for some operations but takes additional time.)

# Arguments
	`x`: control points
	`y`: deformed (warped) control points
	`λ`: stiffness coefficient
	`compute_affine`: computes affine component if `true`

returns a `ThinPlateSpline` structure which can be supplied as an argument to `tps_deform`.

# See almost
	`ThinPlateSpline`
"""
function tps_solve(x,y,λ; compute_affine=true)
	# D: number of dimensions
	# K: number of data points
	K,D = size(x)

	# homogeneous coordinates
	X=cat(dims=2, ones(K,1),x)
	Y=cat(dims=2, ones(K,1),y)

	# compute TPS kernel
	Φ = tps_kernel(x)

	# full QR decomposition
	Q,r = qr(X)
	q1 = Q[:,1:(D+1)]
	q2 = Q[:,(D+2):end]

	# warping coefficients
	c = q2*inv(UniformScaling(λ) + q2'*Φ*q2)*q2'*Y

	# affine component
	d = compute_affine ?  r\(q1'*(Y - Φ*c)) : []
	return ThinPlateSpline(λ,x,Y,Φ,d,c)
end

"""
	tps_energy(tps::ThinPlateSpline)

Thin-plate spline bending energy at the minimum.

# Arguments
	`tps`: `ThinPlateSpline` structure defining the deformation
"""
tps_energy(tps::ThinPlateSpline) = tps.λ*tr(tps.c*tps.Y')

"""
	tps_deform(x2::AbstractMatrix, tps::ThinPlateSpline) 

calculate deformed locations of an input vector of positions`x2` according to thin-plate spline.
Note that in combination with suitable interpolation packages can this be used to warp an image, but this function does not apply a warp by itself.
Yet it can be used to transform other mathematical structures rather than points.

# Arguments
	`x2`: coordinates of points to be deformed.
	`tps`: the thin-plate spline structure of type `ThinPlateSpline` defining the deformation.
"""
function tps_deform(x2::AbstractMatrix{T}, tps::ThinPlateSpline) where {T}
    x1,d,c = tps.x1,tps.d,tps.c
	d==[] && throw(ArgumentError("Affine component not available; run tps_solve with compute_affine=true."))
	D = size(x2, 2)
    all_homo_z = hcat(ones(T, size(x2,1)), x2)
    # calculate sum of squares. Note that the summation is done outside the abs2
	# it may be useful to join the terms below, but this seems a little harder than first thought
    @tullio sumsqr[k,j] := abs2(all_homo_z[k,m+1] .- x1[j,m])
    @tullio yt[i,j] := (tps_basis(sqrt(sumsqr[i,k])) * c[k,j+1]) (j in 1:D)
    @tullio yt[i,j] += d[l,j+1] * all_homo_z[i,l] 

	# it would save some memory and possibly make it a bit faster, if the order of dimensions is changed for x and yt due to cache utilization. 
	# but the advantage is relatively minor (10%?) since @tullio seems to take care of this as well as possible
    return yt
end

"""
	tps_deform(x1, x2, y, λ; compute_affine=true) 

computes a tps deformation based on `x1`, `y` and `λ` and applies it to `x2`.
First `tps_solve` is applied and then `tps_deform`.

# Arguments
	`x1`: control points for the tps generation
	`x2`: coordinates of points to be deformed by the tps.
	`y`: deformed control points for the tpd generation
	`λ`: stiffness coefficient
	`compute_affine`: optional argument defining whether the affine transform is computed.

returned are the deformed coordinates as a matrix
"""
function tps_deform(x1, x2, y, λ; compute_affine=true)
	tps = tps_solve(x1, y, λ; compute_affine=compute_affine)
	tps_deform(x2, tps)
end

end
