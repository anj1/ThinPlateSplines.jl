module ThinPlateSplines

export tps_solve, tps_energy, tps_deform

type ThinPlateSpline
	λ  # Stiffness.
	x1 # control points
	Y  # Homogeneous control point coordinates
	Φ  # TPS kernel
	d  # Affine component
	c  # Non-affine component
end

# Thin-plate splines.
# based on http://en.wikipedia.org/wiki/Thin_plate_spline

tps_basis(r) = abs(r)<eps(r) ? 0.0 : r*r*log(r)

# x: matrix of size KxD
tps_kernel(x) = [tps_basis(norm(x[i,:] - x[j,:])) for i=1:size(x,1),j=1:size(x,1)]

# find solution of tps transformation
# compute_affine: compute affine component
# (required for some operations but takes additional time.)
function tps_solve(x,y,λ,compute_affine=true)
	# D: number of dimensions
	# K: number of data points
	K,D = size(x)

	# homogeneous coordinates
	X=cat(2,ones(K,1),x)
	Y=cat(2,ones(K,1),y)

	# compute TPS kernel
	Φ = tps_kernel(x)

	# full QR decomposition
	Q,r = qr(X,thin=false)
	q1 = Q[:,1:(D+1)]
	q2 = Q[:,(D+2):end]

	# warping coefficients
	c = q2*inv(q2'*Φ*q2 + λ*eye(K-D-1))*q2'*Y

	# affine component
	d = compute_affine ?  r\(q1'*(Y - Φ*c)) : []
	return ThinPlateSpline(λ,x,Y,Φ,d,c)
end

# Thin-plate spline bending energy at minimum
tps_energy(tps::ThinPlateSpline) = tps.λ*trace(tps.c*tps.Y')

# Deform points according to thin-plate spline,
# if affine component d and coefficients c are known.
# x1 are coordinates of control points.
# x2 are coordinates of points to be deformed.
function tps_deform(x2,tps::ThinPlateSpline)
	x1,d,c=tps.x1,tps.d,tps.c
	d==[] && throw(ArgumentError("Affine component not available; run tps_solve with compute_affine=true."))

	# deform
	y2 = zeros(size(x2,1),size(x2,2)+1)
	for i = 1 : size(x2,1)
		z = x2[i,:]
		defc = zeros(1,size(x2,2)+1)
		for j = 1 : size(x1,1)
			defc += tps_basis(norm(z - x1[j,:]))*c[j,:]
		end
		y2[i,:] = cat(2, 1.0, z)*d + defc
	end

	y2[:,2:end]
end

# Deform points according to thin-plate spline,
# solving to find c and d.
# x1 are coordinates of control points.
# y are coordinates of deformed control points;
# x2 are coordinates of points to be deformed.
function tps_deform(x1,x2,y,λ)
	tps = tps_solve(x1,y,λ,compute_affine=true)

	tps_deform_with_param(x1,x2,tps)
end

end