module ThinPlateSplinesGeometryBasicsExt
    using GeometryBasics: Point
    import ThinPlateSplines: tps_solve, tps_deform, ThinPlateSpline

    function tps_solve(
        x::AbstractVector{<:Point},
        y::AbstractVector{<:Point},
        λ;
        compute_affine = true
    )
        x = stack(x; dims=1)
        y = stack(y; dims=1)
        return tps_solve(x, y, λ; compute_affine)
    end

    function tps_deform(x2::AbstractVector{<:Point}, tps::ThinPlateSpline)
        x2 = stack(x2; dims=1)
        deformed = tps_deform(x2, tps)
        return Point{size(deformed,2)}.(eachrow(deformed))
    end
end
