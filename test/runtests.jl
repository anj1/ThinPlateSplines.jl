using Test
using ThinPlateSplines

@testset "tps generation" begin
    x1 = [0.0 1.0 
    1.0 0.0
    1.0 1.0]
    x2 = [0.0 1.0
      1.1 0.0
      1.2 1.5]
    tps = tps_solve(x1, x2, 1.0)
    @test tps.Y == [ 1.0  0.0  1.0
    1.0  1.1  0.0
    1.0  1.2  1.5]
    @test tps.c == zeros((3,3))
    @test tps.d ≈ [  1.0          -0.1  -0.5
    0   1.2   0.5
    0   0.1   1.5]
    @test isapprox(tps.Φ, [ 0.0       0.693147  0.0
    0.693147  0.0       0.0
    0.0       0.0       0.0], atol=1e-5)
end

x1 = [0.0 1.0 
1.0 0.0
1.0 1.0]
x2 = [0.0 1.0
  1.1 0.0
  1.2 1.5]
tps = tps_solve(x1, x2, 1.0)

@testset "tps_deform" begin
    
    x = [1.0 0.0
     2.0 2.0]
    y = tps_deform(x,tps)
    @test y ≈ [ 1.1  0
    2.5   3.5]
    y = tps_deform(x1, x, x2, 1.0)
    @test y ≈ [ 1.1  0
    2.5   3.5]
end

@testset "tps_energy" begin
    @test tps_energy(tps) ≈ 0
end

@testset "Three dimensions" begin
  start_pts = [0 0 0; 0 0 1; 0 1 0; 1 0 0]
  end_pts = [-0.7 -0.7 0; 0 0 1; 0 1 0; 1 0 0]
  tps = tps_solve(start_pts, end_pts, 1.0)
  deformed = tps_deform(start_pts, tps)
  @test size(deformed,1) == size(start_pts, 1)
  @test size(deformed,2) == size(start_pts, 2)
  @test deformed ≈ end_pts
  @test tps_deform([0.5 0.5 0.5], tps) ≈ [0.85 0.85 0.5]
end

@testset "Four dimensions" begin
  start_pts = [0 0 0 0; 0 0 0 1; 0 0 1 0; 0 1 0 0; 1 0 0 0]
  end_pts = [-0.7 -0.7 0 0; 0 0 0 1; 0 0 1 0; 0 1 0 0; 1 0 0 0]
  tps = tps_solve(start_pts, end_pts, 1.0)
  deformed = tps_deform(start_pts, tps)
  @test size(deformed,1) == size(start_pts, 1)
  @test size(deformed,2) == size(start_pts, 2)
  @test deformed ≈ end_pts
  @test tps_deform([0.5 0.5 0.5 0.5], tps) ≈ [1.2 1.2 0.5 0.5]
end