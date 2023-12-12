using Pkg
Pkg.activate()
using GeoParams

CharDim = GEO_units(length=1km, viscosity=1.0e19Pas, stress=100MPa)

L = GeoUnit(500*m)
R = GeoUnit(8.314*J/mol/K)
T = GeoUnit(288*K)
mass = GeoUnit(0.028*kg)
ρa = GeoUnit(1.225*kg/m^3)
ρc = GeoUnit(2700.0*kg/m^3)
ρs = GeoUnit(2800.0*kg/m^3)
βa = GeoUnit(1.0e-6*1/Pa)
βs = GeoUnit(1.0e-11*1/Pa)
ηa = GeoUnit(1.0e17*Pa*s)
ηs = GeoUnit(1.0e21*Pa*s)
μa = GeoUnit(1.0e13*Pa)
μs = GeoUnit(1.0e10*Pa)
g  = GeoUnit(9.81*m/s^2)
P  = GeoUnit(1.0e5*Pa)

Lnon = nondimensionalize(L, CharDim);
Rnon = nondimensionalize(R, CharDim);
Tnon = nondimensionalize(T, CharDim);
massnon = nondimensionalize(mass, CharDim);
ρanon = nondimensionalize(ρa, CharDim);
ρcnon = nondimensionalize(ρc, CharDim);
ρsnon = nondimensionalize(ρs, CharDim);
βanon = nondimensionalize(βa, CharDim);
βsnon = nondimensionalize(βs, CharDim);
ηanon = nondimensionalize(ηa, CharDim);
ηsnon = nondimensionalize(ηs, CharDim);
μanon = nondimensionalize(μa, CharDim);
μsnon = nondimensionalize(μs, CharDim);
gnon  = nondimensionalize(g, CharDim);
Pnon  = nondimensionalize(P, CharDim);

@show ρanon
@show ρcnon
@show ρsnon
@show Pnon
@show βanon
@show βsnon
@show ηanon
@show ηsnon
@show μanon
@show μsnon
@show gnon; 
print("-------------------------------------------\n")