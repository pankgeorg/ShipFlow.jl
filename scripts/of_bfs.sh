set -e
C=/tmp/bfs; rm -rf $C; mkdir -p $C/system $C/constant $C/0; cd $C
# Clean expansion-ratio-2 backward-facing step, H=1, Re_H=5000.
# Inlet channel y∈[1,2] x∈[-4,0]; downstream y∈[0,2] x∈[0,20].

cat > system/blockMeshDict <<'EOF'
FoamFile{version 2.0; format ascii; class dictionary; object blockMeshDict;}
scale 1;
vertices
(
 (-4 1 0) (0 1 0) (0 2 0) (-4 2 0)
 (0 0 0) (20 0 0) (20 1 0) (20 2 0)
 (-4 1 0.1) (0 1 0.1) (0 2 0.1) (-4 2 0.1)
 (0 0 0.1) (20 0 0.1) (20 1 0.1) (20 2 0.1)
);
blocks
(
 hex (0 1 2 3 8 9 10 11) (80 40 1) simpleGrading (1 1 1)   // inlet channel
 hex (4 5 6 1 12 13 14 9) (300 40 1) simpleGrading (1 1 1) // downstream lower
 hex (1 6 7 2 9 14 15 10) (300 40 1) simpleGrading (1 1 1) // downstream upper
);
boundary
(
 inlet { type patch; faces ((0 8 11 3)); }
 outlet { type patch; faces ((5 6 14 13) (6 7 15 14)); }
 upperWall { type wall; faces ((3 11 10 2) (2 10 15 7)); }
 lowerWall { type wall; faces ((0 1 9 8) (4 12 13 5) (4 1 9 12)); }
 frontBack { type empty; faces ((0 3 2 1)(8 9 10 11)(4 5 6 1)(12 9 14 13)(1 6 7 2)(9 10 15 14)); }
);
EOF

cat > constant/transportProperties <<'EOF'
FoamFile{version 2.0; format ascii; class dictionary; object transportProperties;}
transportModel Newtonian;
nu 2e-4;
EOF
cat > constant/turbulenceProperties <<'EOF'
FoamFile{version 2.0; format ascii; class dictionary; object turbulenceProperties;}
simulationType RAS;
RAS { RASModel kOmegaSST; turbulence on; printCoeffs on; }
EOF

cat > 0/U <<'EOF'
FoamFile{version 2.0; format ascii; class volVectorField; object U;}
dimensions [0 1 -1 0 0 0 0]; internalField uniform (1 0 0);
boundaryField{
 inlet{type fixedValue; value uniform (1 0 0);}
 outlet{type inletOutlet; inletValue uniform (0 0 0); value uniform (1 0 0);}
 "(upperWall|lowerWall)"{type noSlip;}
 frontBack{type empty;}
}
EOF
cat > 0/p <<'EOF'
FoamFile{version 2.0; format ascii; class volScalarField; object p;}
dimensions [0 2 -2 0 0 0 0]; internalField uniform 0;
boundaryField{
 inlet{type zeroGradient;}
 outlet{type fixedValue; value uniform 0;}
 "(upperWall|lowerWall)"{type zeroGradient;}
 frontBack{type empty;}
}
EOF
cat > 0/k <<'EOF'
FoamFile{version 2.0; format ascii; class volScalarField; object k;}
dimensions [0 2 -2 0 0 0 0]; internalField uniform 3.75e-3;
boundaryField{
 inlet{type fixedValue; value uniform 3.75e-3;}
 outlet{type inletOutlet; inletValue uniform 3.75e-3; value uniform 3.75e-3;}
 "(upperWall|lowerWall)"{type kqRWallFunction; value uniform 3.75e-3;}
 frontBack{type empty;}
}
EOF
cat > 0/omega <<'EOF'
FoamFile{version 2.0; format ascii; class volScalarField; object omega;}
dimensions [0 0 -1 0 0 0 0]; internalField uniform 3.4;
boundaryField{
 inlet{type fixedValue; value uniform 3.4;}
 outlet{type inletOutlet; inletValue uniform 3.4; value uniform 3.4;}
 "(upperWall|lowerWall)"{type omegaWallFunction; value uniform 3.4;}
 frontBack{type empty;}
}
EOF
cat > 0/nut <<'EOF'
FoamFile{version 2.0; format ascii; class volScalarField; object nut;}
dimensions [0 2 -1 0 0 0 0]; internalField uniform 0;
boundaryField{
 inlet{type calculated; value uniform 0;}
 outlet{type calculated; value uniform 0;}
 "(upperWall|lowerWall)"{type nutkWallFunction; value uniform 0;}
 frontBack{type empty;}
}
EOF

cat > system/controlDict <<'EOF'
FoamFile{version 2.0; format ascii; class dictionary; object controlDict;}
application simpleFoam; startFrom startTime; startTime 0; stopAt endTime;
endTime 5000; deltaT 1; writeControl timeStep; writeInterval 5000;
purgeWrite 1; writeFormat ascii; writePrecision 6; runTimeModifiable true;
EOF
cat > system/fvSchemes <<'EOF'
FoamFile{version 2.0; format ascii; class dictionary; object fvSchemes;}
ddtSchemes{default steadyState;}
gradSchemes{default Gauss linear;}
divSchemes{default none; div(phi,U) bounded Gauss linearUpwind grad(U);
 div(phi,k) bounded Gauss upwind; div(phi,omega) bounded Gauss upwind;
 div((nuEff*dev2(T(grad(U))))) Gauss linear;}
laplacianSchemes{default Gauss linear corrected;}
interpolationSchemes{default linear;}
snGradSchemes{default corrected;}
wallDist{method meshWave;}
EOF
cat > system/fvSolution <<'EOF'
FoamFile{version 2.0; format ascii; class dictionary; object fvSolution;}
solvers{
 p{solver GAMG; tolerance 1e-7; relTol 0.01; smoother GaussSeidel;}
 "(U|k|omega)"{solver smoothSolver; smoother symGaussSeidel; tolerance 1e-7; relTol 0.1;}
}
SIMPLE{nNonOrthogonalCorrectors 0; consistent yes;
 residualControl{p 1e-4; U 1e-4; "(k|omega)" 1e-4;}}
relaxationFactors{equations{U 0.9; "(k|omega)" 0.7;}}
EOF

blockMesh > log.blockMesh 2>&1
grep -i "bounding box" log.blockMesh | head -1
simpleFoam > log.simpleFoam 2>&1
tail -1 log.simpleFoam
grep -m1 "solution converged" log.simpleFoam || echo "(ran to endTime)"
LAST=$(ls -d [0-9]* | sort -n | tail -1); echo "lastTime=$LAST"
# Sample U at first-cell height above the flat lower wall (y≈Δy/2=0.025) for x>0.
cat > system/sdict <<EOF
type sets; libs (sampling); interpolationScheme cellPoint; setFormat raw; fields (U);
sets ( nw { type uniform; axis x; start (0.05 0.026 0.05) ; end (19.5 0.026 0.05); nPoints 390; } );
EOF
postProcess -func sdict -time $LAST > /dev/null 2>&1
F=$(find postProcessing -name "nw_U.xy" | head -1)
# main reattachment: rightmost neg->pos crossing of Ux (skip tiny corner region x<0.2)
awk '{x=$1;ux=$2; if(ux<0)neg=1; if(neg==1&&ux>=0){xr=x;neg=0}} END{printf "x_r=%.3f  x_r/H=%.2f\n",xr,xr}' "$F"
echo "=== Ux(x) near wall (every 20th) ==="
awk 'NR%20==1{printf "x=%.2f Ux=%.4f\n",$1,$2}' "$F" | head -20
