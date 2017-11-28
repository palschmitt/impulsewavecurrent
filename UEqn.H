    MRF.correctBoundaryVelocity(U);

	//vector damp = Zero;
	//damp.replace(vector::X,U.internalField());
    //damp.replace(vector::Y,U.internalField());
    //damp.replace(vector::Z,U.internalField());
	//U.component(3);
	//const dimensioned<vector>& damp;
	//damp = Zero;
	//const dimensionedVector& Damp;// = vector(0., 0., 0.);
	
	//Damp.internalField().replace(vector::X,sand.internalField());
    //Damp.internalField().replace(vector::Y,sand.internalField());
    //Damp.internalField().replace(vector::Z,sand.internalField());
        
        
    fvVectorMatrix UEqn
    (
        fvm::ddt(rho, U) + fvm::div(rhoPhi, U)
      + MRF.DDt(rho, U)
      + turbulence->divDevRhoReff(rho, U)
      //+ sand*rho*U//*damp
      //+ Beach*rho*U
      - Source
     ==
        fvOptions(rho, U)
    );

    UEqn.relax();

    fvOptions.constrain(UEqn);

    if (pimple.momentumPredictor())
    {
        solve
        (
            UEqn
         ==
            fvc::reconstruct
            (
                (
                    mixture.surfaceTensionForce()
                  - ghf*fvc::snGrad(rho)
                  - fvc::snGrad(p_rgh)
                ) * mesh.magSf()
            )
        );

        fvOptions.correct(U);
    }
