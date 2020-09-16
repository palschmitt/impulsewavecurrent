#include "createRDeltaT.H"

Info<< "Reading field p_rgh\n" << endl;
volScalarField p_rgh
(
    IOobject
    (
        "p_rgh",
        runTime.timeName(),
        mesh,
        IOobject::MUST_READ,
        IOobject::AUTO_WRITE
    ),
    mesh
);

Info<< "Reading field U\n" << endl;
volVectorField U
(
    IOobject
    (
        "U",
        runTime.timeName(),
        mesh,
        IOobject::MUST_READ,
        IOobject::AUTO_WRITE
    ),
    mesh
);

#include "createPhi.H"
     Info<< "Reading field sand\n" << endl;
     volScalarField sand
     (
         IOobject
         (
             "sand",
             runTime.timeName(),
             mesh,
             IOobject::MUST_READ,
             IOobject::AUTO_WRITE
         ),
         mesh
     );
     
     Info<< "Reading field Uana\n" << endl;
     volVectorField Uana
     (
         IOobject
         (
             "Uana",
             runTime.timeName(),
             mesh,
             IOobject::MUST_READ,
             IOobject::AUTO_WRITE
         ),
         mesh
     );

Info<< "Reading transportProperties\n" << endl;
immiscibleIncompressibleTwoPhaseMixture mixture(U, phi);

volScalarField& alpha1(mixture.alpha1());
volScalarField& alpha2(mixture.alpha2());

const dimensionedScalar& rho1 = mixture.rho1();
const dimensionedScalar& rho2 = mixture.rho2();


// Need to store rho for ddt(rho, U)
volScalarField rho
(
    IOobject
    (
        "rho",
        runTime.timeName(),
        mesh,
        IOobject::READ_IF_PRESENT
    ),
    alpha1*rho1 + alpha2*rho2
);
rho.oldTime();

scalar rhowater = max(rho1.value(), rho2.value());
// Mass flux
surfaceScalarField rhoPhi
(
    IOobject
    (
        "rhoPhi",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::NO_WRITE
    ),
    fvc::interpolate(rho)*phi
);
     Info<< "Reading field relax\n" << endl;
     volScalarField relax
     (
         IOobject
         (
             "relax",
             runTime.timeName(),
             mesh,
             IOobject::MUST_READ,
             IOobject::AUTO_WRITE
         ),
         mesh
     );

 	volVectorField Source
     (
         IOobject
         (
             "Source",
             runTime.timeName(),
             mesh,
             IOobject::NO_READ,
             IOobject::AUTO_WRITE
             
         ),
 		Uana*rho*relax
     );
     
     volVectorField zVector
		(
	IOobject
	(
	"zVector",
	runTime.timeName(),
	mesh,
	IOobject::NO_READ,
	IOobject::AUTO_WRITE
	),
	mesh,
	Foam::vector(0,0,1)
	);
     
     volVectorField Beach
     (
         IOobject
         (
             "Beach",
             runTime.timeName(),
             mesh,
             IOobject::NO_READ,
             IOobject::AUTO_WRITE
             
         ),
 		sand*zVector*rho
     );
// Construct incompressible turbulence model
autoPtr<incompressible::turbulenceModel> turbulence
(
    incompressible::turbulenceModel::New(U, phi, mixture)
);


#include "readGravitationalAcceleration.H"
#include "readhRef.H"
#include "gh.H"


volScalarField p
(
    IOobject
    (
        "p",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::AUTO_WRITE
    ),
    p_rgh + rho*gh
);

label pRefCell = 0;
scalar pRefValue = 0.0;
setRefCell
(
    p,
    p_rgh,
    pimple.dict(),
    pRefCell,
    pRefValue
);

if (p_rgh.needReference())
{
    p += dimensionedScalar
    (
        "p",
        p.dimensions(),
        pRefValue - getRefCellValue(p, pRefCell)
    );
    p_rgh = p - rho*gh;
}

mesh.setFluxRequired(p_rgh.name());
mesh.setFluxRequired(alpha1.name());

// MULES compressed flux is registered in case scalarTransport FO needs it.
surfaceScalarField alphaPhiUn
(
    IOobject
    (
        "alphaPhiUn",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::NO_WRITE
    ),
    mesh,
    dimensionedScalar(phi.dimensions(), Zero)
);

#include "createMRF.H"
#include "createFvOptions.H"


 //Create pointer to wave and set to appropriate wavetype, better do something like turbulence later
 
 wave *pwave;//Pointer to wave base class, use autoPointer!!!
 
 IOdictionary waveDict
 (
 IOobject
 (
 "waveDict",
 U.time().constant(),//registry?
 U.db(),
 IOobject::MUST_READ,
 IOobject::NO_WRITE
 )
 );
 
 word Type   //ugly workaround, correct later
 (
 waveDict.lookup("Type")
 );
 	
 //create object of wavetype and point wave pointer to it
 if (Type == word("sinewave"))
 {
 	pwave = new sinewave(waveDict,g.value(),rhowater);//Assign pointer to wave of chosen and initialized type
 	}
 else if (Type == word("stokes2wave"))
 {
 	pwave = new stokes2wave(waveDict,g.value(),rhowater);//Assign pointer to wave of chosen and initialized type
 	}
 else if (Type == word("stokes5wave"))
 {
 	pwave = new stokes5wave(waveDict,g.value(),rhowater);//Assign pointer to wave of chosen and initialized type
 	}
 //else if (Type == word("stokes5waveseries"))
 //{
 	//pwave = new stokes5waveseries(waveDict,g.value(),rhowater);//Assign pointer to wave of chosen and initialized type
 	//}
 else if (Type == word("cnoidal2wave"))
 {
 	pwave = new cnoidal2wave(waveDict,g.value(),rhowater);//Assign pointer to wave of chosen and initialized type
 	}
 else if (Type == word("genericwave"))
 {
 	pwave = new genericwave(waveDict,g.value(),rhowater);//Assign pointer to wave of chosen and initialized type
 	}
 else if (Type == word("velocityseries"))
 {
 	pwave = new velocityseries(waveDict,g.value(),rhowater);//Assign pointer to wave of chosen and initialized type
 	}    
    
 else
 	Info<<"wave type not implemented yet!";
 
 relaxation  relaxo(pwave,mesh.C(),U, Uana,relax);

 
