/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 1991-2010 OpenCFD Ltd.
     \\/     M anipulation  |
-------------------------------------------------------------------------------
License
    This file is part of OpenFOAM.

    OpenFOAM is free software: you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    OpenFOAM is distributed in the hope that it will be useful, but WITHOUT
    ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
    FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
    for more details.

    You should have received a copy of the GNU General Public License
    along with OpenFOAM.  If not, see <http://www.gnu.org/licenses/>.

Class
    Foam::cnoidal2wave::sinecnoidal2wave
    virtual cnoidal2wave class, different theories should be implemented in derivations
\*---------------------------------------------------------------------------*/
#include "cnoidal2wave.H"

namespace Foam
 {
  ////* * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

scalar Foam::cnoidal2wave::evalwavenumber()
{
	//According to Hunt J.N. 1979 Direct Solution of cnoidal2wave dispersion equation . Proc. ASCE,Jour. Waterway,Port,Coastal and Ocean Eng.,105:457-459
	scalar alpha 		= Foam::pow( (Foam::constant::mathematical::twoPi/Period_ ) , 2 )*Depth_/mag(maggravity_);
	scalar Palpha 		= 1+0.6522*alpha+0.4622*Foam::pow(alpha,2)+0.0864*Foam::pow(alpha,4)+0.0675*Foam::pow(alpha,5);
	return((Foam::constant::mathematical::twoPi/Period_/( Foam::sqrt( mag(maggravity_)*Depth_) )*Foam::sqrt(alpha+1/Palpha)));
}


vector Foam::cnoidal2wave::velocity(const vector & p)
{
	vector U;
	if ( p.component(2) <= Xi(p))
	{
		U = direction_*//horizontal Velocity
		maggravity_*wavenumber_*waveHeight_/(2.*omega_)*Foam::cosh(wavenumber_*(p.component(2)))/Foam::cosh(wavenumber_*Depth_)*Foam::cos(wavenumber_*xstar(p)-omega_*Time_) +
		-gravityUnitV_*//vertical Velocity
		maggravity_*wavenumber_*waveHeight_/(2.*omega_)*Foam::sinh(wavenumber_*(p.component(2)))/Foam::cosh(wavenumber_*Depth_)*Foam::sin(wavenumber_*xstar(p)-omega_*Time_);
	}
	else
	{
		U = vector(0.,0.,0.);
		}
	return(U);
	}


scalar Foam::cnoidal2wave::Xi(const vector &p)
{
	return(Depth_ +
	waveHeight_/2.*Foam::cos(wavenumber_*xstar(p)-omega_*Time_ ));
	}
	
scalar Foam::cnoidal2wave::p(const vector &p)
{
	return(rho_*maggravity_*(Xi(p)-p.component(2)));
	}
	
void Foam::cnoidal2wave::setTime(const scalar &Time)
{
	Time_ = Time;
	}


//Validity checks and info  for all monochromatic wave types derived from cnoidal2wave!!!
void Foam::cnoidal2wave::checkUrsell()
{

	Info<<"####################  Wave data       ####################### "<<endl;
	Info<<"Waterlevel :	"<<Depth_<<endl;
	Info<<"Period 	  :	"<<Period_<<endl;
	Info<<"WaveHeight :	"<<waveHeight_<<endl;
	Info<<"WaveLength :	"<<lambda_<<endl;
	Info<<"Ursell 	  :	"<<Ursell()<<endl;
	Info<<"rho water  :	"<<rho_<<endl;
	Info<<"############################################################# "<<endl;
	
	
	if ( (waveHeight_/lambda_ > 0.01 || Depth_/lambda_ > 1 ) )
	{
		Info<<"WARNING: Linear theory not suitable!!"<<endl;	
		}	
	}
					
	
//***************************Constructor********************************
Foam::cnoidal2wave::cnoidal2wave(IOdictionary waveDict, vector gravity, scalar waterdensity):wave(waveDict, gravity)
{
	rho_ = waterdensity;	
	waveHeight_= readScalar(waveDict.lookup("WaveHeight")); 
	Period_= readScalar(waveDict.lookup("Period")); 
	omega_	=  Foam::constant::mathematical::twoPi/Period_;	
	wavenumber_ = evalwavenumber();
	lambda_ = Foam::constant::mathematical::twoPi/wavenumber_;
	celerity_  = (Foam::sqrt(maggravity_/wavenumber_*Foam::tanh(wavenumber_*Depth_) ));
	checkUrsell();
	}

				

Foam::cnoidal2wave::cnoidal2wave(scalar waterdensity):wave(waterdensity)
{
	rho_ = waterdensity;	
	}

}

// End namespace Foam
 
// ************************************************************************* //
