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
    Foam::stokes2wave::sinestokes2wave
    virtual stokes2wave class, different theories should be implemented in derivations
\*---------------------------------------------------------------------------*/
#include "stokes2wave.H"

namespace Foam
 {
  ////* * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

scalar Foam::stokes2wave::evalwavenumber()
{
	//According to Hunt J.N. 1979 Direct Solution of stokes2wave dispersion equation . Proc. ASCE,Jour. Waterway,Port,Coastal and ocena Eng.,105:457-459
	scalar alpha 		= Foam::pow( (Foam::constant::mathematical::twoPi/Period_ ) , 2 )*Depth_/mag(maggravity_);
	scalar Palpha 		= 1+0.6522*alpha+0.4622*Foam::pow(alpha,2)+0.0864*Foam::pow(alpha,4)+0.0675*Foam::pow(alpha,5);
	return((Foam::constant::mathematical::twoPi/Period_/( Foam::sqrt( maggravity_*Depth_) )*Foam::sqrt(alpha+1/Palpha)));
}


vector Foam::stokes2wave::velocity(const vector & p)
{
	vector U;
	scalar kh = wavenumber_*Depth_;
	if ( p.component(2) <= Xi(p))
	{
		U = direction_*//horizontal Velocity
		(waveHeight_/2.*maggravity_*wavenumber_/omega_*Foam::cosh(wavenumber_*(p.component(2)))/cosh(kh)*
		Foam::cos(wavenumber_*xstar(p)-omega_*Time_) + 
		3./16.*Foam::pow(waveHeight_,2)*omega_*wavenumber_*Foam::cosh(2.*kh*(p.component(2)))/Foam::pow(sinh(kh),4)*
		Foam::cos(2*(wavenumber_*xstar(p)-omega_*Time_)))+

		-gravityUnitV_*//vertical Velocity
		(waveHeight_/2.*maggravity_*wavenumber_/omega_*Foam::sinh(wavenumber_*(p.component(2)))/cosh(kh)*
		Foam::sin(wavenumber_*xstar(p)-omega_*Time_) +
		3./16.*Foam::pow(waveHeight_,2)*omega_*wavenumber_*Foam::sinh(2.*kh*(p.component(2)))/Foam::pow(sinh(kh),4)*
		Foam::sin(2.*(wavenumber_*xstar(p)-omega_*Time_)));
	}
	else
	{
		U = vector(0.,0.,0.);
		}
	return(U);
	}


scalar Foam::stokes2wave::Xi(const vector &p)
{
	return(
	Depth_ + 
	waveHeight_/2.*Foam::cos( wavenumber_*xstar(p)-omega_*Time_)+
	Foam::pow(waveHeight_,2)*wavenumber_/16.*Foam::cosh(wavenumber_*Depth_)/Foam::pow(Foam::sinh(wavenumber_*Depth_),3)*
	(2.+Foam::cosh(2*wavenumber_*Depth_)*
	Foam::cos( wavenumber_*xstar(p)-omega_*Time_))
	);
	}
	
scalar Foam::stokes2wave::p(const vector &p)
{
	//only linear approximation for pressure
	return(rho_*maggravity_*(Xi(p)-p.component(2)));
	}
	
void Foam::stokes2wave::setTime(const scalar &Time)
{
	Time_ = Time;
	}


//Validity checks and info  for all monochromatic wave types derived from stokes2wave!!!
void Foam::stokes2wave::checkUrsell()
{

	Info<<"####################  Wave data       ####################### "<<endl;
	Info<<"Waterlevel :	"<<Depth_<<endl;
	Info<<"Period 	  :	"<<Period_<<endl;
	Info<<"WaveHeight :	"<<waveHeight_<<endl;
	Info<<"WaveLength :	"<<lambda_<<endl;
	Info<<"Ursell 	  :	"<<Ursell()<<endl;
	Info<<"rho water  :	"<<rho_<<endl;
	Info<<"Type		  :	 stokes2wave "<<endl;
	Info<<"############################################################# "<<endl;
	
	
	if ( (Ursell() > 10 ) )
	{
		Info<<"WARNING: Ursell parameter is bigger than 10!"<<endl;	
		Info<<"Why don't You rather use stokes5waves?"<<endl;	
		}	
	}
					
	
//***************************Constructor********************************
Foam::stokes2wave::stokes2wave(IOdictionary waveDict, vector gravity, scalar waterdensity):wave(waveDict, gravity)
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

}						

// End namespace Foam
 
// ************************************************************************* //
