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
    Foam::genericwave::sinegenericwave
    virtual genericwave class, different theories should be implemented in derivations
\*---------------------------------------------------------------------------*/
#include "genericwave.H"

namespace Foam
 {
  ////* * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //




vector Foam::genericwave::velocity(const vector & p)
{
	vector U;
	if ( p.component(2) <= Depth_)
	{
		U = direction_*imposedvel_*Foam::sin(omega_*Time_); 
	}
	else
	{
		U = vector(0.,0.,0.);
		}
	return(U);
	}


void Foam::genericwave::setTime(const scalar &Time)
{
	Time_ = Time;
	}



					
	
//***************************Constructor********************************
Foam::genericwave::genericwave(IOdictionary waveDict, vector gravity, scalar waterdensity):wave(waveDict, gravity)
{
	rho_ = waterdensity;	
	imposedvel_= readScalar(waveDict.lookup("imposedVel")); 
	Period_= readScalar(waveDict.lookup("Period")); 
	omega_	=  Foam::constant::mathematical::twoPi/Period_;	
	}

scalar Foam::genericwave::Xi(const vector &p)
{
	return(Depth_);
	}			

Foam::genericwave::genericwave(scalar waterdensity):wave(waterdensity)
{
	rho_ = waterdensity;	
	}

}

// End namespace Foam
 
// ************************************************************************* //
