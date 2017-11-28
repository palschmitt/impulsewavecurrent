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
    Foam::wave::cnoidal2wave

Description
    Linear theory waves

SourceFiles
    cnoidal2wave.C

\*---------------------------------------------------------------------------*/
#include "volFields.H"
#include "typeInfo.H"
#include "IOdictionary.H"
#include "autoPtr.H"
#include "mathematicalConstants.H"
#include "wave.H"
// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
namespace Foam
{
#ifndef cnoidal2wave_H
#define cnoidal2wave_H
class cnoidal2wave: public wave
	{
	protected:
	scalar waveHeight_;
	scalar wavenumber_;
	scalar lambda_;
	scalar omega_;
	scalar celerity_;
	scalar rho_;
	scalar Period_;
	scalar evalwavenumber();
	
	inline scalar Ursell(){return(waveHeight_*pow(lambda_,2)/pow(Depth_,3));}
	void checkUrsell();
	inline scalar xstar(vector p){return(p&direction_);};//Get total distance in prop direction 
	
	public:		
	vector velocity(const vector &c);
	scalar Xi(const vector & c);
	scalar p(const vector & c);
	void setTime(const scalar & Time);
	
	        //- Select null constructed
    static autoPtr<cnoidal2wave> New();
	
	cnoidal2wave(IOdictionary waveDict, vector gravity, scalar waterdensity);
	cnoidal2wave(scalar waterdensity);//to be called from derived classes
	~cnoidal2wave(){};
	};
#endif
}
