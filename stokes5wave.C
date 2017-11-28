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
    Foam::stokes5wave::sinestokes5wave
    virtual stokes5wave class, different theories should be implemented in derivations
\*---------------------------------------------------------------------------*/
#include "stokes5wave.H"

namespace Foam
 {
  ////* * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

scalar Foam::stokes5wave::evalwavenumber()
{
	
	scalar omegaest = Foam::constant::mathematical::twoPi/Period_;
	scalar kestimate = Foam::pow(omegaest,scalar(2.))/maggravity_*
					   Foam::pow( Foam::pow( coth(omegaest*Foam::sqrt(Depth_/maggravity_ )),scalar(3./2.) ),scalar(2./3.));
	//print 'kestimate ',kestimate
	scalar right = kestimate*1.3;
	scalar left = kestimate*0.7;
	scalar rl,rr,middle,rm;
	while (fabs(left-right) > 0.001)
	{
		rl = K(left);
		rr = K(right);
		middle = 0.5*(left+right);
		rm = K(middle);
		//print 'rm ',rm
		if (rm*rr < 0.0){
			left = middle;
		}
		else
		{
			right = middle;
		}
		if(rr*rl > 0.)
		{
			FatalErrorIn("Foam::stokes5wave::evalwavenumber")
			<< "ERROR: Iteration of k failed! This is either a severe bug or a very unsuitbale use of stokes theory! " << middle
			<< abort(FatalError);
		}
	}
	return(middle);
}


vector Foam::stokes5wave::velocity(const vector & p)
{
	scalar  u = celerity_ - U_ + Foam::sqrt(maggravity_/Foam::pow(wavenumber_,scalar(3.)))*C0(wavenumber_)*
		(epsilon5_*
		(5.*A55_*wavenumber_*Foam::cos(5.*wavenumber_*(xstar(p)-celerity_*Time_))*Foam::cosh(5.*wavenumber_*p.component(2))+
		3.*A53_*wavenumber_*Foam::cos(3.*wavenumber_*(xstar(p)-celerity_*Time_))*Foam::cosh(3.*wavenumber_*p.component(2))+
		A51_*wavenumber_*Foam::cos(wavenumber_*(xstar(p)-celerity_*Time_))*Foam::cosh(wavenumber_*p.component(2)))+
		epsilon4_*
		(4.*A44_*wavenumber_*Foam::cos(4.*wavenumber_*(xstar(p)-celerity_*Time_))*Foam::cosh(4.*wavenumber_*p.component(2))+
		2.*A42_*wavenumber_*Foam::cos(2.*wavenumber_*(xstar(p)-celerity_*Time_))*Foam::cosh(2.*wavenumber_*p.component(2)))+
		epsilon3_*
		(3.*A33_*wavenumber_*Foam::cos(3.*wavenumber_*(xstar(p)-celerity_*Time_))*Foam::cosh(3.*wavenumber_*p.component(2))+
		A31_*wavenumber_*Foam::cos(wavenumber_*(xstar(p)-celerity_*Time_))*Foam::cosh(wavenumber_*p.component(2)))+
		epsilon2_*
		(2.*A22_*wavenumber_*Foam::cos(2.*wavenumber_*(xstar(p)-celerity_*Time_))*Foam::cosh(2.*wavenumber_*p.component(2)))+
		epsilon_*
		A11_*wavenumber_*Foam::cos(wavenumber_*(xstar(p)-celerity_*Time_))*Foam::cosh(wavenumber_*p.component(2)));

	scalar w = Foam::sqrt(maggravity_/Foam::pow(wavenumber_,scalar(3.)))*C0(wavenumber_)*
		(epsilon5_*
		(5.*A55_*wavenumber_*Foam::sin(5.*wavenumber_*(xstar(p)-celerity_*Time_))*Foam::sinh(5.*wavenumber_*p.component(2))+
		3.*A53_*wavenumber_*Foam::sin(3.*wavenumber_*(xstar(p)-celerity_*Time_))*Foam::sinh(3.*wavenumber_*p.component(2))+
		A51_*wavenumber_*Foam::sin(wavenumber_*(xstar(p)-celerity_*Time_))*Foam::sinh(wavenumber_*p.component(2)))+
		epsilon4_*
		(4.*A44_*wavenumber_*Foam::sin(4.*wavenumber_*(xstar(p)-celerity_*Time_))*Foam::sinh(4.*wavenumber_*p.component(2))+
		2.*A42_*wavenumber_*Foam::sin(2.*wavenumber_*(xstar(p)-celerity_*Time_))*Foam::sinh(2.*wavenumber_*p.component(2)))+
		epsilon3_*
		(3.*A33_*wavenumber_*Foam::sin(3.*wavenumber_*(xstar(p)-celerity_*Time_))*Foam::sinh(3.*wavenumber_*p.component(2))+
		A31_*wavenumber_*Foam::sin(wavenumber_*(xstar(p)-celerity_*Time_))*Foam::sinh(wavenumber_*p.component(2)))+
		epsilon2_*
		(2.*A22_*wavenumber_*Foam::sin(2.*wavenumber_*(xstar(p)-celerity_*Time_))*Foam::sinh(2.*wavenumber_*p.component(2)))+
		epsilon_*
		A11_*wavenumber_*Foam::sin(wavenumber_*(xstar(p)-celerity_*Time_))*Foam::sinh(wavenumber_*p.component(2)));


	return(direction_*u-gravityUnitV_*w);
	}


scalar Foam::stokes5wave::Xi(const vector &p)
{
	return(	Depth_+(
	B55_*epsilon5_*Foam::cos(5.*wavenumber_*xstar(p)-5.*celerity_*wavenumber_*Time_)+
	(B44_*epsilon4_)*Foam::cos(4*wavenumber_*xstar(p)-4*celerity_*wavenumber_*Time_)+
	(B53_*epsilon5_+B33_*epsilon3_*Foam::cos(3.*wavenumber_*xstar(p)-3.*celerity_*wavenumber_*Time_)+
	(B42_*epsilon4_+B22_*epsilon2_*Foam::cos(2.*wavenumber_*xstar(p)-2.*celerity_*wavenumber_*Time_)+
	(B51_*epsilon5_+B31_*epsilon3_+B11_*epsilon_)*Foam::cos(wavenumber_*xstar(p)-celerity_*wavenumber_*Time_))/(wavenumber_))));
	}
	
scalar Foam::stokes5wave::p(const vector &p)
{
	//only linear approximation for pressure
	return(rho_*maggravity_*(Xi(p)-p.component(2)));
	}
	

//Validity checks and info  for all monochromatic wave types derived from stokes5wave!!!
void Foam::stokes5wave::checkUrsell()
{

	Info<<"####################  Wave data       ####################### "<<endl;
	Info<<"Waterlevel :	"<<Depth_<<endl;
	Info<<"Period 	  :	"<<Period_<<endl;
	Info<<"WaveHeight :	"<<waveHeight_<<endl;
	Info<<"WaveLength :	"<<lambda_<<endl;
	Info<<"Ursell 	  :	"<<Ursell()<<endl;
	Info<<"rho water  :	"<<rho_<<endl;
	Info<<"Type  	  :	stokes5wave	"<<endl;
	Info<<"############################################################# "<<endl;
	
	
	if ( (Ursell() > 40 ) )
	{
		Info<<"WARNING: Ursell parameter is bigger than 40!"<<endl;	
		Info<<"Why don't You use cnoidal wave theory?"<<endl;	
		}	
	}

void Foam::stokes5wave::setCoeffs()
{

	S_ = Sf(wavenumber_);
	S2_ = Foam::pow(S_,scalar(2.));
	S3_ = Foam::pow(S_,scalar(3.));
	S4_ = Foam::pow(S_,scalar(4.));
	S5_ = Foam::pow(S_,scalar(5.));
	S6_ = Foam::pow(S_,scalar(6.));
	S7_ = Foam::pow(S_,scalar(7.));
	S8_ = Foam::pow(S_,scalar(8.));
	oneminusS_ = scalar(1.)-S_;
	oneminusS2_ = Foam::pow(oneminusS_,scalar(2.));
	oneminusS3_ = Foam::pow(oneminusS_,scalar(3.));
	oneminusS4_ = Foam::pow(oneminusS_,scalar(4.));
	oneminusS5_ = Foam::pow(oneminusS_,scalar(5.));
	oneminusS6_ = Foam::pow(oneminusS_,scalar(6.));
	
	epsilon_ = wavenumber_*waveHeight_/scalar(2.);
	epsilon2_ = Foam::pow(epsilon_,scalar(2.));
	epsilon3_ = Foam::pow(epsilon_,scalar(3.));
	epsilon4_ = Foam::pow(epsilon_,scalar(4.));
	epsilon5_ = Foam::pow(epsilon_,scalar(5.));
	
	A11_ = A11();
	A22_ = A22();
	A31_ = A31();
	A33_ = A33();
	A42_ = A42();
	A44_ = A44();
	A51_ = A51();
	A53_ = A53();
	A55_ = A55();
	B11_ = B11();
	B22_ = B22();
	B31_ = B31();
	B33_ = B33();
	B42_ = B42();
	B44_ = B44();
	B51_ = B51();
	B53_ = B53();
	B55_ = B55();
}					
	
//***************************Constructor********************************
Foam::stokes5wave::stokes5wave(IOdictionary waveDict, vector gravity, scalar waterdensity):sinewave( waterdensity)
{
	waveHeight_= readScalar(waveDict.lookup("WaveHeight")); 
	vector dir  
	(
		waveDict.lookup("direction")
	);
	direction_ = dir/mag(dir);
	maggravity_= mag(gravity);
	gravityUnitV_ = gravity/maggravity_;//Normalise gravity vector 
	
	Period_= readScalar(waveDict.lookup("Period"));
	Depth_= readScalar(waveDict.lookup("waterLevel")); 
	omega_	=  Foam::constant::mathematical::twoPi/Period_;	
	wavenumber_ = evalwavenumber();
	setCoeffs();
	
	lambda_ = Foam::constant::mathematical::twoPi/wavenumber_;
	U_ = (1./Foam::sqrt(wavenumber_/maggravity_))*(C0(wavenumber_)+epsilon2_*C2(wavenumber_)+epsilon4_*C4(wavenumber_));
	Q_ = (1./Foam::sqrt(Foam::pow(wavenumber_,scalar(3.))/maggravity_))*(C0(wavenumber_)*wavenumber_*Depth_+epsilon2_*(C2(wavenumber_)*wavenumber_*Depth_+D2(wavenumber_))+epsilon4_*(C4(wavenumber_)*wavenumber_*Depth_+D4(wavenumber_)));
	celerity_ = Q_/Depth_;
	checkUrsell();
	}

}						




// End namespace Foam
 
// ************************************************************************* //
