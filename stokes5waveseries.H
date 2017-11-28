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
    Foam::wave::stokes5waveseries

Description
    5th order stokes waves, according to
    * Fenton 1985

SourceFiles
    stokes5waveseries.C

\*---------------------------------------------------------------------------*/
#include "volFields.H"
#include "typeInfo.H"
#include "IOdictionary.H"
#include "scalarListIOList.H"
#include "autoPtr.H"
#include "mathematicalConstants.H"
#include "sinewave.H"
// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
namespace Foam
{
#ifndef stokes5waveseries_H
#define stokes5waveseries_H
class stokes5waveseries: public sinewave
	{
	protected:
	int waveindex_;
	scalar evalwavenumber();
	scalar Q_,U_;

	//Fenton Coeffs
	scalar S_ ,	S2_, S3_, S4_,S5_,S6_,S7_,S8_;
	scalar oneminusS_,oneminusS2_,oneminusS3_,oneminusS4_,oneminusS5_,oneminusS6_;
	scalar 	epsilon_,epsilon2_,	epsilon4_,epsilon3_, epsilon5_;
	scalar	A11_,A22_,A31_,A33_,A42_,A44_,A51_,	A53_,A55_;
	scalar 	B11_ ,B22_,B31_,B33_,B42_,B44_,B51_,B53_,B55_;
	scalarList PeriodList_;
	scalarList WaveHeightList_;
	void checkUrsell();
	void setTime(const scalar & Time);
	inline scalar xstar(vector p){return(p&direction_);};//Get total distance in prop direction 
	
	//Eval coefficients 

	
	//Quick and dirty fix
	inline scalar coth(scalar x)
	{
	return(Foam::cosh(x)/Foam::sinh(x));}
	
	inline scalar A11(){
		return(scalar(1.)/Foam::sinh(wavenumber_*Depth_));};

	inline scalar A22(){
		return(scalar(3.)*S2_/(scalar(2.)*oneminusS2_));};
	
	inline scalar A31(){
		return((-scalar(4.)-scalar(20.)*S_+scalar(10.)*S2_-scalar(13.)*S3_ )/(scalar(8.)*Foam::sinh(wavenumber_*Depth_)*oneminusS3_));};
	
	inline scalar A33(){
		return((-scalar(2.)*S2_+scalar(11.)*S3_)/(scalar(8.)*Foam::sinh(wavenumber_*Depth_)*oneminusS3_));};
	
	inline scalar A42(){
		return((scalar(12.)*S_-scalar(14.)*S2_-scalar(264.)*S3_-scalar(45.)*S4_-scalar(13.)*S5_)/(scalar(24.)*oneminusS5_));};
	
	inline scalar A44(){
		return((10.*S3_-174.*S4_+291.*S5_+278.*S6_)/(48.*(3.+2.*S_)*oneminusS5_));};
	
	inline scalar A51(){
		return((-1184.+32.*S_+13232.*S2_+21712.*S3_+20940.*S4_+12554.*S5_-500.*S6_-3341.*S7_-670.*S8_)/(64.*Foam::sinh(wavenumber_*Depth_)*(3.+2.*S_)*(4.+S_)*(oneminusS6_)) );};
	
	inline scalar A53(){
		return((4.*S_+105.*S2_+198.*S3_-1376.*S4_-1302.*S5_-117.*S6_+58.*S7_)/(32.*Foam::sinh(wavenumber_*Depth_)*(3.+2.*S_)*(oneminusS6_)));};
	
	inline scalar A55(){
		return((-6.*S3_+272.*S4_-1552.*S5_+852.*S6_+2029.*S7_+430.*S8_)/(64.*Foam::sinh(wavenumber_*Depth_)*(3.+2.*S_)*(4.+S_)*oneminusS6_));};
	
	inline scalar B11(){
		return(1.);};
	
	inline scalar B22(){
		return(coth(wavenumber_*Depth_)*(1.+2.*S_)/(2.*oneminusS_));};
		
	inline scalar B31(){
		return(-3.*(1.+3.*S_+3.*S2_+2.*S3_)/(8.*oneminusS3_));};
		
	inline scalar B33(){
		return(3.*(1.+3.*S_+3.*S2_+2.*S3_)/(8.*oneminusS3_));};
		
	inline scalar B42(){
		return(coth(wavenumber_*Depth_)*(6.-26.*S_-182.*S2_-204.*S3_-25.*S4_+26.*S5_)/(6.*(3.+2.*S_)*oneminusS4_));}
	
	inline scalar B44(){
		return(coth(wavenumber_*Depth_)*(24.+92.*S_+122.*S2_+66.*S3_+67.*S4_+34.*S5_)/(24.*(3.+2.*S_)*oneminusS4_));};
		
	inline scalar B51(){
		return(-B53()+B55());};
		
	inline scalar B53(){
		return(9.*(132.+17.*S_-2216.*S2_-5897.*S3_-6292.*S4_-2687.*S5_+194.*S6_+467.*S7_+82.*S8_)/(128.*(3.+2.*S_)*(4.+S_)*oneminusS6_));};
	
	inline scalar B55(){
		return(5.*(300.+1579.*S_-3176.*S2_-2949.*S3_-1188.*S4_-675.*S5_+1326.*S6_+827.*S7_+130.*S8_)/(384.*(3.+2.*S_)*(4.+S_)*oneminusS6_));};
	
	inline scalar E2(){
		return(Foam::tanh(wavenumber_*Depth_)*(2.+2.*S_+5.*S2_)/(4.*oneminusS2_));};
		
	inline scalar E4(){
		return(Foam::tanh(wavenumber_*Depth_)*(8.+12.*S_-152.*S2_-308.*S3_-42.*S4_+77.*S5_)/(32.*oneminusS5_));};
		
	inline scalar Sf(scalar k){
		return(scalar(1.)/Foam::cosh(2.*k*Depth_));};
	
	inline scalar C0(scalar k){
		return(Foam::sqrt(Foam::tanh(k*Depth_)));};
	
	inline scalar C2(scalar k){
		return(Foam::sqrt(Foam::tanh(k*Depth_))*(2.+7.*Foam::pow(Sf(k),2))/(4.*Foam::pow(1-Sf(k),2)));};
	
	inline scalar C4(scalar k){
		return(Foam::sqrt(Foam::tanh(k*Depth_))*(4.+32.*Sf(k)-116.*Foam::pow(Sf(k),2.)-400.*Foam::pow(Sf(k),3.)-71.*Foam::pow(Sf(k),4.)+146.*Foam::pow(Sf(k),5.))/(32.*Foam::pow((1.-Sf(k)),5.)));};
	
	inline scalar D2(scalar k){
		return(-Foam::sqrt(coth(k*Depth_))/2.);};
	
	inline scalar D4(scalar k){
		return(Foam::sqrt(coth(k*Depth_))*(2.+4.*Sf(k)+Foam::pow(Sf(k),4)+2.*Foam::pow(Sf(k),scalar(3.)))/(8.*Foam::pow((1.-Sf(k)),scalar(3.))));};
	
	
	
	inline scalar K(scalar k)
	{
		return(-Foam::constant::mathematical::twoPi/(Period_*Foam::sqrt(maggravity_*k))+C0(k)+Foam::pow((k*waveHeight_/scalar(2.)),scalar(2.))*(C2(k)+D2(k)/(k*Depth_))+Foam::pow((k*waveHeight_/scalar(2.)),scalar(4.))*(C4(k)+D4(k)/(k*Depth_)));
		}
	
	void setCoeffs();
	
	public:		
	vector velocity(const vector &c);
	scalar Xi(const vector & c);
	scalar p(const vector & c);

	        //- Select null constructed
    static autoPtr<stokes5waveseries> New();
	
	stokes5waveseries(IOdictionary waveDict, vector gravity, scalar waterdensity);
	~stokes5waveseries(){};
	};
	
#endif
}
