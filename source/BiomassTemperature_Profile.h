
#ifndef OpenSMOKE_BiomassTemperature_Profile_H
#define	OpenSMOKE_BiomassTemperature_Profile_H

#include "math/OpenSMOKEVector.h"

namespace OpenSMOKE
{

	//!  A class for setting profiles along a plug flow reactor
	/*!
	The purpose of this class is to set a profile along a plug flow reactor
	*/

	class BiomassTemperature_Profile
	{
	public:

		BiomassTemperature_Profile(OpenSMOKE::OpenSMOKEVectorDouble& x, OpenSMOKE::OpenSMOKEVectorDouble& y);

		double Get(const double x) const;
		double Get_increase(const double x) const;		
		double Get_temperature(const double x) const;
		double Get_Tfinal() const;

	private:

		unsigned int number_of_points_;
		bool time_independent_;
		double x0_;
		double xf_;
		double y0_;
		double yf_;

		OpenSMOKE::OpenSMOKEVectorDouble x_;
		OpenSMOKE::OpenSMOKEVectorDouble y_;
		OpenSMOKE::OpenSMOKEVectorDouble m_;
	};
}

#include "BiomassTemperature_Profile.hpp"

#endif	/* OpenSMOKE_BiomassTemperature_Profile_H */

