#include "math/OpenSMOKEVector.h"
#include "BiomassTemperature_Profile.h"

namespace OpenSMOKE
{
	BiomassTemperature_Profile::BiomassTemperature_Profile(OpenSMOKE::OpenSMOKEVectorDouble& x, OpenSMOKE::OpenSMOKEVectorDouble& y)
	{
		x_ = x;
		y_ = y;
		number_of_points_ = x_.Size();

		x0_ = x[1];
		y0_ = y[1];
		xf_ = x[number_of_points_];
		yf_ = y[number_of_points_];

		ChangeDimensions(number_of_points_ - 1, &m_, true);

		for (unsigned int i = 2; i <= number_of_points_; i++)
		{
			m_[i - 1] = (y_[i] - y_[i - 1]) / (x_[i] - x_[i - 1]);
			//std::cout << m_[i - 1] << " " << std::endl;
		}
	}

	double BiomassTemperature_Profile::Get(const double x) const
	{
		if (x0_ > 0 && ((x - x0_) / x0_ > -1.e6))
			OpenSMOKE::FatalErrorMessage("Profile class: the required point is outside the domain (too small)");
		//else if ((x - xf_) / xf_ > 1.e-6)
		//   OpenSMOKE::FatalErrorMessage("Profile class: the required point is outside the domain. Please, increase the interval where the profile is defined");

		if (number_of_points_ == 2)
			return y0_ + m_[1] * (x - x0_);

		for (unsigned int i = 2; i <= number_of_points_; i++)
			if (x <= x_[i])
				return y_[i - 1] + m_[i - 1] * (x - x_[i - 1]);

		// In case of small excess
		return y_[number_of_points_];
	}

	double BiomassTemperature_Profile::Get_increase(const double x) const
	{
		for (unsigned int i = 2; i <= number_of_points_; i++)
			if (x <= x_[i])
				return m_[i - 1];

		// In case of small excess
		return m_[number_of_points_ - 1];
	}

	double BiomassTemperature_Profile::Get_temperature(const double x) const
	{
		for (unsigned int i = 2; i <= number_of_points_; i++)
			if (x <= x_[i])
				return m_[i - 1] * (x - x_[i - 1]) + y_[i - 1];

		// In case of small excess
		int i = number_of_points_;
		return m_[i - 1]* (x - x_[i - 1]) + y_[i - 1];
	}

	double BiomassTemperature_Profile::Get_Tfinal() const
	{
		return y_[number_of_points_];
	}
}

