/*----------------------------------------------------------------------*\
|    ___                   ____  __  __  ___  _  _______                  |
|   / _ \ _ __   ___ _ __ / ___||  \/  |/ _ \| |/ / ____| _     _         |
|  | | | | '_ \ / _ \ '_ \\___ \| |\/| | | | | ' /|  _| _| |_ _| |_       |
|  | |_| | |_) |  __/ | | |___) | |  | | |_| | . \| |__|_   _|_   _|      |
|   \___/| .__/ \___|_| |_|____/|_|  |_|\___/|_|\_\_____||_|   |_|        |
|        |_|                                                              |
|                                                                         |
|   Author: Alberto Cuoci <alberto.cuoci@polimi.it>                       |
|   CRECK Modeling Group <http://creckmodeling.chem.polimi.it>            |
|   Department of Chemistry, Materials and Chemical Engineering           |
|   Politecnico di Milano                                                 |
|   P.zza Leonardo da Vinci 32, 20133 Milano                              |
|                                                                         |
|-------------------------------------------------------------------------|
|                                                                         |
|   This file is part of OpenSMOKE++ framework.                           |
|                                                                         |
|   License                                                               |
|                                                                         |
|   Copyright(C) 2014, 2013, 2012  Alberto Cuoci                          |
|   OpenSMOKE++ is free software: you can redistribute it and/or modify   |
|   it under the terms of the GNU General Public License as published by  |
|   the Free Software Foundation, either version 3 of the License, or     |
|   (at your option) any later version.                                   |
|                                                                         |
|   OpenSMOKE++ is distributed in the hope that it will be useful,        |
|   but WITHOUT ANY WARRANTY; without even the implied warranty of        |
|   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the         |
|   GNU General Public License for more details.                          |
|                                                                         |
|   You should have received a copy of the GNU General Public License     |
|   along with OpenSMOKE++. If not, see <http://www.gnu.org/licenses/>.   |
|                                                                         |
\*-----------------------------------------------------------------------*/

#ifndef OpenSMOKE_ODESystem_Interface_H
#define OpenSMOKE_ODESystem_Interface_H

#include "math/OpenSMOKEVector.h"

typedef void(*odefunction)(const double* y, const double t, double* dy, void* args);
typedef void(*printfunction)(const double t, const Eigen::VectorXd& Y);

namespace OpenSMOKE 
{
  class ODESystem_Interface
  {
  public:

    //ODESystem_Interface() { args_ = NULL; }
    ODESystem_Interface() {}

    ~ODESystem_Interface () {}

    void SetSystemOfEquations(odefunction odefun)
    {
        odefun_ = odefun;
    }

    void SetNumberOfEquations(unsigned int ne)
    {
      ne_ = ne;
    }

    void SetPrintFunction(printfunction printfun)
    {
        printfun_ = printfun;
    }

    void SetUserArgs(void * args)
    {
      args_ = args;
    }

  
  protected:

    unsigned int ne_;

    void MemoryAllocation ()
    {
    }

    virtual void Equations(const Eigen::VectorXd &Y, const double t, Eigen::VectorXd &DY)
    {
      odefun_(Y.data(), t, DY.data(), args_);
    }

    void Jacobian(const Eigen::VectorXd &Y, const double t, Eigen::MatrixXd &J)
    {
    }

    void Jacobian(const Eigen::VectorXd &Y, const double t, Eigen::SparseMatrix<double> &J)
    {
    }

    void Print(const double t, const Eigen::VectorXd &Y)
    {
        printfun_(t, Y);
    }


  private:

    void * args_ = NULL;
    odefunction odefun_ = NULL;
    printfunction printfun_ = NULL;
  };
}

#endif  // ODESystem_Interface
