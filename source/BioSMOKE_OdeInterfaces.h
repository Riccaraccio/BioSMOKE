#ifndef BioSMOKE_OdeInterfaces_H
#define BioSMOKE_OdeInterfaces_H

#include <math/external-ode-solvers/OpenSMOKE_OdeSystemObject.h>

#include "TGAnalysis.h"

namespace BioSMOKE
{
    class ODESystem_BioSMOKE_TGAnalysis
    {
        public:

            ODESystem_BioSMOKE_TGAnalysis() {};

            void SetReactor(TGAnalysis* reactor)
            {
                reactor_ = reactor;
            }

        protected:

            unsigned int ne_;

            void MemoryAllocation()
            {
                y_.resize(ne_);
                dy_.resize(ne_);
            }

            virtual void Equations(const Eigen::VectorXd &Y, const double t, Eigen::VectorXd &DY)
            {
                y_.assign(Y.data(), Y.data() + Y.size());
                reactor_->Equations(t, y_, dy_);
                DY = Eigen::Map<Eigen::VectorXd>(dy_.data(), dy_.size());
            }

            void Print(const double t, const Eigen::VectorXd &Y)
            {
                y_.assign(Y.data(), Y.data() + Y.size());
                reactor_->Print(t, y_);
            }

            void Jacobian(const Eigen::VectorXd &Y, const double t, Eigen::SparseMatrix<double> &J)
            {
                y_.assign(Y.data(), Y.data() + Y.size());
                reactor_->SparseAnalyticalJacobian(t, y_, J);
            }
            void Jacobian(const Eigen::VectorXd &Y, const double t, Eigen::MatrixXd &J)
            {
                y_.assign(Y.data(), Y.data() + Y.size());
                reactor_->DenseAnalyticalJacobian(t, y_, J);
            }

        private:

            TGAnalysis* reactor_;
            std::vector<double> y_;
            std::vector<double> dy_;
    };
}

#endif // BioSMOKE_OdeInterfaces_H