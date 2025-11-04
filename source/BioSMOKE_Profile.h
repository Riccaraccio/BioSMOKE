#ifndef OpenSMOKE_BioSMOKE_Profile_H
#define OpenSMOKE_BioSMOKE_Profile_H

// This is just a link to the Plug Flow Reactor Profile, further options could be added here
#include <idealreactors/plugflow/PlugFlowReactor_Profile.h>

namespace BioSMOKE
{
class BioSMOKE_Profile : public OpenSMOKE::PlugFlowReactor_Profile
{
  public:
    using OpenSMOKE::PlugFlowReactor_Profile::PlugFlowReactor_Profile;
};
} // namespace BioSMOKE

#endif /* OpenSMOKE_BioSMOKE_Profile_H */