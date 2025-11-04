#pragma once

// This is just a link to the Batch Reactor Options, further options could be added here
#include <idealreactors/batch/BatchReactor_Options.h>

namespace BioSMOKE
{
class BioSMOKE_Options : public OpenSMOKE::BatchReactor_Options
{
    // Default constructor, destructor, and all base class methods are taken form the batch class
};
} // namespace BioSMOKE
