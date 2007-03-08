// Here are the necessary incantations to declare your module to the
// framework, so it can be referenced in a cmsRun file.
//
#include "PluginManager/ModuleDef.h"
#include "FWCore/Framework/interface/MakerMacros.h"

#include "HLTanalyzers/Jets/interface/HLTJetAnalyzer.h"

DEFINE_SEAL_MODULE();
DEFINE_ANOTHER_FWK_MODULE(HLTJetAnalyzer);

