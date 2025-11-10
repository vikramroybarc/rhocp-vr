// File for Coupled Raccoon and RhoCp Simulations

#include "RhocpApp.h"
#include "Moose.h"
#include "AppFactory.h"
#include "ModulesApp.h"
#include "MooseSyntax.h"

// Try to include Raccoon's app header if available
#include "raccoonApp.h"
#define HAVE_RACCOON_APP 1

InputParameters
RhocpApp::validParams()
{
  InputParameters params = MooseApp::validParams();
  params.set<bool>("use_legacy_material_output") = false;
  params.set<bool>("use_legacy_initial_residual_evaluation_behavior") = false;
  return params;
}

RhocpApp::RhocpApp(const InputParameters & parameters) : MooseApp(parameters)
{
  RhocpApp::registerAll(_factory, _action_factory, _syntax);
}

RhocpApp::~RhocpApp() {}

void
RhocpApp::registerAll(Factory & f, ActionFactory & af, Syntax & syntax)
{
  // Register all MOOSE module objects under this app
  ModulesApp::registerAllObjects<RhocpApp>(f, af, syntax);

  // Register this app's objects/actions
  Registry::registerObjectsTo(f, {"RhocpApp"});
  Registry::registerActionsTo(af, {"RhocpApp"});

  // --- Bring in Raccoon's objects/actions when linked as a library ---
#ifdef HAVE_RACCOON_APP
  // Option A: Call Raccoon's own registerAll (ensures any custom syntax/flags are added)
  raccoonApp::registerAll(f, af, syntax);
#else
  // Option B (header not found): If you still link Raccoon as a lib, this maps by app-name.
  // Leave these two lines if you want to rely on the registry string mapping only.
  Registry::registerObjectsTo(f, {"raccoonApp"});
  Registry::registerActionsTo(af, {"raccoonApp"});
#endif

  /* register custom execute flags, action syntax, etc. here */
}

void
RhocpApp::registerApps()
{
  // Always register this app
  registerApp(RhocpApp);

  // Also register the Raccoon application if its header is available
#ifdef HAVE_RACCOON_APP
  raccoonApp::registerApps();
#endif
}

/***************************************************************************************************
 *********************** Dynamic Library Entry Points - DO NOT MODIFY ******************************
 **************************************************************************************************/
extern "C" void
RhocpApp__registerAll(Factory & f, ActionFactory & af, Syntax & s)
{
  RhocpApp::registerAll(f, af, s);
}

extern "C" void
RhocpApp__registerApps()
{
  RhocpApp::registerApps();
}
