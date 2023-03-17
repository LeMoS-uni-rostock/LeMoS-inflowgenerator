
#include "base/tools.h"
#include "base/casedirectory.h"
#include "openfoam/ofes.h"
#include "openfoam/openfoamcase.h"
#include "openfoam/solveroutputanalyzer.h"
#include "openfoam/blockmesh_templates.h"
#include "openfoam/caseelements/numerics/unsteadyincompressiblenumerics.h"
#include "openfoam/caseelements/basic/singlephasetransportmodel.h"
#include "openfoam/caseelements/turbulencemodels/smagorinsky_lesmodel.h"
#include "openfoam/caseelements/boundaryconditions/velocityinletbc.h"
#include "openfoam/caseelements/boundaryconditions/pressureoutletbc.h"
#include "openfoam/caseelements/boundaryconditions/wallbc.h"

using namespace insight;

int main(int argc, char*argv[])
{
    try
    {
        insight::assertion(argc==2, "expected exactly one command line argument");

        CaseDirectory wd(false);

        OpenFOAMCase ofc(OFEs::get(argv[1]));

        const std::string inlet="inlet", outlet="outlet";
        int np=4;

        bmd::blockMeshDict_Box::Parameters mp;
        mp.mesh.resolution=
                bmd::blockMeshDict_Box::Parameters::mesh_type::resolution_individual_type
                (1,10,10);
        mp.mesh.XmPatchName=inlet;
        mp.mesh.XpPatchName=outlet;
        auto bm = ofc.insert(new bmd::blockMeshDict_Box(ofc, mp));

        unsteadyIncompressibleNumerics::Parameters nump;
        nump.set_np(np);
        nump.endTime=1.;
        nump.deltaT=1.;
        nump.writeControl=unsteadyIncompressibleNumerics::Parameters::timeStep;
        nump.writeInterval=1.;
        nump.time_integration.timestep_control=
                unsteadyIncompressibleNumerics::Parameters::time_integration_type::timestep_control_fixed_type();
        nump.time_integration.pressure_velocity_coupling=
                unsteadyIncompressibleNumerics::Parameters::time_integration_type::pressure_velocity_coupling_PISO_type(2);
        ofc.insert(new unsteadyIncompressibleNumerics(ofc, nump));

        ofc.insert(new singlePhaseTransportProperties(ofc));
        ofc.insert(new Smagorinsky_LESModel(ofc));

        ofc.createOnDisk(wd);
        ofc.runBlockMesh(wd, bm->nBlocks(), &consoleProgressDisplayer);

        OFDictData::dict boundaryDict;
        ofc.parseBoundaryDict(wd, boundaryDict);

        ofc.insert(new VelocityInletBC(ofc, inlet, boundaryDict));
        ofc.insert(new PressureOutletBC(ofc, outlet, boundaryDict));
        ofc.addRemainingBCs<WallBC>(boundaryDict, WallBC::Parameters() );


        std::shared_ptr<OFdicts> dicts=ofc.createDictionaries();


        auto& cd = dicts->lookupDict("system/controlDict");
        cd.getList("libs").push_back( OFDictData::data("\"libinflowGeneratorBC.so\"") );

        auto& U = dicts->lookupDict("0/U");
        auto& Uinlet = U.subDict("boundaryField").subDict(inlet);
        Uinlet["type"]="inflowGenerator<hatSpot>";
        Uinlet["UMeanSource"]="uniform (1 0 0)";
        Uinlet["RSource"]="uniform (1 0 0 1 0 1)";
        Uinlet["LSource"]="uniform (0.5 0.5 0.5)";
        Uinlet["calibrationFactorSource"]="uniform 1";

        ofc.createOnDisk(wd, dicts);
        ofc.modifyCaseOnDisk ( wd );

        // serial run
        {
            insight::CurrentExceptionContext ex("performing serial run");
            SolverOutputAnalyzer analyzer(consoleProgressDisplayer);
            ofc.runSolver(wd, analyzer, "pimpleFoam");
        }

        // parallel run
        {
            insight::CurrentExceptionContext ex("performing parallel run");
            SolverOutputAnalyzer analyzer(consoleProgressDisplayer);
            ofc.executeCommand(wd, "decomposePar");
            ofc.runSolver(wd, analyzer, "pimpleFoam", np);
            ofc.executeCommand(wd, "reconstructPar", {"-latestTime"});
        }

        return 0;
    }
    catch (std::exception& e)
    {
        std::cerr<<e.what()<<std::endl;
        return -1;
    }
}
