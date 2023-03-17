/*
 * This file is part of Insight CAE, a workbench for Computer-Aided Engineering 
 * Copyright (C) 2014  Hannes Kroeger <hannes@kroegeronline.net>
 *
 * This program is free software; you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation; either version 2 of the License, or
 * (at your option) any later version.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License along
 * with this program; if not, write to the Free Software Foundation, Inc.,
 * 51 Franklin Street, Fifth Floor, Boston, MA 02110-1301 USA.
 *
 */


#include "fvCFD.H"
#include "OFstream.H"
#include "transformField.H"
#include "transformGeometricField.H"
#include "fixedGradientFvPatchFields.H"
#include "addToRunTimeSelectionTable.H"
#include "inflowGeneratorFvPatchVectorField.H"
#include "wallDist.H"
#include "interpolationTable.H"

#include <boost/concept_check.hpp>
#include <boost/assign.hpp>


#include "uniof.h"

#include "vtkSmartPointer.h"
#include "vtkUnstructuredGridReader.h"
#include "vtkUnstructuredGrid.h"
#include "vtkPolyDataReader.h"
#include "vtkGenericDataObjectReader.h"
#include "vtkPolyData.h"

#include "structurebasedinflowgenerator.h"
#include "fielddataprovider.h"
#include "turbulentStructures/hatspot.h"
#include "turbulentStructures/anisotropicvorton.h"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

using namespace Foam;


template<class ...Args>
autoPtr<InflowGeneratorBase<face, List, Foam::pointField> >
createStructureBasedInflowGenerator(
        const word& spotTypeName,
        Args&&... addArgs
        )
{
    autoPtr<InflowGeneratorBase<face, List, Foam::pointField> > generator;

    if (spotTypeName=="hatSpot")
    {
        generator.reset(
                new StructureBasedInflowGenerator<hatSpot, face, List, Foam::pointField>(
                      std::forward<Args>(addArgs)...
                        )
                );
    }
    else if (spotTypeName=="anisotropicVortonPseudoInverse")
    {
        typedef anisotropicVorton<
                    PseudoInvsereAnisotropicVortonParametersComputer<
                        MaxLengthResolvedByMaxEdge> >
                Vorton;

        generator.reset(
                new StructureBasedInflowGenerator<Vorton, face, List, Foam::pointField>(
                      std::forward<Args>(addArgs)...
                        )
                );
    }
    else if (spotTypeName=="anisotropicVortonAnalytic")
    {
        typedef anisotropicVorton<
                    AnalyticAnisotropicVortonParametersComputer<
                        MaxLengthResolvedByMaxEdge> >
                Vorton;

        generator.reset(
                new StructureBasedInflowGenerator<Vorton, face, List, Foam::pointField>(
                      std::forward<Args>(addArgs)...
                        )
                );
    }
    else if (spotTypeName=="anisotropicVortonNumerical")
    {
        typedef anisotropicVorton<
                    NumericalAnisotropicVortonParametersComputer<
                        MaxLengthResolvedByMaxEdge> >
                Vorton;

        generator.reset(
                new StructureBasedInflowGenerator<Vorton, face, List, Foam::pointField>(
                      std::forward<Args>(addArgs)...
                        )
                );
    }
    return generator;
}


int main(int argc, char *argv[])
{

    UNIOF_ADDOPT(argList, "source", "string", "how to get generator mesh: may be 'mesh' or 'vtk'");

    UNIOF_ADDOPT(argList, "vtkFile", "string", "path to vtk file with mesh and possibly data");
    UNIOF_ADDOPT(argList, "UMeanSource", "string", "FieldDataProvider string for UMean");
    UNIOF_ADDOPT(argList, "RSource", "string", "FieldDataProvider string for Reynolds stress");
    UNIOF_ADDOPT(argList, "LSource", "string", "FieldDataProvider string for length scale data");
    UNIOF_ADDOPT(argList, "ySource", "string", "FieldDataProvider string for wall distance data. LSource should point to a scalar quantity, if this is provided.");
    UNIOF_ADDOPT(argList, "calibrationFactorSource", "string", "FieldDataProvider string for calibration factor data. calibrationFactorSource should point to a scalar quantity, if this is provided.");
    UNIOF_ADDOPT(argList, "vortonDensitySource", "string", "FieldDataProvider string for vorton density data. vortonDensitySource should point to a scalar quantity, if this is provided.");
    UNIOF_ADDOPT(argList, "outputFile", "string", "prefix of output vtm file");
    UNIOF_ADDOPT(argList, "vortonType", "string", "type of vorton, hatSpot (default), anisotropicVortonPseudoInverse, anisotropicVortonAnalytic");

    UNIOF_ADDOPT(argList, "doFullTimeLoop", "bool", "dry-run fully configured openfoam case");

    UNIOF_ADDOPT(argList, "deltaT", "number", "time step");
    UNIOF_ADDOPT(argList, "writeInterval", "number", "number of time steps between output");
    UNIOF_ADDOPT(argList, "nSteps", "number", "number of time steps to run");

    UNIOF_ADDOPT(argList, "writeCalibrationFactor", "string", "file name for output of calibration factor");

    using Foam::vector; // fix issue in OF5
#   include "setRootCase.H"


    word source="mesh";
    if (args.optionFound("source"))
    {
        source=UNIOF_OPTION(args, "source");
    }

    scalar deltaT = 1.0;
    label writeInterval = 1000;
    label nsteps = 10*writeInterval;
    string outputFile = "inflowGenerator";
    if (args.optionFound("outputFile"))
    {
        outputFile=UNIOF_OPTION(args, "outputFile");
    }

    if (source=="mesh")
    {

#       include "createTime.H"

        bool do_fulltimeloop=false;

        autoPtr<UNIOF_WORDRELIST> patchNames;

        IOobject dicthead
                (
                    "checkInflowGeneratorStatisticsDict",
                    runTime.system(),
                    runTime,
                    IOobject::MUST_READ,
                    IOobject::NO_WRITE
                    );
        if (UNIOF_HEADEROK(dicthead,IOdictionary))
        {
            IOdictionary dict(dicthead);

            if (dict.found("writeInterval"))
                writeInterval=readLabel(dict.lookup("writeInterval"));

            if (dict.found("nSteps"))
                nsteps=readLabel(dict.lookup("nSteps"));

            if (dict.found("doFullTimeLoop"))
                do_fulltimeloop=Switch(dict.lookup("doFullTimeLoop"));

            patchNames.reset(new UNIOF_WORDRELIST(dict.lookup("patches")));
        }

        // command line beats configuration file...
        if (args.optionFound("writeInterval"))
        {
            writeInterval=readLabel(IStringStream(UNIOF_OPTION(args, "writeInterval"))());
        }

        if (args.optionFound("nSteps"))
        {
            nsteps=readLabel(IStringStream(UNIOF_OPTION(args, "nSteps"))());
        }

        if (args.optionFound("doFullTimeLoop"))
        {
            do_fulltimeloop=Switch(IStringStream(UNIOF_OPTION(args, "doFullTimeLoop"))());
        }

#       include "createMesh.H"

        Info << "Reading field U\n" << endl;
        volVectorField U
                (
                    IOobject
                    (
                        "U",
                        runTime.timeName(),
                        mesh,
                        IOobject::MUST_READ,
                        IOobject::AUTO_WRITE
                        ),
                    mesh
                    );

        if (do_fulltimeloop)
        {
            runTime.functionObjects().off();
#           include "createPhi.H"
            while (runTime.run())
            {
                runTime++;
                U.correctBoundaryConditions();
                runTime.write();
            }
        }
        else
        {
            if (args.optionFound("patches"))
            {
                patchNames.reset(new UNIOF_WORDRELIST(IStringStream(UNIOF_OPTION(args, "patches"))()));
            }

            if (!patchNames.valid())
            {
                FatalErrorIn("main")
                        <<"patches to test are not specified!\n"
                       <<"Either provide a list on the command line or in the system/checkInflowGeneratorStatisticsDict!\n"
                      <<abort(FatalError);
            }

            labelHashSet patches =
                    mesh.boundaryMesh().patchSet(patchNames());

            forAllConstIter(labelHashSet, patches, iter)
            {
                label patchI=iter.key();
//                if (isA<inflowGeneratorFvPatchVectorField>(U.boundaryField()[patchI]))
//                {
//                    inflowGeneratorFvPatchVectorField &ifpf =
//                            refCast<inflowGeneratorFvPatchVectorField>( UNIOF_BOUNDARY_NONCONST(U)[patchI] );

//                    ifpf.computeConditioningFactor(writeInterval, nsteps);
//                }
            }
        }


        Info<< "ExecutionTime = " << runTime.elapsedCpuTime() << " s"
            << "  ClockTime = " << runTime.elapsedClockTime() << " s"
            << nl << endl;

    }
    else if (source=="vtk")
    {

        if (args.optionFound("writeInterval"))
        {
            writeInterval=readLabel(IStringStream(UNIOF_OPTION(args, "writeInterval"))());
        }
        if (args.optionFound("nSteps"))
        {
            nsteps=readLabel(IStringStream(UNIOF_OPTION(args, "nSteps"))());
        }
        if (args.optionFound("deltaT"))
        {
            deltaT=readScalar(IStringStream(UNIOF_OPTION(args, "deltaT"))());
        }

        word vortonType="hatSpot";
        if (args.optionFound("vortonType"))
        {
            vortonType = UNIOF_OPTION(args, "vortonType");
        }

        fileName vtkfn(UNIOF_OPTION(args, "vtkFile"));
        string UMeanSource=UNIOF_OPTION(args, "UMeanSource");
        string RSource=UNIOF_OPTION(args, "RSource");
        string LSource=UNIOF_OPTION(args, "LSource");
        string calibrationFactorSource = "uniform 1";
        string vortonDensitySource = "uniform 1";

        if (args.optionFound("calibrationFactorSource"))
            calibrationFactorSource = UNIOF_OPTION(args, "calibrationFactorSource");
        if (args.optionFound("vortonDensitySource"))
            vortonDensitySource = UNIOF_OPTION(args, "vortonDensitySource");

        //auto r = vtkSmartPointer<vtkUnstructuredGridReader>::New();
//        auto r = vtkSmartPointer<vtkPolyDataReader>::New();
        auto r = vtkSmartPointer<vtkGenericDataObjectReader>::New();
        r->SetFileName(vtkfn.c_str());
        r->Update();

        Foam::pointField pts;
        faceList faces;
        if (r->IsFilePolyData())
        {
            VTKMeshToOF(
                        r->GetPolyDataOutput(),
                        pts, faces );
        }
        else if (r->IsFileUnstructuredGrid())
        {
            VTKMeshToOF(
                        r->GetUnstructuredGridOutput(),
                        pts, faces );
        }
        else
        {
            FatalErrorIn("main")
                    <<"Unhandled grid VTK file type!\n"
                  <<abort(FatalError);
        }


        scalar T=0.0;

        auto UMeanInput = FieldDataProvider<vector>::New(
                    IStringStream(UMeanSource)());

        auto RInput = FieldDataProvider<symmTensor>::New(
                    IStringStream(RSource)());

        auto calibrationFactorInput = FieldDataProvider<scalar>::New(
                    IStringStream(calibrationFactorSource)());

        auto vortonDensityInput = FieldDataProvider<scalar>::New(
                    IStringStream(vortonDensitySource)());

        typedef InflowGeneratorBase<face, List, Foam::pointField> Generator;
        autoPtr<Generator> generator;

        if (args.optionFound("ySource"))
        {
            auto LInput = FieldDataProvider<scalar>::New(
                        IStringStream(LSource)());
            auto yInput = FieldDataProvider<scalar>::New(
                        IStringStream(UNIOF_OPTION(args, "ySource"))());

            generator=createStructureBasedInflowGenerator(
                    vortonType,
                    T,
                    UMeanInput(), RInput(),
                    Generator::WallDistanceAndLengthScaleInput::input_type
                                { yInput(), LInput() },
                    calibrationFactorInput(),
                    vortonDensityInput(),
                    nullptr,
                    faces, pts);
        }
        else
        {
            auto LInput = FieldDataProvider<vector>::New(
                        IStringStream(LSource)());

            generator=createStructureBasedInflowGenerator(
                    vortonType,
                    T,
                    UMeanInput(), RInput(), LInput(),
                    calibrationFactorInput(),
                    vortonDensityInput(),
                    nullptr,
                    faces, pts);
        }


        Foam::vectorField UMeanResult(generator->size(), pTraits<vector>::zero);
        Foam::symmTensorField RResult(generator->size(), pTraits<symmTensor>::zero);

        for (long int i=0; i<=nsteps; ++i)
        {
            T+=deltaT;
            generator->advance(T);

            Foam::vectorField u( generator->turbulentField() );
            Foam::vectorField uPrime( u - generator->UMean(T) );

            UMeanResult += u*deltaT;
            RResult += sqr(uPrime)*deltaT;
            Foam::symmTensorField R( RResult/T );

            scalar relax=0.7;
            Foam::scalarField calib(
                        generator->calibrationFactor(T) // currently in use
                        *
                        ( 1. - relax*( 1. - Foam::tr(generator->R(T))/Foam::tr(R) )) // desired
                        );

            if (i % writeInterval == 0)
            {
                generator->writeVTK(
                            str(boost::format(outputFile+"_%d.vtm") % i),
                            {
                                {"calibrationFactor", calib },
                                {"vortonDensity", generator->vortonDensity(T) }
                            },
                            {
                                { "u", u },
                                { "uPrime", uPrime },
                                { "UMeanInput", generator->UMean(T) },
                                { "LInput", generator->L(T) },
                                { "UMean", UMeanResult/T }
                            },
                            {
                                { "RInput", generator->R(T) },
                                { "R", RResult/T }
                            }
                            );

                if (args.optionFound("writeCalibrationFactor"))
                {
                    OFstream outf(UNIOF_OPTION(args, "writeCalibrationFactor"));
                    outf<<calib;
                }
            }

            Info<<"R mean/max/min = "<<average(R)<<"/"<<max(R)<<"/"<<min(R)<<endl;
        }
    }


    Info<< "End\n" << endl;

    return 0;
}

// ************************************************************************* //
