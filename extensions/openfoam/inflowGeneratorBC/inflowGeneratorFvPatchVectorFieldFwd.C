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

#include "inflowGeneratorFvPatchVectorField.H"
#include "fixedValueFvPatchFields.H"
#include "transform.H"
#include "transformField.H"
#include "volFields.H"
#include "typeInfo.H"
#include "addToRunTimeSelectionTable.H"

#include "turbulentStructures/hatspot.h"
#include "turbulentStructures/anisotropicvorton.h"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#define makeInflowGeneratorFvPatchField(spotType)                   	\
                                                                        \
typedef inflowGeneratorFvPatchVectorField<spotType>                     \
inflowGeneratorFvPatchVectorField##spotType;                            \
                                                                        \
defineTemplateTypeNameAndDebugWithName(                                 \
    inflowGeneratorFvPatchVectorField##spotType,                        \
    "inflowGenerator<"#spotType">", 0);                                 \
                                                                        \
addToRunTimeSelectionTable                                              \
(                                                                       \
    fvPatchVectorField,                                                 \
    inflowGeneratorFvPatchVectorField##spotType,                        \
    patch                                                               \
);                                                                      \
                                                                        \
addToRunTimeSelectionTable                                              \
(                                                                       \
    fvPatchVectorField,                                                 \
    inflowGeneratorFvPatchVectorField##spotType,                        \
    dictionary                                                          \
);                                                                      \
                                                                        \
addToRunTimeSelectionTable                                              \
(                                                                       \
    fvPatchVectorField,                                                 \
    inflowGeneratorFvPatchVectorField##spotType,                        \
    patchMapper                                                         \
)                                                                       \
    

makeInflowGeneratorFvPatchField(hatSpot);


typedef anisotropicVorton<
            PseudoInvsereAnisotropicVortonParametersComputer<
                MaxLengthResolvedByMaxEdge> >
        anisotropicVortonPseudoInverse;

makeInflowGeneratorFvPatchField(anisotropicVortonPseudoInverse);


typedef anisotropicVorton<
            AnalyticAnisotropicVortonParametersComputer<
                MaxLengthResolvedByMaxEdge> >
        anisotropicVortonAnalytic;

makeInflowGeneratorFvPatchField(anisotropicVortonAnalytic);


typedef anisotropicVorton<
            NumericalAnisotropicVortonParametersComputer<
                MaxLengthResolvedByMaxEdge> >
        anisotropicVortonNumerical;

makeInflowGeneratorFvPatchField(anisotropicVortonNumerical);


} // End namespace Foam

// ************************************************************************* //
