#ifndef FOAM_INFLOWGENERATORBASE_H
#define FOAM_INFLOWGENERATORBASE_H

#include "autoPtr.H"
#include "fielddataprovider.h"
#include "globalPatch.H"

#include "boost/variant.hpp"

#include "vtkSmartPointer.h"
#include "vtkPolyData.h"
#include "vtkPolyDataWriter.h"
#include "vtkDoubleArray.h"
#include "vtkCellData.h"
#include "vtkMultiBlockDataSet.h"
#include "vtkXMLMultiBlockDataWriter.h"

namespace Foam {




template
<
    class Face,
    template<class> class FaceList,
    class PointField
>
class InflowGeneratorBase
 : public UniPrimitivePatch<Face, FaceList, PointField>
{

public:
    typedef FixedSizeFieldDataProvider<
                scalar,
                UniPrimitivePatch<Face, FaceList, PointField> >
            scalarFixedSizeFieldDataProvider;

    typedef FixedSizeFieldDataProvider<
                vector,
                UniPrimitivePatch<Face, FaceList, PointField> >
            vectorFixedSizeFieldDataProvider;

    typedef FixedSizeFieldDataProvider<
                symmTensor,
                UniPrimitivePatch<Face, FaceList, PointField> >
            symmTensorFixedSizeFieldDataProvider;

private:
    vectorFixedSizeFieldDataProvider UMeanInput_;
    symmTensorFixedSizeFieldDataProvider RInput_;

public:
    typedef vectorFixedSizeFieldDataProvider DirectLengthScaleInput;

    struct WallDistanceAndLengthScaleInput
    {
        scalarFixedSizeFieldDataProvider wallDistance;
        scalarFixedSizeFieldDataProvider lengthScale;
        mutable vectorField L;

        typedef std::tuple<
            const typename scalarFixedSizeFieldDataProvider::input_type&,
            const typename scalarFixedSizeFieldDataProvider::input_type&
          > input_type;
    };

    boost::variant<DirectLengthScaleInput, WallDistanceAndLengthScaleInput> LInput_;

    scalarFixedSizeFieldDataProvider calibrationFactorInput_;
    scalarFixedSizeFieldDataProvider vortonDensityInput_;

    typedef boost::variant<
                const typename DirectLengthScaleInput::input_type&,
                const typename WallDistanceAndLengthScaleInput::input_type&
                > LInput_input_type;

    boost::variant<DirectLengthScaleInput, WallDistanceAndLengthScaleInput> getLInput(
         const LInput_input_type& linput )
    {
        if (const auto* fl =
                boost::get<const typename DirectLengthScaleInput::input_type&>(&linput))
        {
            return DirectLengthScaleInput(*fl, *this);
        }
        else if (const auto* fl =
                boost::get<typename WallDistanceAndLengthScaleInput::input_type>(&linput))
        {
            return WallDistanceAndLengthScaleInput{
                scalarFixedSizeFieldDataProvider(std::get<0>(*fl), *this),
                scalarFixedSizeFieldDataProvider(std::get<1>(*fl), *this)
             };
        }
        else
        {
            FatalErrorIn("InflowGeneratorBase::getLInput")
                    << "internal error: unhandled selection"
                    << abort(FatalError);
        }
    }

protected:
    scalar currentTime_;

protected:
    const polyMesh* mesh_;

public:
    //- Construct from components
    template<class ...Args>
    InflowGeneratorBase
    (
        scalar t0,
        const FieldDataProvider<vector>& UMeanInput,
        const FieldDataProvider<symmTensor>& RInput,
        const LInput_input_type& LInput,
        const FieldDataProvider<scalar>& calibrationFactorInput,
        const FieldDataProvider<scalar>& vortonDensityInput,
        const polyMesh* mesh,
        Args&&... addArgs
    )
        : UniPrimitivePatch<Face, FaceList, PointField>(
              std::forward<Args>(addArgs)... ),
          UMeanInput_(UMeanInput, *this),
          RInput_(RInput, *this),
          LInput_( getLInput(LInput) ),
          calibrationFactorInput_(calibrationFactorInput, *this),
          vortonDensityInput_(vortonDensityInput, *this),
          currentTime_(t0),
          mesh_(mesh)
    {}


    const FieldDataProvider<vector>& UMeanInput() const
    {
        return UMeanInput_.fieldDataProvider();
    }

    const FieldDataProvider<symmTensor>& RInput() const
    {
        return RInput_.fieldDataProvider();
    }

    LInput_input_type LInput() const
    {
        if (const auto* fl =
                boost::get<DirectLengthScaleInput>(&LInput_))
        {
            return fl->fieldDataProvider();
        }
        else if (const auto* fl =
                boost::get<WallDistanceAndLengthScaleInput>(&LInput_))
        {
            return typename WallDistanceAndLengthScaleInput::input_type(
                fl->wallDistance.fieldDataProvider(),
                fl->lengthScale.fieldDataProvider()
            );
        }
        else
        {
            FatalErrorIn("InflowGeneratorBase::LInput()")
                    << "internal error: unhandled selection"
                    << abort(FatalError);
        }
    }


    const FieldDataProvider<scalar>& calibrationFactorInput() const
    {
        return calibrationFactorInput_.fieldDataProvider();
    }

    const FieldDataProvider<scalar>& vortonDensityInput() const
    {
        return vortonDensityInput_.fieldDataProvider();
    }

    const vectorField& UMean(scalar t) const
    {
        return UMeanInput_(t);
    }

    const symmTensorField& R(scalar t) const
    {
        return RInput_(t);
    }

    const vectorField& L(scalar t) const
    {
        if (auto *di = boost::get<DirectLengthScaleInput>(&LInput_))
        {
            return (*di)(t);
        }
        else if (auto *wl = boost::get<WallDistanceAndLengthScaleInput>(&LInput_))
        {
            wl->L =
                    vector(1., 0.216/0.512, 0.064/0.512)
                    * min(0.41*wl->wallDistance(t), wl->lengthScale(t));

            return wl->L;
        }
        else
        {
            FatalErrorIn("InflowGeneratorBase::L(scalar t)")
                    << "internal error: unhandled selection" << endl
                    << abort(FatalError);
        }
    }

    const scalarField& calibrationFactor(scalar t) const
    {
        return calibrationFactorInput_(t);
    }

    const scalarField& vortonDensity(scalar t) const
    {
        return vortonDensityInput_(t);
    }

    virtual void addToVTK(vtkMultiBlockDataSet* ds) const
    {}


    virtual void writeVTK(/*ptf.curTimeIndex_*/
            const fileName& outputFile,
            const std::vector<std::pair<std::string, const Field<scalar>& > >& scalarFields,
            const std::vector<std::pair<std::string, const Field<vector>& > >& vectorFields,
            const std::vector<std::pair<std::string, const Field<symmTensor>& > >& symmTensorFields
            ) const
    {
        auto baseSurface = vtkSmartPointer<vtkPolyData>::New();
        OFPrimitivePatchToVTK(this->localPoints(), this->localFaces(), baseSurface);


        for (const auto& f: scalarFields)
        {
            addFieldToVTK(f, baseSurface->GetCellData());
        }
        for (const auto& f: vectorFields)
        {
            addFieldToVTK(f, baseSurface->GetCellData());
        }
        for (const auto& f: symmTensorFields)
        {
            addFieldToVTK(f, baseSurface->GetCellData());
        }

        auto mbds = vtkSmartPointer<vtkMultiBlockDataSet>::New();
        mbds->SetBlock(0, baseSurface);
        addToVTK(mbds);

        auto w = vtkSmartPointer<vtkXMLMultiBlockDataWriter>::New();
        w->SetFileName(outputFile.c_str());
        w->SetInputData(mbds);
        w->Write();
    }


    virtual void advance(scalar toTime)
    {
        currentTime_=toTime;
    }

    virtual tmp<vectorField> fluctuation() const =0;

    tmp<vectorField> turbulentField() const
    {
        return this->UMean(currentTime_) + fluctuation();
    }
};

} // namespace Foam

#endif // FOAM_INFLOWGENERATORBASE_H
