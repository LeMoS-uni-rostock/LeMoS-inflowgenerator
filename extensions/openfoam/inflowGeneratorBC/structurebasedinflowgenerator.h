#ifndef FOAM_STRUCTUREBASEDINFLOWGENERATOR_H
#define FOAM_STRUCTUREBASEDINFLOWGENERATOR_H

#include "globalPatch.H"
#include "recursiveapply.h"

#include "inflowgeneratorbase.h"
#include "vortonqueue.h"

#include "vtkPolyData.h"
#include "vtkPointData.h"

#include "PstreamReduceOps.H"

namespace Foam {

template
<
    class TurbulentStructure,
    class Face,
    template<class> class FaceList,
    class PointField
>
class StructureBasedInflowGenerator
        : public InflowGeneratorBase<Face, FaceList, PointField>
{

protected:
    autoPtr<globalPatch> globalPatch_;

public:
    PtrList<vortonQueue<TurbulentStructure> > queues_;

    StructureBasedInflowGenerator(
            const StructureBasedInflowGenerator& other
            )
        : InflowGeneratorBase<Face, FaceList, PointField>( other )
    {
        queues_.setSize(other.queues_.size());
        forAll(queues_, qi)
        {
            queues_.set(
                        qi,
                        new vortonQueue<TurbulentStructure>(
                            other.queues_[qi]
                            )
                        );
        }
    }

    template<class ...Args>
    StructureBasedInflowGenerator(
            Args&&... addArgs
            )
        : InflowGeneratorBase<Face, FaceList, PointField>(
              std::forward<Args>(addArgs)... )
    {
        queues_.setSize(this->size());
        forAll(*this, fi)
        {
            queues_.set(
                        fi,
                        new vortonQueue<TurbulentStructure>(
                            this->localFaces()[fi],
                            this->localPoints(),
                            this->UMean(this->currentTime_)[fi],
                            this->R(this->currentTime_)[fi],
                            this->L(this->currentTime_)[fi],
                            this->vortonDensity(this->currentTime_)[fi]
                            )
                        );
        }
    }


    void readVortons(Istream& is)
    {
        forAll(*this, fi)
        {
            queues_[fi].read(is);
        }
    }




    Ostream& writeVortons(Ostream& os) const
    {
        forAll(*this, fi)
        {
            if (fi>0) os << endl;
            queues_[fi].write(os);
        }
        return os;
    }




    size_t nVortons() const
    {
        size_t n=0;
        for (auto& q: queues_)
            n+=q.size();
        return n;
    }




    void advance(scalar toTime) override
    {
        label nBefore = nVortons();
        reduce(nBefore, sumOp<label>());

        for (auto& q: queues_)
        {
            q.advance(toTime - this->currentTime_);
        }
        label nAfter = nVortons();
        reduce(nAfter, sumOp<label>());

        InflowGeneratorBase<Face, FaceList, PointField>::advance(toTime);

        Info<<"t="<<this->currentTime_<<": vortons before/after update: "<<nBefore<<"/"<<nAfter<<endl;
    }




    tmp<vectorField> fluctuation() const override
    {
        if (!globalPatch_.valid())
        {
            const_cast<StructureBasedInflowGenerator*>(this)
                -> globalPatch_.reset(
                    new globalPatch(
                        globalPatch::createGlobalPatch
                            <Face, FaceList, PointField>(
                            static_cast<const UniPrimitivePatch<Face, FaceList, PointField>&>(
                                *this ),
                            this->mesh_)(),
                        this->size()
                        )
                    );
        }

        vectorField gfluc(globalPatch_->size(), vector::zero);

        RecursiveApply<TurbulentStructure, globalPatch> apl(globalPatch_(), gfluc);

        labelList n_affected(nVortons(), 0);
        size_t j=0;

        forAll(queues_, lfi)
        {
            auto& q = queues_[lfi];
            for (const auto& vorton: q)
            {
              label gfi = globalPatch_->toGlobalFaceI(lfi); // global face index
//              Pout<<"lfi/gfi="<<lfi<<"/"<<gfi<<endl;
              n_affected[j++] =
                      apl.apply
                      (
                          *vorton,
                          gfi,
                          this->calibrationFactor(this->currentTime_)[lfi],
                          q.properties().influenceLength()
                      );
            }
        }
        Info<<"n_affected: min="<<gMin(n_affected)<<" / max="<<gMax(n_affected)<<" / avg="<<gAverage(n_affected)<<endl;

        // Make fluctuations global
        reduce(gfluc, sumOp<vectorField>());

        return globalPatch_->extractLocalFaceValues(gfluc);
    }




    void addToVTK(vtkMultiBlockDataSet* ds) const override
    {
        auto cloud = vtkSmartPointer<vtkPolyData>::New();
        auto pts = vtkSmartPointer<vtkPoints>::New();
        auto fif = vtkSmartPointer<vtkDoubleArray>::New();
        fif->SetNumberOfComponents(1);
        fif->SetName("faceIndex");
        auto vl = vtkSmartPointer<vtkDoubleArray>::New();
        vl->SetNumberOfComponents(9);
        vl->SetName("L");
        forAll(*this, fi)
        {
            const auto& q = queues_[fi];

            vector l = q.properties().L;
            const vector& L1 = l[0]*q.properties().e[0];
            const vector& L2 = l[1]*q.properties().e[1];
            const vector& L3 = l[2]*q.properties().e[2];
            tensor tt ( L1/ ( l[0]+SMALL ), L2/ ( l[1]+SMALL ), L3/ ( l[2]+SMALL ) ); // vectors are rows!
            tensor Lt = tt.T() &diagTensor( l[0],l[1],l[2] ) &tt;

            for (const auto& vorton: q)
            {
                const point& p = *vorton;
                pts->InsertNextPoint(p.x(), p.y(), p.z());
                fif->InsertNextTuple1(fi);

                vl->InsertNextTuple9(
                        Lt.xx(), Lt.xy(), Lt.xz(),
                        Lt.yx(), Lt.yy(), Lt.yz(),
                        Lt.zx(), Lt.zy(), Lt.zz() );
            }
        }
        cloud->SetPoints(pts);
        cloud->GetPointData()->AddArray(fif);
        cloud->GetPointData()->AddArray(vl);
        ds->SetBlock(1, cloud);


        auto *surfaceDS = static_cast<vtkPolyData*>(ds->GetBlock(0));
        auto disableReason = vtkSmartPointer<vtkIntArray>::New();
        disableReason->SetName("disableReason");
        disableReason->SetNumberOfComponents(1);
        auto minEdgeL = vtkSmartPointer<vtkDoubleArray>::New();
        minEdgeL->SetName("minEdgeLen");
        minEdgeL->SetNumberOfComponents(3);
        for (const auto& q: queues_)
        {
            disableReason->InsertNextTuple1(q.properties().disableReason);
            auto el = q.edgeLength();
            minEdgeL->InsertNextTuple3( el.x(), el.y(), el.z() );
        }
        surfaceDS->GetCellData()->AddArray(disableReason);
        surfaceDS->GetCellData()->AddArray(minEdgeL);
    }

};

} // namespace Foam

#endif // FOAM_STRUCTUREBASEDINFLOWGENERATOR_H
