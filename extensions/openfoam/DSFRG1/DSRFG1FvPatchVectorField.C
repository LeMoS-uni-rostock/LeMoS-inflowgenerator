#define PI 3.14159
#include "DSRFG1FvPatchVectorField.H"
#include "DLList.H"
#include "point.H"
#include "Time.H"

namespace Foam
{


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //
	//- Construct from patch and internal field


DSRFG1FvPatchVectorField::DSRFG1FvPatchVectorField
(
  const fvPatch& p,
  const DimensionedField<vector, volMesh>& iF
)
:
fixedValueFvPatchField<vector>(p,iF),
rand_(label(0)),
SpectrumParts_(10),
NModes_(300),
curTimeIndex_(-1),
DataCalIndex_(label(0))
{
}


//- Construct from patch, internal field and dictionary
	DSRFG1FvPatchVectorField::DSRFG1FvPatchVectorField
	(
    		const fvPatch& p,
    		const DimensionedField<vector, volMesh>& iF,
    		const dictionary& dict
	)
	:
	fixedValueFvPatchField<vector>(p,iF),
	rand_(label(0)),
	SpectrumParts_(readScalar(dict.lookup("SpectrumParts"))),
	NModes_(readScalar(dict.lookup("NModes"))),
	Ufile(dict.lookup("U")),
	curTimeIndex_(-1),
	Lfile(dict.lookup("L")),
	Rfile(dict.lookup("R")),
	Theta1(readScalar(dict.lookup("Theta1"))),
        Theta2(readScalar(dict.lookup("Theta2"))),
	DataCalIndex_(label(0))
	{

		if (dict.found("value"))
    		{
        		fixedValueFvPatchField<vector>::operator==
        		(
            		  Field<vector>("value", dict, p.size())
        		);
    		}

	}



	//- Construct by mapping given DSRFG1FvPatchVectorField onto a new patch
	DSRFG1FvPatchVectorField::DSRFG1FvPatchVectorField
         (
             const DSRFG1FvPatchVectorField& ptf,
             const fvPatch& p,
             const DimensionedField<vector, volMesh>& iF,
             const fvPatchFieldMapper& mapper
         )
	 :
	 fixedValueFvPatchField<vector>(ptf,p,iF,mapper),
	 rand_(label(0)),
	 SpectrumParts_(ptf.SpectrumParts_),
	 NModes_(ptf.NModes_),
	 Ufile(ptf.Ufile),
	 curTimeIndex_(-1),
	 Lfile(ptf.Lfile),
	 Rfile(ptf.Rfile),
	 Theta1(ptf.Theta1),
         Theta2(ptf.Theta2),
	 DataCalIndex_(label(0))
	 {}


	//-construct as copy
	DSRFG1FvPatchVectorField::DSRFG1FvPatchVectorField
         (
             const DSRFG1FvPatchVectorField& ptf
         )
	 :
	 fixedValueFvPatchField<vector>(ptf),
	 rand_(label(0)),
	 SpectrumParts_(ptf.SpectrumParts_),
	 NModes_(ptf.NModes_),
	 Ufile(ptf.Ufile),
	 curTimeIndex_(-1),
	 Lfile(ptf.Lfile),
	 Rfile(ptf.Rfile),
	 Theta1(ptf.Theta1),
         Theta2(ptf.Theta2),
	 DataCalIndex_(label(0))
	 {}


	//- Construct as copy setting internal field reference
	DSRFG1FvPatchVectorField::DSRFG1FvPatchVectorField
         (
            const DSRFG1FvPatchVectorField& ptf,
             const DimensionedField<vector, volMesh>& iF
        )
	:
	fixedValueFvPatchField<vector>(ptf,iF),
	rand_(label(0)),
	SpectrumParts_(ptf.SpectrumParts_),
	NModes_(ptf.NModes_),
	Ufile(ptf.Ufile),
	curTimeIndex_(-1),
	Lfile(ptf.Lfile),
	Rfile(ptf.Rfile),
	Theta1(ptf.Theta1),
        Theta2(ptf.Theta2),
	DataCalIndex_(label(0))
	{}

// * * * * * * * * Definitions of static Members  * * * * * * * * * * * //

	 const scalar DSRFG1FvPatchVectorField::C_L = 3.0;
    	 const scalar DSRFG1FvPatchVectorField::C_tau = 2.0;
   	 const scalar DSRFG1FvPatchVectorField::Beta = 0.5;
    	 const scalar DSRFG1FvPatchVectorField::Alpha = 0.03;
    	 const scalar DSRFG1FvPatchVectorField::L_eta = 0.00001;


// * * * * * * * * * * * Member Functions * * * * * * * * * * * * * * * * * //

// "Calculating Deltay and Deltaz and Distance ..."
void DSRFG1FvPatchVectorField::CalculateDeltas()
{

	const polyPatch& Ppatch = patch().patch(); // alias to the polypatch class (Object)
	const labelListList& FaceEdges = Ppatch.faceEdges(); // list of labels for faces(inlet patch) and correspoding labels for edges
	const edgeList& Edges = Ppatch.edges(); // list of edges, edge is a class whose private members are of class 'point'.
	const Field<point>& lpoints = Ppatch.localPoints(); // list of points, point is a typedef of vector class stores 3 scalars as a private members

 // Declaration of deltaz, deltay, calling size() member function inherited from fvPatch class
 // size() function gives the number of faces for a patch, Here our patch is inlet.
	deltay = scalarField(size(),0);
	deltaz = scalarField(size(),0);

		//scalar ymin;
// Declaration of distance from the bottom wall to the its face center
	Distancey = scalarField(size(),0);
	point pMin; // closet point from the wall
///////////////////////////////////////////////////////////////
// We are calculating deltay and deltaz for each face of the inlet patch
// Logic:
//-------
// 'FaceEdge' is of type labelListList, each face for a patch has a label(number), Each List(array) number is an address to the face label(number)
// Each List(array) index points to a list of edge labels. Therefore finally we are getting label list of edges for a corresponding face.
// In inner loop, edge magnitude is calculated.
// Idea: In a yz plane, The y coordinates are different at starting and end points of the face edge, difference is assignes to deltay.
// Same procedure follows for deltaz
// CurdeltayIndex, CurdeltazIndex are the index not to calculate deltay or deltaz twice for a face, because face has 4 edges.
//drawback 1:
//--------
// The inlet need to be in yz plane
//solution for drawback 1:
//-----------------------
// Change lpoints[Start].z(), lpoints[Start].y() according to your plane.
//second solution, add extra constraint any two coordinates are same for a edge of a face on a inlet patch (may not work inclined inlet patches)
// drawback 2:
//-----------
// Doesnot work for circular cross-section inlet for curve edges
//solution for drawback 2:
//-----------------------
// Need to find a way to get lenghts of  curved edges, rest is fine.

	forAll(FaceEdges,i) // forAll is the macro for the for loop
	{

		label CurdeltayIndex =0;
		label CurdeltazIndex =0;

   	   forAll(FaceEdges[i],j)
   	   {
		const label& L = FaceEdges[i][j];
		const edge& e = Edges[L];
		const label& Start = e.start();
		const label& End = e.end();


		if(lpoints[Start].y() != lpoints[End].y() && CurdeltayIndex == 0)
		{
		  deltay[i] = mag(lpoints[End].y()-lpoints[Start].y());
		  CurdeltayIndex = 1;
		}
		else if(lpoints[Start].z() != lpoints[End].z() && CurdeltazIndex == 0)
		{
		  deltaz[i] = mag(lpoints[End].z()-lpoints[Start].z());
		  CurdeltazIndex = 1;
		}
   	   }
	}

	//Calculating  lower wall y coordinate
	/*ymin= lpoints[0].y();

	for(label pi=1;pi<lpoints.size();pi++)
	{
		if(lpoints[pi].y()<=ymin)
		{
			ymin = lpoints[pi].y();
		}

	}*/

	pMin = gMin(lpoints); //finiding least point in a parallel programming environment

	//Calculating  distance of face centre from lower wall

	const vectorField& faceCenter = patch().Cf();
	forAll(faceCenter,fCi)
	{

		//Distancey[fCi] = faceCenter[fCi].y() - ymin;
	        Distancey[fCi] = faceCenter[fCi].y() - pMin.y();
	}
}

///////////////////////////////////////////////////////////////////
//In this method we are calculating k0 which is important to determine the range of Energy Spectrum
void DSRFG1FvPatchVectorField::CalcK0AndLs()
{

	Info << "Calculating Le_ " <<endl;
	Le_ = scalarField(size(),0); // which stores resultant of Length scales(L11, L22, L33) at each grid point
	IFstream fin(Lfile); // providing input stream to file having lenght scales

	DLList<scalar> Dy; //Dynamic list stores y coordinate (we took dynamic list because we do not the length of data in the file)
	DLList<vector> DFieldL; //Dynamic list stores L11, L22, L33 values as a vector, if have only one length scale data L11 change 'vector' to 'scalar'.
	scalar a;
	vector L(0,0,0);
	label sizeF=0;

	label index =0;
	direction d=0;

// This do loop is coded considering my input file has data y, L11, L22, L33, If you have input data only one length scale L11 (or) L22 (or) L33
// Just change index%4 to index%2

	do{
		if(index%4 == 0)
		{

			if(index!=0)
			{
			//Assigns vector Length scales to Dynamic List
			  DFieldL.append(L);
			}
			fin.read(a); // Assigns a value to variable 'a' read from Length scales file.
			Dy.append(a); // Appending value of y position to Dynamic List
			index = index+1;
			sizeF = sizeF+1;
			d =0;
		}
		else
		{
			       fin.read(a); // Assigns a value to variable 'a' read from Length scales file.
			       L.replace(d,a); //d is a component of a vector, a assigns value to that vector component
			       index = index+1;
			       d =d+1;
		}
	}while(fin.peek()!=EOF);

/////////////////////////////////////////////////////////////
// From here I transfered data from dynamic list to normal array, because few operations I cannot do by dynamic List

	scalarField yy(sizeF);
	vectorField FieldLL(sizeF);

	label yl = 0;
	forAllIter(DLList<scalar>,Dy,iter)
	{
		yy[yl] = *iter;
		yl = yl+1;
	}

	label Ll = 0;
	forAllIter(DLList<vector>,DFieldL,iter1)
	{
		FieldLL[Ll] = *iter1;
		Ll = Ll+1;
	}

        scalarField y(sizeF-1);
        vectorField FieldL(sizeF-1);

        for(label i=0;i<sizeF-1;i++)
        {
                y[i] = yy[i];
                FieldL[i] = FieldLL[i];
        }

////////////////////////////////////////////////////////////////////
// Assigning Length scale (resultant of L11, L22, L33) to each face center of inlet patch

	const vectorField& faceCenter = patch().Cf();
	vector L_;

// Logic:
//------
// If face center 'y' coordinate falls between y value addressed in length scale file, linear interpolation of Length scale is assigned to face center.
	forAll(faceCenter, fi)
	{
		for(label yi=1;yi<y.size();yi++)
		{
		   if( (faceCenter[fi].y() > y[yi-1] && faceCenter[fi].y() < y[yi])
			|| (faceCenter[fi].y() < y[yi-1] && faceCenter[fi].y() > y[yi]) )
		   {
			L_ = (FieldL[yi]+FieldL[yi-1])/2.0;
			scalar L = sqrt(sqr(L_.x())+sqr(L_.y())+sqr(L_.z()));
			Le_[fi] = L;
			break;
		   }
		   else if(faceCenter[fi].y() == y[yi-1])
		   {
			L_ = FieldL[yi-1];
			scalar L = sqrt(sqr(L_.x())+sqr(L_.y())+sqr(L_.z()));
                        Le_[fi] = L;
			break;
		   }
		   else if(faceCenter[fi].y() == y[yi])
		   {
			L_ = FieldL[yi];
			scalar L = sqrt(sqr(L_.x())+sqr(L_.y())+sqr(L_.z()));
                        Le_[fi] = L;
			break;
		   }
		}

	}

    	  Info << "Calculating k0 and Ls ..." << endl;
		Ls_ = gMax(Le_); // finding Maximum value from a scalar Field (used in parallel programming)
		k0_ = 2.0*PI/Ls_;

}


//In this method, The energy spectrum is calculated for Adamjan method
void DSRFG1FvPatchVectorField::EnergySpectrum()
{

	Info << "Calculating Spectrum ..." << endl;
	const polyPatch& Ppatch = patch().patch();

//'&' is the inner product
//nf() class method returns vector normal face
//Cf() face Center (vector)
	deltax = 2.0*(patch().Cf()-Ppatch.faceCellCentres()) & patch().nf();

	scalar factor00,l_cut,factor01;
	scalarField k_cut(size(),0);
	scalar H_max=0.0,k_max = 0.0;
	scalar k_min = Beta*k0_;
	scalar DeltaK = 0.0;

	for(label xi=0;xi<size();xi++)
	{

		if(deltax[xi]>=deltay[xi] && deltax[xi]>=deltaz[xi] )
		{
        		H_max = deltax[xi];
    		}
    		else if(deltay[xi]>=deltax[xi] && deltay[xi]>=deltaz[xi])
		{
        		H_max = deltay[xi];
    		}
    		else if(deltaz[xi]>=deltax[xi] && deltaz[xi]>=deltay[xi])
		{
        		H_max = deltaz[xi];
    		}
		factor01 = max(deltay[xi],deltaz[xi]);
		factor00 = max(factor01,0.3*H_max)+0.1*Distancey[xi];
		l_cut = 2.0*min(mag(factor00),H_max);
		scalar k_cut1 = 2.0*PI/l_cut;
		k_cut[xi] = k_cut1;

	}

	k_max = 1.5* gMax(k_cut);

	k_m = scalarField(SpectrumParts_,0);

	forAll(k_m,ki)
	{
		k_m[ki] = (((k_max-k_min)/(SpectrumParts_-1.0))*ki)+k_min;
		DeltaK = (k_max - k_min)/(SpectrumParts_-1.0);
	}




	scalar ke,v,fcut,feta;
	scalar k_eta = 2.0*PI/L_eta;

// Calculating Energy Spectrum
	E_ = Field<scalarField>(size(),scalarField(SpectrumParts_,0));


	for(label p=0;p<size();p++)
	{
		for(label i=0;i<SpectrumParts_;i++)
		{
			ke = 2.0*PI/Le_[p];
			v = k_m[i]/ke;
			scalar k_cut1 = k_cut[p];
			fcut = exp(-::pow(4.0*max(k_m[i]-0.9*k_cut1,0.0)/k_cut1,3));
           		feta = exp(-::pow(12.0*k_m[i]/k_eta,2));
			E_[p][i] = (pow4(v)/::pow((1+2.4*::pow(v,2)),(17.0/6.0)))*fcut*feta;
		}
		scalar Denominator =0.0;
		for(label m =0;m<SpectrumParts_;m++)
		{
            		Denominator = Denominator + E_[p][m]*DeltaK;
        	}

        	for(label mm = 0;mm<SpectrumParts_;mm++)
		{
            		E_[p][mm] = DeltaK*E_[p][mm]/Denominator;
        	}
	}

}

//Calculating Random coefficients
void DSRFG1FvPatchVectorField::RandomValuesToPQK()
{

	Info << "Calculating Random Values ..." << endl;
	scalar k_i[3]= {0.0,0.0,0.0};

	scalar a;
        scalar** zeta;
	scalar** zi;
  	zeta = new scalar*[NModes_];
  	  zi = new scalar*[NModes_];

	for(label i = 0 ;i < NModes_ ;i++)
	{
    		zeta[i] = new scalar[3];
    		  zi[i] = new scalar[3];
  	}


	scalar kmodul;
	scalar pmodul;
	scalar qmodul;
	scalar ttttt;

	Info<<"Spectrum parts: "<<SpectrumParts_<<"Modes: "<<NModes_<< endl;

	Omega =Field<scalarField>(SpectrumParts_,scalarField(NModes_,0));


Info << "Omega is defined"<<endl;

  for(label i=0;i<3;i++)
  {
    p[i].mn = Field<scalarField>(SpectrumParts_,scalarField(NModes_,0));
    q[i].mn = Field<scalarField>(SpectrumParts_,scalarField(NModes_,0));
    k[i].mn = Field<scalarField>(SpectrumParts_,scalarField(NModes_,0));
  }


Info << "Intialization of Random values are fiished" << endl;
	for(label nn =0;nn<NModes_;nn++)
	{
		 rand_ = Random(nn);
        	zeta[nn][0] = rand_.GaussNormal();

		 rand_ = Random(nn+1);
        	zeta[nn][1] = rand_.GaussNormal();

		rand_ = Random(nn+2);
        	zeta[nn][2] = rand_.GaussNormal();

		rand_ = Random(nn+3);
        	zi[nn][0] = rand_.GaussNormal();

		rand_ = Random(nn+4);
        	zi[nn][1] = rand_.GaussNormal();

		rand_ = Random(nn+5);
        	zi[nn][2] = rand_.GaussNormal();

    	}
Info << "Assigned Random values to zeta and zi" <<endl;

	for(label m=0;m<SpectrumParts_;m++){
                        ttttt = 2.0*PI*k_m[m]*UMean_;
    		for(label n=0;n<NModes_;n++){

			rand_= Random(n+m+6);
        		k_i[0] = rand_.GaussNormal();

			rand_= Random(n+m+7);
        		k_i[1] = rand_.GaussNormal();

			rand_= Random(n+m+8);
        		k_i[2] = rand_.GaussNormal();

        		kmodul = sqrt(sqr(k_i[0])+sqr(k_i[1])+sqr(k_i[2]));


        		k_i[0] = k_i[0]/kmodul;
        		k_i[1] = k_i[1]/kmodul;
        		k_i[2] = k_i[2]/kmodul;
        		/**computation of k vector**/
        		k[0].mn[m][n] = k_i[0]*k_m[m];
        		k[1].mn[m][n] = k_i[1]*k_m[m];
        		k[2].mn[m][n] = k_i[2]*k_m[m];

			rand_ = Random(m+n);
            		 a = rand_.scalar01();

			/**computation of p vector**/
        		p[0].mn[m][n] = zeta[n][1]*k[2].mn[m][n]-zeta[n][2]*k[1].mn[m][n];
        		p[1].mn[m][n] = zeta[n][2]*k[0].mn[m][n]-zeta[n][0]*k[2].mn[m][n];
        		p[2].mn[m][n] = zeta[n][0]*k[1].mn[m][n]-zeta[n][1]*k[0].mn[m][n];

        		pmodul = sqrt(sqr(p[0].mn[m][n])+sqr(p[1].mn[m][n])+sqr(p[2].mn[m][n]));


        		p[0].mn[m][n] = (p[0].mn[m][n]/pmodul)*sqrt(4.0*a/(NModes_*1.0));
        		p[1].mn[m][n] = (p[1].mn[m][n]/pmodul)*sqrt(4.0*a/(NModes_*1.0));
        		p[2].mn[m][n] = (p[2].mn[m][n]/pmodul)*sqrt(4.0*a/(NModes_*1.0));
        			/**computation of q vector**/

        		q[0].mn[m][n] = zi[n][1]*k[2].mn[m][n]-zi[n][2]*k[1].mn[m][n];
        		q[1].mn[m][n] = zi[n][2]*k[0].mn[m][n]-zi[n][0]*k[2].mn[m][n];
        		q[2].mn[m][n] = zi[n][0]*k[1].mn[m][n]-zi[n][1]*k[0].mn[m][n];

        		qmodul = sqrt(sqr(q[0].mn[m][n])+sqr(q[1].mn[m][n])+sqr(q[2].mn[m][n]));

        		q[0].mn[m][n] = (q[0].mn[m][n]/qmodul)*sqrt(4.0*(1-a)/(NModes_*1.0));
        		q[1].mn[m][n] = (q[1].mn[m][n]/qmodul)*sqrt(4.0*(1-a)/(NModes_*1.0));
        		 q[2].mn[m][n] = (q[2].mn[m][n]/qmodul)*sqrt(4.0*(1-a)/(NModes_*1.0));

			 /**computation of Omega**/
	               rand_= Random(n+m+9);
        	       Omega[m][n] = rand_.GaussNormal()*ttttt;

    		}
	}


for (label i = 0; i < NModes_; i++) {
  delete[] zeta[i];
  delete[] zi[i];
}

  delete[] zeta;
  delete[] zi;

}

//Data input from Reynolds stress files
void DSRFG1FvPatchVectorField::InputR()
{
	R_ = symmTensorField(size());
	IFstream fin(Rfile);

//The Reynolds stress input need to be  provided (y R11 R12 R13 R22 R23 R33) fields
	DLList<scalar> Dy;
	DLList<symmTensor> DFieldR;
	scalar a;
	symmTensor R__(0,0,0,0,0,0);
	label sizeR = 0;

	label index =0;
	direction d=0;

	do{
		if(index%7 == 0)
		{

			if(index!=0)
			{
			  DFieldR.append(R__);
			}
			fin.read(a);
			Dy.append(a);
			index = index+1;
			sizeR = sizeR+1;
			d =0;
		}
		else
		{
			       fin.read(a);
			       R__.replace(d,a);
			       index = index+1;

			       d =d+1;
		}
	}while(fin.peek()!=EOF);

	scalarField y(sizeR);
	symmTensorField FieldR(sizeR);


	label yl = 0;
	forAllIter(DLList<scalar>,Dy,iter)
	{
		if(yl < sizeR)
               {
                  y[yl] = *iter;
              }
                yl = yl+1;

	}

	label Rl = 0;
	forAllIter(DLList<symmTensor>,DFieldR,iter1)
	{
		 if(Rl < sizeR)
                {
                  FieldR[Rl] = *iter1;
                }
                Rl = Rl+1;

	}

	const vectorField& faceCenter = patch().Cf();

	forAll(faceCenter, fi)
	{
		for(label yi=1;yi<y.size();yi++)
		{
		   if( (faceCenter[fi].y() > y[yi-1] && faceCenter[fi].y() < y[yi])
			|| (faceCenter[fi].y() < y[yi-1] && faceCenter[fi].y() > y[yi]) )
		   {
			R_[fi] = (FieldR[yi]+FieldR[yi-1])/2.0;
			break;
		   }
		   else if(faceCenter[fi].y() == y[yi-1])
		   {
			R_[fi] = FieldR[yi-1];
			break;
		   }
		   else if(faceCenter[fi].y() == y[yi])
		   {
			R_[fi] = FieldR[yi];
			break;
		   }
		}

	}


}

void DSRFG1FvPatchVectorField::InputU()
{
	U_ = vectorField(size());
	IFstream fin(Ufile);

	DLList<scalar> Dy;
	DLList<vector> DFieldU;
	scalar a;
	vector U(0,0,0);
	label sizeU = 0;

	label index =0;
	direction d=0;

	do{
		if(index%4 == 0)
		{

			if(index!=0)
			{
			  DFieldU.append(U);
			}
			fin.read(a);
			Dy.append(a);
			index = index+1;
			sizeU = sizeU+1;
			d =0;
		}
		else
		{
			       fin.read(a);
			       U.replace(d,a);
			       index = index+1;

			       d =d+1;
		}
	}while(fin.peek()!=EOF);

	scalarField y(sizeU);
	vectorField FieldU(sizeU);

	label yl = 0;
	forAllIter(DLList<scalar>,Dy,iter)
	{
		if(yl < sizeU)
                {
                   y[yl] = *iter;
                }
                yl = yl+1;
	}

	label Ul = 0;
	forAllIter(DLList<vector>,DFieldU,iter1)
	{
		if(Ul < sizeU)
                {
                   FieldU[Ul] = *iter1;
                }
                Ul = Ul+1;
	}


	const vectorField& faceCenter = patch().Cf();

	forAll(faceCenter, fi)
	{
		for(label yi=1;yi<y.size();yi++)
		{
		   if( (faceCenter[fi].y() > y[yi-1] && faceCenter[fi].y() < y[yi])
			|| (faceCenter[fi].y() < y[yi-1] && faceCenter[fi].y() > y[yi]) )
		   {
			U_[fi] = (FieldU[yi]+FieldU[yi-1])/2.0;
			break;
		   }
		   else if(faceCenter[fi].y() == y[yi-1])
		   {
			U_[fi] = FieldU[yi-1];
			break;
		   }
		   else if(faceCenter[fi].y() == y[yi])
		   {
			U_[fi] = FieldU[yi];
			break;
		   }
		}

	}

// Calculating Bulk velocity from the inlet patch
	const vectorField& FaceAreaVector = patch().Sf();
        const scalarField& faceArea = patch().magSf();
        vector UMean_vec = gSumCmptProd(U_,FaceAreaVector)/gSum(faceArea);
        UMean_ = mag(UMean_vec);
        Info << "UMean: " << UMean_ << endl;

}


 void DSRFG1FvPatchVectorField::autoMap
 (
     const fvPatchFieldMapper& m
 )
 {

     fixedValueFvPatchField<vector>::autoMap(m);

 }


 void DSRFG1FvPatchVectorField::rmap
 (
     const fvPatchField<vector>& ptf,
     const labelList& addr
 )
 {

     fixedValueFvPatchField<vector>::rmap(ptf, addr);

 }

void DSRFG1FvPatchVectorField::updateCoeffs()
{
   if (this->updated())
    {
        return;
    }
    scalar f,t;



  if (curTimeIndex_ != this->db().time().timeIndex())
  {

//All input data is calaculated initially
	if(DataCalIndex_ == 0)
  	{

		CalculateDeltas();
		InputU();
		InputR();
		CalcK0AndLs();
		EnergySpectrum();
		RandomValuesToPQK();
		DataCalIndex_ = 1;

	}


	Info << "Calculating U for Inlet patch ..." << endl;


	vectorField& patchField = *this;
	vectorField Signal(size());
  	scalarField rho11(size(),0);
  	scalarField rho12(size(),0);
  	scalarField rho13(size(),0);
  	scalarField rho22(size(),0);
  	scalarField rho23(size(),0);
  	scalarField rho33(size(),0);
	t = (this->db().time().timeOutputValue());
	//scalar timeIndex = this->db().time().timeIndex();

  const vectorField& cellFaceCenter = patch().Cf();

  forAll(patchField,si){
	  const scalar& x = cellFaceCenter[si].x();
  	  const scalar& y = cellFaceCenter[si].y();
  	  const scalar& z = cellFaceCenter[si].z();
	  scalar u0_m=0.0,u1_m=0.0,u2_m=0.0;
	  Le_[si] = Theta1*Le_[si];
          scalar Tau = Theta2*Le_[si]/UMean_;
	for(label m=0;m<SpectrumParts_;m++)
	{

	    scalar u0_n=0.0,u1_n=0.0,u2_n=0.0;

	    for(label n=0;n<NModes_;n++)
	    {
		f = sqrt(E_[si][m]);
		   scalar argument =
((k[0].mn[m][n]*x+k[1].mn[m][n]*y+k[2].mn[m][n]*z)/(k0_*Le_[si]))+(Omega[m][n]*t/Tau);
          u0_n = u0_n +( p[0].mn[m][n]*cos(argument)*f) +
                  (q[0].mn[m][n]*sin(argument)*f);
          u1_n = u1_n + (p[1].mn[m][n]*cos(argument)*f)+
                  (q[1].mn[m][n]*sin(argument)*f);
          u2_n = u2_n + (p[2].mn[m][n]*cos(argument)*f)+
                  (q[2].mn[m][n]*sin(argument)*f);

          rho11[si] = rho11[si] + sqr(p[0].mn[m][n]*f) +
                      sqr(q[0].mn[m][n]*f);
          rho12[si] = rho12[si] + p[0].mn[m][n]*p[1].mn[m][n]*sqr(f) +
                      q[0].mn[m][n]*q[1].mn[m][n]*sqr(f);
          rho13[si] = rho13[si] + p[0].mn[m][n]*p[2].mn[m][n]*sqr(f) +
                      q[0].mn[m][n]*q[2].mn[m][n]*sqr(f);
          rho22[si] = rho22[si] + sqr(p[1].mn[m][n]*f) +
                      sqr(q[1].mn[m][n]*f);
          rho23[si] = rho23[si] + p[1].mn[m][n]*p[2].mn[m][n]*sqr(f) +
                      q[1].mn[m][n]*q[2].mn[m][n]*sqr(f);
          rho33[si] = rho33[si] + sqr(p[2].mn[m][n]*f) +
                      sqr(q[2].mn[m][n]*f);
	    }
		u0_m = u0_m +u0_n;
		u1_m = u1_m +u1_n;
		u2_m = u2_m +u2_n;
	}
		Signal[si] = vector(u0_m,u1_m,u2_m);

}//End for Calculating Signal

Info << "Calculating Primary Field" << endl;

rho11 = 0.5*rho11;
rho12 = 0.5*rho12;
rho13 = 0.5*rho13;
rho22 = 0.5*rho22;
rho23 = 0.5*rho23;
rho33 = 0.5*rho33;

Info << "Anisotropic Tensor" <<endl;

tensorField a_(size(),tensor::zero);

scalarField R11(R_.component(symmTensor::XX));
scalarField R12(R_.component(symmTensor::XY));
scalarField R13(R_.component(symmTensor::XZ));
scalarField R22(R_.component(symmTensor::YY));
scalarField R23(R_.component(symmTensor::YZ));
scalarField R33(R_.component(symmTensor::ZZ));

scalarField a11(sqrt(R11/rho11));
scalarField a22(sqrt((R22-sqr(R12)/R11)/(rho22-sqr(rho12)/rho11)));
scalarField a21(R12/sqrt(R11*rho11)- a22*rho12/rho11);
scalarField denom(sqr(rho13)*rho22-2.0*rho13*rho12*rho23+rho11*sqr(rho23)
            +sqr(rho12)*rho33-rho11*rho22*rho33);

scalarField a33
(
    sqrt(
      (
       (R33*sqr(a11)*sqr(a22)*(sqr(rho12)-rho11*rho22)+sqr(R13)*
       (sqr(a21)*rho11+sqr(a22)*rho22+2.0*a21*a22*rho12))/denom+
       (-2.0*R13*R23*a11*(a22*rho12+a21*rho11)+sqr(R23)*sqr(a11)*rho11)/denom
      )/sqr(a11*a22)
     )
);

scalarField a32( (a33*(rho13*rho12-rho23*rho11)+(R23*rho11/a22)
-(R13*(a21*rho11+a22*rho12))/(a11*a22) )/(rho22*rho11-sqr(rho12)));
scalarField a31(
(a33*(rho23*rho12-rho13*rho22)-(R23*rho12/a22)
+(R13*(a21*rho12+a22*rho22))/(a11*a22) )/(rho22*rho11-sqr(rho12)));


a_.replace(tensor::XX,a11);
a_.replace(tensor::YX,a21);
a_.replace(tensor::YY,a22);
a_.replace(tensor::ZX,a31);
a_.replace(tensor::ZY,a32);
a_.replace(tensor::ZZ,a33);

	patchField = U_ + (a_ & Signal);
	const vectorField& AreaVector = patch().Sf();
	const scalarField& PatchArea = patch().magSf();
	vector Ub_T = gSumCmptProd(patchField,AreaVector)/gSum(PatchArea);
        //scalar Ub = 17.56;
 	Info << "Instantaneous Bulk Velocity: " << Ub_T << endl;
	scalar magUb_T = sqrt(sqr(Ub_T.x())+sqr(Ub_T.y())+sqr(Ub_T.z()));
	patchField = patchField * (UMean_/magUb_T);

	curTimeIndex_ = this->db().time().timeIndex();

  Info <<"End of AssignMent"<<endl;

}//End of time loop

    fixedValueFvPatchField<vector>::updateCoeffs();

}//End of updateCoeff Method


void DSRFG1FvPatchVectorField::write(Ostream& os) const
{

   Info << "Writing to dictionary ..." << endl;
    fvPatchField<vector>::write(os);
    os.writeKeyword("SpectrumParts")
        << SpectrumParts_ << token::END_STATEMENT << nl;
    os.writeKeyword("NModes") << NModes_
	     << token::END_STATEMENT << nl;
    os.writeKeyword("U") << Ufile
	     << token::END_STATEMENT << nl;
    os.writeKeyword("R") << Rfile
	     << token::END_STATEMENT << nl;
    os.writeKeyword("L") << Lfile
             << token::END_STATEMENT << nl;
    os.writeKeyword("Theta1") << Theta1
             << token::END_STATEMENT << nl;
    os.writeKeyword("Theta2") << Theta2
             << token::END_STATEMENT << nl;

 this->writeEntry("value", os);

}

makePatchTypeField
(
  fvPatchVectorField,
  DSRFG1FvPatchVectorField
);
}
