#include <libplugin/plugin.h>
#include <psifiles.h>
#include <psi4-dec.h>
#include <physconst.h>
#include <libparallel/parallel.h>
#include <liboptions/liboptions.h>
#include <libmints/mints.h>
#include <libpsio/psio.h>
#include <libmints/rel_potential.h>
#include <libmints/integral.h>

#include "x2cint.h"

INIT_PLUGIN

namespace psi{ namespace x2c {

void ComputeX2CIntegrals(Options& options);

extern "C"
int read_options(std::string name, Options &options)
{
    if (name == "X2C"|| options.read_globals()) {
        /*- The amount of information printed
            to the output file -*/
        options.add_int("PRINT", 1);

        /*- Whether to uncontract the basis set in a dual basis calculation -*/
        options.add_str("X2C_BASIS","");

        /*- The uncontracted basis set to use in a dual basis calculation -*/
        options.add_str("BASIS_UN","");

        /*- The uncontracted basis set to use in a dual basis calculation -*/
        options.add_str("X2C_BASIS","");

        /*- The uncontracted basis set to use in a dual basis calculation -*/
        options.add_str("X2CMODE","CLASS","CLASS CODE");
    }

    return true;
}

extern "C"
PsiReturnType x2c(Options &options)
{
    if (options.get_str("X2CMODE") == "CLASS"){
        X2CInt x2cint;
        SharedMatrix T,V;
        x2cint.compute(T,V,options);
    }else{
        ComputeX2CIntegrals(options);
    }
    return Success;
}

void ComputeX2CIntegrals(Options& options)
{
    outfile->Printf("\n\n         ---------------------------------------------------------");
    outfile->Printf("\n                                   X2C");
    outfile->Printf("\n                by Prakash Verma and Francesco Evangelista");
    outfile->Printf("\n         ---------------------------------------------------------\n");

    int print = options.get_int("PRINT");

    // Read the molecule information
    boost::shared_ptr<Molecule> molecule = Process::environment.molecule();

    // Form basis object:
    // Create a basis set parser object.
    boost::shared_ptr<BasisSetParser> parser(new Gaussian94BasisSetParser());
    // Construct a new basis set.
    boost::shared_ptr<BasisSet> aoBasis = BasisSet::pyconstruct_orbital(molecule, "BASIS",options.get_str("BASIS"));

//    // If the uncontracted basis is not provided,
//    if (not options["BASIS_UN"].has_changed()){
//        aoBasis = BasisSet::construct(parser, molecule, "BASIS"); //pv
//    }else{
//        aoBasis = BasisSet::construct(parser, molecule, "BASIS_UN"); //pv
//    }

    //lets another BasisSet object
    boost::shared_ptr<BasisSet> aoBasis_cont = BasisSet::pyconstruct_orbital(molecule, "BASIS",options.get_str("BASIS")); //pv

    boost::shared_ptr<IntegralFactory> integral_cont(new IntegralFactory
           (aoBasis_cont, aoBasis, aoBasis, aoBasis)); //pv


    // The integral factory oversees the creation of integral objects
    boost::shared_ptr<IntegralFactory> integral(new IntegralFactory
            (aoBasis, aoBasis, aoBasis, aoBasis));
    
    // Check the point group of the molecule. If it is not set, set it.
    if (!molecule->point_group()) {
        molecule->set_point_group(molecule->find_point_group());
    }

    // Create an SO basis...we need the point group for this part.
    boost::shared_ptr<SOBasisSet> soBasis(new SOBasisSet(aoBasis, integral));


    boost::shared_ptr<SOBasisSet> soBasis_cont(new SOBasisSet(aoBasis_cont, integral_cont));

//    boost::shared_ptr<PetiteList> pet(new PetiteList(basisset_, integral_));
//    AO2SO_ = pet->aotoso();

    // Obtain the dimension object to initialize the factory.
    const Dimension nsopi = soBasis->dimension();

     const Dimension nsopi_cont = soBasis_cont->dimension();

    // Create a Dimension object for the spinors
    Dimension nsspi = nsopi + nsopi;

    // Matrix factory for matrices of dimension nbf x nbf
    boost::shared_ptr<MatrixFactory> soFactory(new MatrixFactory);
    soFactory->init_with(nsopi,nsopi);


    boost::shared_ptr<MatrixFactory> soFactory_cont(new MatrixFactory);
    soFactory_cont->init_with(nsopi_cont,nsopi);

    boost::shared_ptr<MatrixFactory> soFactory_cont_cont(new MatrixFactory);
    soFactory_cont_cont->init_with(nsopi_cont,nsopi_cont);

    // Matrix factory for matrices of dimension 2 nbf x 2 nbf
    boost::shared_ptr<MatrixFactory> ssFactory(new MatrixFactory);
    ssFactory->init_with(nsspi,nsspi);

    // Get the number of irreducible representations
    int nirrep = nsopi.n();

    nsopi.print();
    nsspi.print();

    // N.B. This should be called after the basis has been built, because the geometry has not been
    // fully initialized until this time.
    boost::shared_ptr<OneBodySOInt> sOBI(integral->so_overlap());  
    boost::shared_ptr<OneBodySOInt> tOBI(integral->so_kinetic());
    boost::shared_ptr<OneBodySOInt> vOBI(integral->so_potential());
    boost::shared_ptr<OneBodySOInt> wOBI(integral->so_rel_potential());

    // Overlap integral in the contracted basis
    boost::shared_ptr<OneBodySOInt> sOBI_cont(integral_cont->so_overlap());

    // Form the one-electron integral matrices from the matrix factory
    SharedMatrix sMat(soFactory->create_matrix("Overlap"));
    SharedMatrix tMat(soFactory->create_matrix("Kinetic"));
    SharedMatrix vMat(soFactory->create_matrix("Potential"));
    SharedMatrix wMat(soFactory->create_matrix("Relativistic Potential"));

    // Compute the one electron integrals, telling each object where to store the result
    sOBI->compute(sMat);
    tOBI->compute(tMat);
    vOBI->compute(vMat);
    wOBI->compute(wMat);

    SharedMatrix sMat_cont(soFactory_cont->create_matrix("Overlap"));
    sOBI_cont->compute(sMat_cont);

#if X2CDEBUG
    sMat->print();
    tMat->print();
    vMat->print();
    wMat->print();
    sMat_cont->print();
#endif



    // Form the Dirac Hamiltonian
    //      | V       T       |
    //  d = |                 |
    //      | T  1/4c^2 W - T |

    SharedMatrix dMat(ssFactory->create_matrix("Dirac Hamiltonian"));
    SharedMatrix SXMat(ssFactory->create_matrix("SX Hamiltonian"));

    for (int h = 0; h < nirrep; ++h){
        for (int p = 0; p < nsopi[h]; ++p){
            for (int q = 0; q < nsopi[h]; ++q){
                double Vpq = vMat->get(h,p,q);
                double Tpq = tMat->get(h,p,q);
                double Wpq = wMat->get(h,p,q);
                double Spq = sMat->get(h,p,q);

                /*         | S        0 |
              * SXMat=  |            |
              *         |0   T/2c**2 |
              */

                SXMat->set(h,p,q,Spq);
                SXMat->set(h,p + nsopi[h],q + nsopi[h],0.5 *Tpq / (pc_c_au * pc_c_au));

                // Set the large-large component block
                dMat->set(h,p,q,Vpq);
                // Set the small-large component block
                dMat->set(h,p + nsopi[h],q,Tpq);
                // Set the large-small component block
                dMat->set(h,p,q + nsopi[h],Tpq);
                // Set the small-small component block
                dMat->set(h,p + nsopi[h],q + nsopi[h],0.25 * Wpq/ (pc_c_au * pc_c_au)  - Tpq);
            }
        }
    }

#if X2CDEBUG
    dMat->print();
#endif

/* Prakash
 * C_LS_Mat   Eigenvector that contains both large and small component
 * E_LS_Mat   EigenValues  if Dirac Hamiltonian
 */

     SharedMatrix C_LS_Mat(ssFactory->create_matrix("Dirac EigenVectors"));
     SharedVector E_LS_Mat(new Vector("Dirac EigenValues",nsspi));
     SharedMatrix dtmpMat(ssFactory->create_matrix("Dirac tmp Hamiltonian"));


 /* Prakash
  * Step 3.
  * Diagonalize the Dirac Hamiltonian with X2C 1e Hamiltonian
  */

     SXMat->power(-1.0/2.0);
     dMat->transform(SXMat);
     dMat->diagonalize(dtmpMat,E_LS_Mat);                  //diagonalize Dfock
     C_LS_Mat->gemm(false, false, 1.0, SXMat, dtmpMat, 0.0 );  // C = X C'

#if X2CDEBUG
     C_LS_Mat->print();    
#endif

    E_LS_Mat->print();

#if X2CDEBUG
     double pv;
     outfile->Printf("\n 4c-Dirac Eigenvalues");
     double *pv_e = E_LS_Mat->pointer();
     for (int h = 0; h < nirrep; ++h){
         for (int p = 0; p < 2*nsopi[h]; ++p){
             pv = *(pv_e+p);
             outfile->Printf("\n %5d  %28.18F", p+1, pv);
         }
     }
#endif


/* Prakash
 * Step 4.
 * Divide positive energy part C_LS_Mat into two parts C_Large and C_Small
 * X can be constructed as X * C_Large = C_Small
 * more appropriately as   X = C_Small *  C_Large {^-1}
 */

     SharedMatrix clMat(soFactory->create_matrix("Large EigenVectors"));
     SharedMatrix csMat(soFactory->create_matrix("Small EigenVectors"));
     SharedMatrix  xMat(soFactory->create_matrix("Xmatrix"));


 /*
  * collect the correct Matrix element from Eivenvectos
  */

     for (int h = 0; h < nirrep; ++h){
         for (int p = 0; p < nsopi[h]; ++p){
             for (int q = 0; q < nsopi[h]; ++q){
                 double Lpq = C_LS_Mat->get(h,p,q+nsopi[h]);
                 double Spq = C_LS_Mat->get(h,p+nsopi[h],q+nsopi[h]);
                 // Set the large-large component block
                 clMat->set(h,p,q, Lpq);
                 // Set the small-small component block
                 csMat->set(h,p,q, Spq);
             }
         }
     }

/*
 * Take the inverse of large or upper component
 */

     clMat->general_invert();           // Find C_Large inverse


/*
 * FORM X = C_small * (C_large)^{-1}
 */

     xMat->gemm(false, false, 1.0, csMat, clMat, 0.0 );
#if X2CDEBUG
     xMat->print();
#endif

/* Prakash
 * step 5
 * S_{tilda} = S + 1/2c**2  X^{dagger} T X
 * since we need X^{dagger} T later so we keep that around for little longer
 */


/*
 * FORM X^ T X
 * where X = C_small * ( C_large ) ^{-1}
 * and  T is the kinetic energy
 */


      clMat->transform(xMat,tMat,xMat);



/*
 * FORM  S_{tilda} = S + (X^ T X) / 2c**2
 * S is the overlap matrix
 */

     clMat->scale(1.0/(2.0*pc_c_au * pc_c_au));             // scale by 1/2c**2
     csMat->copy(sMat);                                 // S_tilda = S + X^ T X
     csMat->add(clMat);

#if X2CDEBUG
     csMat->print();
#endif

 /*Prakash
  * Step 6
  * R = S^{-1/2} (S ^{-1/2} S_tilda S^{-1/2})^{-1/2} S^{1/2}
  */

     SharedMatrix  Evec(soFactory->create_matrix("Eigenvector S matrix"));
     SharedVector  Eval(new Vector("Eigenvalues S Matrix", nsopi));
     SharedMatrix  sTmp(soFactory->create_matrix("S tmp matrix"));
     SharedMatrix  S_inv(soFactory->create_matrix("Eigenvector S matrix"));

     SharedMatrix  sTmp_cont(soFactory_cont->create_matrix("S_cont tmp matrix"));
     SharedMatrix  sTmp1_cont(soFactory_cont_cont->create_matrix("S_cont tmp1 matrix"));

/*
 * FORM  S^{-1/2}
 * Matrix = U diag(Matrix) U^
 */

      S_inv->copy(sMat);
      S_inv->power(-1.0/2.0);
#if X2CDEBUG
      S_inv->print();
#endif

/*
 * FORM  (S ^{-1/2} S_tilda S^{-1/2})^{-1/2}
 */

      clMat->gemm(false, false, 1.0, S_inv, csMat, 0.0 );  // S^{-1/2} S_tilda
       sTmp->gemm(false, false, 1.0, clMat, S_inv,0.0);      // S^{-1/2} S_tilda S^{-1/2}
       sTmp->power(-1.0/2.0);



/*
 * S^{-1/2} (S ^{-1/2} S_tilda S^{-1/2})^{-1/2} S^{1/2}
 */

       clMat->gemm(false, false,  1.0, S_inv, sTmp, 0.0);    // S^{-1/2} * (sTmp)
       S_inv->general_invert();                              // S^{1/2}
       //lets use csMat instead of rMat
       csMat->gemm(false, false, 1.0, clMat, S_inv, 0.0 );   // R=S^{-1/2} (S_inv S_tilda S_inv)^{-1/2}S^{1/2}

#if X2CDEBUG
       csMat->print();
#endif

/*
 * FORM XR matrix
 */

      // instead of xrMat I am going to use S_inv

       S_inv->gemm(false, false, 1.0, xMat, csMat, 0.0 );     // XR = X R matrix

#if X2CDEBUG
       S_inv->print();
#endif


/*Prakash
 * step 7
 * construct h^{FW}_{+}
 *                   =   R^ T  XR +
 *                       R^ V  R  +
 *                    (XR)^ T  R -
 *                    (XR)^ T (XR) +
 *                    (XR)^ W' XR
 *
 * where W' is the scaled version of W i.e W/4c**2
 */


      SharedMatrix  HMat(soFactory->create_matrix("Hfw matrix"));
      SharedMatrix  TfwMat(soFactory->create_matrix("Tfw matrix"));
      SharedMatrix  VfwMat(soFactory->create_matrix("Vfw matrix"));


      SharedMatrix  SfwMat_cont(soFactory_cont_cont->create_matrix(PSIF_SO_S));
      SharedMatrix  TfwMat_cont(soFactory_cont_cont->create_matrix(PSIF_SO_T));
      SharedMatrix  VfwMat_cont(soFactory_cont_cont->create_matrix(PSIF_SO_V));


/*
 *   Save the projected overlap matrix
*/
      clMat->copy(sMat);
      clMat->general_invert();
      sTmp_cont->gemm(false,false,1.0,sMat_cont,clMat,0.0);
      SfwMat_cont->gemm(false,true,1.0,sTmp_cont,sMat_cont,0.0);

      // Save the modified kinetic integrals to disk
      SfwMat_cont->save(_default_psio_lib_,PSIF_OEI);

/*
 *   R^ T XR
 */


       clMat->transform(csMat,tMat,S_inv);
       HMat->copy(clMat);
       TfwMat->copy(clMat);


#if X2CDEBUG
        HMat->print();
#endif

 /*
  *  (XR)^ T R  - (XR)^ T (XR)
  */


        sTmp->gemm(true, false, 1.0, S_inv, tMat, 0.0);
        clMat->gemm(false, false, 1.0, sTmp, csMat,0.0);
        HMat->add(clMat);
        TfwMat->add(clMat);
        clMat->gemm(false,false,1.0,sTmp,S_inv,0.0);
        HMat->subtract(clMat);
        TfwMat->subtract(clMat);
        /*
         * S^{-1}HS^{-1} same basis set
         */

        clMat->copy(sMat);
        clMat->general_invert();
        TfwMat->transform(clMat);

        /*
          * S_p Htran S_p
          */


        sTmp_cont->gemm(false,false,1.0,sMat_cont,TfwMat,0.0);
        TfwMat_cont->gemm(false,true,1.0,sTmp_cont,sMat_cont,0.0);

        // Save the modified kinetic integrals to disk
        TfwMat_cont->save(_default_psio_lib_,PSIF_OEI);

#if X2CDEBUG
        HMat->print();
        TfwMat->print();
#endif


/*
 *   R^ V R
 */



       vMat->transform(csMat);
       HMat->add(vMat);
       VfwMat->add(vMat);

#if X2CDEBUG
        HMat->print();
#endif

/*
 *   (XR)^ W' XR
 */

       wMat->transform(S_inv);
       wMat->scale(1.0/(4.0*pc_c_au * pc_c_au));
       HMat->add(wMat);  // NOTE: does HMat get transformed?
       VfwMat->add(wMat);


       /*
        * S^{-1}HS^{-1} same basis set
        */

      clMat->copy(sMat);
      clMat->general_invert();

      VfwMat->transform(clMat);

        /*
         * S_p Htran S_p
         */


        sTmp_cont->gemm(false,false,1.0,sMat_cont,VfwMat,0.0);
        VfwMat_cont->gemm(false,true,1.0,sTmp_cont,sMat_cont,0.0);

       // Save the modified potential integrals to disk
       VfwMat_cont->save(_default_psio_lib_,PSIF_OEI);

#if X2CDEBUG
       HMat->print();
#endif

/*
 * Diagonalize the Hamiltonian
 */

//       sMat->power(-0.5);
//       HMat->transform(sMat);
//       HMat->diagonalize(Evec,Eval);

//       outfile->Printf("\n  Eigenvalues of the positive-energy states from H^FQ_++");
//       Eval->print();


       SharedMatrix  Evec_c(soFactory_cont_cont->create_matrix("Eigenvector S1 matrix"));
       SharedVector  Eval_c(new Vector("Eigenvalues S Matrix", nsopi_cont));
       SharedMatrix  sTmp_c(soFactory_cont_cont->create_matrix("S1 tmp matrix"));
       SharedMatrix  S_inv_c(soFactory_cont_cont->create_matrix("Eigenvector S1 matrix"));
       SharedMatrix clMat_c(soFactory_cont_cont->create_matrix("Large1 EigenVectors"));
       SharedMatrix HMat_c(soFactory_cont_cont->create_matrix("Large2 EigenVectors"));

       boost::shared_ptr<IntegralFactory> integral_cont_cont(new IntegralFactory
                                                             (aoBasis_cont, aoBasis_cont, aoBasis, aoBasis));
       boost::shared_ptr<OneBodySOInt> sOBI_cont_c(integral_cont_cont->so_overlap());
       SharedMatrix sMat_c(soFactory_cont_cont->create_matrix("Overlap"));
       sOBI_cont_c->compute(sMat_c);

       HMat_c->zero();
       HMat_c->add(VfwMat_cont);
       HMat_c->add(TfwMat_cont);
       sMat_c->power(-0.5);
       HMat_c->transform(sMat_c);
       HMat_c->diagonalize(Evec_c,Eval_c);
       outfile->Printf("\n  Eigenvalues of the positive-energy states from H^FQ_++ in the contracted basis");
       Eval_c->print();
}

}} // End Namespaces

