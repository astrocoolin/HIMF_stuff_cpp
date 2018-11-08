#include <iostream>
#include <fstream>
#include <thread>
#include <CCfits/CCfits>
#include <cmath>
#include <iomanip>
using namespace CCfits;
using namespace std;

int deffile_print(double *radi,double *vrot,double *sbr,double *z,int radi_len, double inclination) {
  ofstream deffile;
  deffile.open("cube_input.def",ios::out);
  deffile << "LOGNAME=Logfile.log\n";
  deffile << "ACTION=1\n";
  deffile << "PROMPT=1\n";
  unsigned NCORES = std::thread::hardware_concurrency();
  deffile << "NCORES=" << NCORES;
  deffile << "\n";

  deffile <<"\nINSET=empty.fits";//+insetd;
  deffile <<"\nOUTSET=Cube_base.fits";//+outset;
  deffile <<"\nOUTCUBUP=10000000";
  deffile <<"\n";
  deffile <<"\nPROGRESSLOG=";
  deffile <<"\nTEXTLOG=";
  deffile <<"\n";
  deffile <<"\nRMS= 1.00000e-10";
  deffile <<"\nBMAJ=0";
  deffile <<"\nBMIN=0";
  deffile <<"\nBPA=0";
  deffile <<"\n";
  deffile <<"\nNDISKS=";
  deffile <<"\nNUR=" << radi_len+1;
  deffile <<"\n\nRADI=";
  for (int i=0;i<=radi_len;i++){
  	deffile <<" +"  << std::scientific<< std::setprecision(3)<<radi[i];
  }
  deffile <<"\n\nVROT=";
  for (int i=0;i<=radi_len;i++){
  	deffile <<" +"  << std::scientific<< std::setprecision(3)<<vrot[i];
  }
  deffile <<"\n\nSBR=";
  for (int i=0;i<=radi_len;i++){
  	deffile <<" +"  << std::scientific<< std::setprecision(3)<<sbr[i];
  }
  deffile <<"\n\nZ=";
  for (int i=1;i<=radi_len;i++){
	if (i==0) {
  	deffile <<" +"  << std::scientific<< std::setprecision(3)<<z[1];
	}
  	deffile <<" +"  << std::scientific<< std::setprecision(3)<<z[i];
  }

  deffile <<"\nINCL= +" << inclination;
  deffile <<"\nPA= +4.50000E+01";
  deffile <<"\nXPOS= +2.77675430E+02";
  deffile <<"\nYPOS= +7.34348280E+01";
  deffile <<"\nVSYS= +1403.93164636";
  deffile <<"\nSDIS=2.0";
  deffile <<"\nCONDISP=0";
  deffile <<"\nLTYPE= 3";
  deffile <<"\n";
  deffile <<"\nCFLUX=1.000e-06";
  deffile <<"\nPENALTY= 0";
  deffile <<"\nWEIGHT= 0 ";
  deffile <<"\nRADSEP= 0.1 ";
  deffile <<"\nINIMODE= 0 ";
  deffile <<"\nISEED= 8981 ";
  deffile <<"\n";
  deffile <<"\nFITMODE= 1";
  deffile <<"\nLOOPS=0";
  deffile <<"\nMAXITER=";
  deffile <<"\nCALLITE=";
  deffile <<"\nSIZE=";
  deffile <<"\nINTY=";
  deffile <<"\nINDINTY=";
  deffile <<"\nPSSE=";
  deffile <<"\nPSNP=";
  deffile <<"\nPSCO=";
  deffile <<"\nPSSO=";
  deffile <<"\nPSMV=";
  deffile <<"\nPSNF=";
  deffile <<"\nPSII=";
  deffile <<"\nPSFI=";
  deffile <<"\nPSID=";
  deffile <<"\nPSDD=";
  
  if (radi_len > 6) {
        deffile <<"\n\nVARY= INCL 1:5, !INCL 6:8,  PA 1:5 , !PA 6:, !VROT 2:8, !SBR 2:8,  Z0 1:8";
  } else {
        deffile <<"\n\nVARY= INCL 1:5, !INCL 6:" <<radi_len << ",  PA 1:5 , !PA 6:" <<radi_len <<", !VROT 2:" <<radi_len <<", !SBR 2:" <<radi_len <<",  Z0 1:" <<radi_len <<"";
  } 

  deffile <<"\nVARINDX=";
  deffile <<"\nPARMAX=   90  90  360  360  500    1  30";
  deffile <<"\nPARMIN=    0   0 -360 -360   20    0   0";
  deffile <<"\nMODERATE=  3   3    3    3    3    3   3";
  deffile <<"\nDELSTART=0.5   2  0.5    2   16 1E-5   1";
  deffile <<"\nDELEND=  0.1 0.5  0.1  0.5    8 1E-6 0.1";
  deffile <<"\nITESTART= 70  70   70   70   70   70  70";
  deffile <<"\nITEEND=   70  70   70   70   70   70  70";
  deffile <<"\nSATDELT=   2   4    2    4   10 5E-6   1";
  deffile <<"\nMINDELTA=0.1 0.5  0.1  0.5    8 5E-7 0.1";
  deffile <<"\n";



  deffile <<"\nREGPARA=";
  deffile <<"\nTIRSMO=";
  deffile <<"\nCOOLGAL=";
  deffile <<"\nCOOLBEAM=";
  deffile <<"\nTILT=";
  deffile <<"\nBIGTILT=\n";
  
  deffile.close();
  return 0;

}

int emptyFits(){
    // Create a FITS primary array containing a 2-D image
    // declare axis arrays.
    long naxis    =   3;
    long naxes[3] = { 400, 400 ,200};

        // declare auto-pointer to FITS at function scope. Ensures no resources
    // leaked if something fails in dynamic allocation.
    std::auto_ptr<FITS> pFits(0);

    try
    {
        // overwrite existing file if the file already exists.

        const std::string fileName("!empty.fits");

        // Create a new FITS object, specifying the data type and axes for the primary
        // image. Simultaneously create the corresponding file.

        // this image is unsigned short data, demonstrating the cfitsio extension
        // to the FITS standard.

        pFits.reset( new FITS(fileName , DOUBLE_IMG , naxis , naxes ) );
    }
    catch (FITS::CantCreate)
    {
          // ... or not, as the case may be.
          return -1;
    }

    // references for clarity.
    long& vectorLength = naxes[0];
    long& numberOfRows = naxes[1];
    long& channels     = naxes[2];
    long nelements(1);

    // Find the total size of the array.
    // this is a little fancier than necessary ( It's only
    // calculating naxes[0]*naxes[1]) but it demonstrates  use of the
    // C++ standard library accumulate algorithm.
    nelements = std::accumulate(&naxes[0],&naxes[naxis],1,std::multiplies<long>());

    // create a new image extension with a 300x300 array containing float data.
    std::vector<long> extAx(2,300);
    //cout << extAx << "HUH" << endl;
    string newName ("NEW-EXTENSION");
    ExtHDU* imageExt = pFits->addImage(newName,FLOAT_IMG,extAx);

    // create a dummy row with a ramp. Create an array and copy the row to
    // row-sized slices. [also demonstrates the use of valarray slices].
    // also demonstrate implicit type conversion when writing to the image:
    // input array will be of type float.

    std::valarray<float> row(vectorLength);
    std::valarray<float> col(numberOfRows);

    std::valarray<float> array(nelements);
    array = 0.0;
    long  fpixel(1);

    // write the image extension data: also demonstrates switching between
    // HDUs.

    pFits->pHDU().addKey("CDELT1", -0.00111111,"Degrees per Pixel");
    pFits->pHDU().addKey("CRPIX1", 200.00,"Reference Pixel");
    pFits->pHDU().addKey("CRVAL1", 277.675431576,"Reference Pixel Value");
    pFits->pHDU().addKey("CTYPE1", "RA---TAN" ,"Axis Name");
    pFits->pHDU().addKey("CUNIT1", "DEGREE" ,"Axis Units");

    pFits->pHDU().addKey("CDELT2",  0.00111111,"Degrees per Pixel");
    pFits->pHDU().addKey("CRPIX2", 200.00,"Reference Pixel");
    pFits->pHDU().addKey("CRVAL2", 73.4348279512,"Reference Pixel Value");
    pFits->pHDU().addKey("CTYPE2", "DEC---TAN","Axis Name");
    pFits->pHDU().addKey("CUNIT2", "DEGREE" ,"Axis Units");
    
    pFits->pHDU().addKey("CDELT3", 5000.00,"M/S per voxel");
    pFits->pHDU().addKey("CRPIX3", 100.00,"Reference voxel");
    pFits->pHDU().addKey("CRVAL3", 1403931.64636,"Reference Pixel Value");
    pFits->pHDU().addKey("CTYPE3", "VELO-LSR","Axis Name");
    pFits->pHDU().addKey("CUNIT3", "M/S" ,"Axis Units");


    pFits->pHDU().addKey("EPOCH", 2.0000E3,"EPOCH");
    pFits->pHDU().addKey("MAPTYP","MAP","");

    pFits->pHDU().addKey("BUNIT","Jy/Beam","");
    pFits->pHDU().addKey("BMAJ","0.00","");
    pFits->pHDU().addKey("BMIN","0.00","");
    pFits->pHDU().addKey("BPA","0.00","");


    // The function PHDU& FITS::pHDU() returns a reference to the object representing
    // the primary HDU; PHDU::write( <args> ) is then used to write the data.

    pFits->pHDU().write(fpixel,nelements,array);
    return 0;
}


