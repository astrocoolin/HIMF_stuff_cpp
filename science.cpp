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

        const std::string fileName("empty.fits");
        // Create a new FITS object, specifying the data type and axes for the primary
        // image. Simultaneously create the corresponding file.

        // this image is unsigned short data, demonstrating the cfitsio extension
        // to the FITS standard.

        pFits.reset( new FITS(fileName , FLOAT_IMG , naxis , naxes ) );
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

    //for (long j = 0; j < vectorLength; ++j) row[j] = j;
    //for (long j = 0; j < numberOfRows; ++j) col[j] = j;

    std::valarray<float> array(nelements);
    array = 0.0;
    //array[std::slice(vectorLength*numberOfRows*static_cast<int>(k),vectorLength*numberOfRows,1)] = k;
    long  fpixel(1);
        pFits->pHDU().write(fpixel,nelements,array);


    // PHDU's friend ostream operator. Doesn't print the entire array, just the
    // required & user keywords, and is provided largely for testing purposes [see
    // readImage() for an example of how to output the image array to a stream].

    // std::cout << pFits->pHDU() << std::endl;

    return 0;
