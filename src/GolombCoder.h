//
//  GolombCoder.h
//  compresoreeg
//
//  Created by Ignacio Capurro on 12/13/13.
//  Copyright (c) 2013 Ignacio Capurro. All rights reserved.
//

#ifndef __compresoreeg__GolombCoder__
#define __compresoreeg__GolombCoder__

//#include <iostream>
//#include "BinaryFileWriter.h"
#include "Golomb.h"
#include <cmath>

using namespace std;

class GolombCoder : public Golomb {
public:
//    GolombCoder(BinaryFileWriter *);
 GolombCoder():Golomb(),bitcount() {}
    void codeSample(unsigned);
    long bitcount;
private:
//    BinaryFileWriter * file;
    void binaryEncode(unsigned, unsigned);
};

#endif /* defined(__compresoreeg__GolombCoder__) */
