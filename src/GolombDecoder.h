//
//  GolombDecoder.h
//  Predictor RLS
//
//  Created by Ignacio Capurro on 2/7/14.
//  Copyright (c) 2014 Ignacio Capurro. All rights reserved.
//

#ifndef Predictor_RLS_GolombDecoder_h
#define Predictor_RLS_GolombDecoder_h

#include <iostream>
#include "BinaryFileReader.h"
#include "Golomb.h"
#include <cmath>

using namespace std;

class GolombDecoder : public Golomb {
public:
    GolombDecoder(BinaryFileReader * );
    int decodeSample();
private:
    BinaryFileReader * file;
    int binaryDecode(unsigned);
};

#endif
