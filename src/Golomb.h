//
//  Golomb.h
//  Predictor RLS
//
//  Created by Ignacio Capurro on 2/7/14.
//  Copyright (c) 2014 Ignacio Capurro. All rights reserved.
//

#ifndef Predictor_RLS_Golomb_h
#define Predictor_RLS_Golomb_h

class Golomb {
public:
    Golomb(){
        accumulatedError = 0;
 //       negativeSamples2 = 0;
        samples = 0;
        k = 1;
    };
protected:
    unsigned accumulatedError;
//	unsigned negativeSamples2;
	unsigned samples;
	unsigned k;
    
//    unsigned rice(int n) { return (n < 0) ? (unsigned(-n) << 1) - 1 : (unsigned)n << 1;};
//	int unRice(unsigned n) { return (n & 1)? -((int)((n+1) >> 1)) : n >> 1;};
    
};

#endif
