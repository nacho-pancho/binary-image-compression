//
//  GolombDecoder.cpp
//  Predictor RLS
//
//  Created by Ignacio Capurro on 2/7/14.
//  Copyright (c) 2014 Ignacio Capurro. All rights reserved.
//

#include "GolombDecoder.h"

GolombDecoder::GolombDecoder(BinaryFileReader * file) {
    this->file = file;
}

int GolombDecoder::binaryDecode(unsigned k){
	
	unsigned binary = file->readBits(k);
	unsigned unary = file->countZeros();
	file->readBits(1);
    
	unsigned sample = (unary << k) | binary;
	return unRice(sample);
}

int GolombDecoder::decodeSample(){
    int sample = 0;
    if (negativeSamples2 > samples) {
        sample = binaryDecode(k);
        sample = -sample-1;
    } else
        sample = binaryDecode(k);
    
    samples++;
    if (sample < 0)
	{
        negativeSamples2 += 2;
		accumulatedError -= sample;
	}
	else
		accumulatedError += sample;
    for(k=0 ; (samples<<k) < accumulatedError; k++ );
	return sample;
}