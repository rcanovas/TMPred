/* TMPred - Transmembrane Region Predictor
 * Copyright (C) 2018 Rodrigo Canovas
 *
 * This program is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with this program.  If not, see http://www.gnu.org/licenses/ .
 * */


#ifndef FDGPI_DEPENDENT_PARTS_H
#define FDGPI_DEPENDENT_PARTS_H

#include <stdlib.h>
#include <iostream>
#include <fstream>
#include <sstream>

namespace fdgpi {

    //tmpred constants
    const int MAXSEQLEN = 6000;
    const int MAXHELIX = 150;

    //global variable
    std::string depend_err = "";

    // There are 9 valid format types:
    const uint8_t IG = 1;
    const uint8_t GenBank = 2;
    const uint8_t NBRF = 3;
    const uint8_t EMBL = 4;
    const uint8_t GCG = 5;
    const uint8_t Strider = 6;
    const uint8_t Fitch = 7;
    const uint8_t Pearson = 8;
    const uint8_t Zuker = 9;
    const uint8_t Unknown = 10;

    // Some constant parameters
    const int8_t NNB = 18;
    const std::string aminos = "ABCDEFGHIKLMNPQRSTVWXYZ*";
    const std::string iubbase = "ACGTUMRWSYKVHDBXN.";
    const std::string igbase = "ACGTUJRLMYKNNNNNN?";
    const std::string primenuc = "ACGTU";
    const std::string protonly = "EFIPQZ";
    const double pi = 3.1415926;


    //Digits ['0'..'9'] are NOT allowed here.
    inline bool isSeqChar(char c) {
        if (c <= ' ')
            return false;
        else if (((c >= 'A') and (c <= 'Z'))
                 or ((c >= 'a') and (c <= 'z'))
                 or (c = '_') or (c = '?') or (c = '*') or (c = '.') or (c = '-')) {
            return true;
        }
        else
            return false;
    }

    inline void open_w(std::string name, std::ofstream &f) {
        std::ofstream out_tmp(name);
        f.swap(out_tmp);
    }

}


#endif //FDGPI_DEPENDENT_PARTS_H
