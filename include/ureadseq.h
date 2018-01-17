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


#ifndef FDGPI_UREADSEQ_H
#define FDGPI_UREADSEQ_H

#include <vector>
#include "./dependent_parts.h"

namespace fdgpi {


    class rwSeq {

    private:
        int seqlen; //sequence size
        int max_seq; //Maximum sequence size
        std::string seqId; //seq information
        std::string seq; //Sequence or list storage

        //auxiliary variables
        std::string line;
        int err; //0 if all goes ok, otherwise the error number

    public:

        rwSeq() {

        }

        /*
         * input:
         * mseq: Maximum sequence size
         * tseq: Sequence
         * header: Sequence header
         * first: true if we are reading the first sequence
         */
        int readNextSeq(std::istream &in, int &mseq, std::string &tseq, std::string &header, bool &first) {
            int len_line;
            max_seq = mseq;
            err = 0;
            if (first) {
                if (!std::getline(in, line))
                    return 0;
                first = false;
            }
            if (in.eof() or line.length() == 0) {
                //std::cout << "Finished reading the file" << std::endl;
                return 0;
            }
            header = line; //it was read before
            seqlen = 0;
            seq = "";
            while (std::getline(in, line)) {
                if (line.find(">") == 0)
                    break;
                len_line = line.length();
                while (len_line > 0 and line[len_line - 1] == ' ')
                    len_line--;
                if (len_line == 0)
                   break;   //we assume that the file can have empty lines at the end
                addseq(line);
            }
            mseq = seqlen;
            tseq = seq;
            return 1;
        }

    private:

        void addseq(std::string s) {
            for (int i = 0; i < s.length(); i++) {
                if (isSeqChar(s[i])) {  //note that this erase any not valid character from the sequence
                    if (seqlen >= max_seq)
                        err = -3;  //not enough space error
                    else {
                        seq.push_back(s[i]);
                        seqlen++;
                    }
                }
            }
        }

    }; // end CLASS

}


#endif //FDGPI_UREADSEQ_H
