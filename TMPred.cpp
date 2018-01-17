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

#include <unistd.h>
#include <iostream>
#include "./include/tmpred.h"
#include <boost/program_options.hpp>

int main(int argc, char* argv[]) {

    std::string file, out_file;
    int pr = 1, m = 17, M = 35;
    int t1 = 80, t2 = 200, t3 = 500, t4 = 80;
    bool g = false;
    namespace po = boost::program_options;
    po::options_description desc("Options");
    desc.add_options()
            // First parameter describes option name/short name
            // The second is parameter to option
            // The third is description
            ("help,h", "print usage message")
            ("output,o", po::value<std::string>(), "pathname for output. Default: file_name.tmpred")
            ("print-mode,p", po::value<int>(), "data to be printed: \n   0 - empty output\n   1 - all TMPred output\n   2 - only predicted models information\nDefault: 1")
            ("min-len,m", po::value<int>(), "minimal length of transmembrane sequence. Default: 17")
            ("max-len,M", po::value<int>(), "maximal length of transmembrane sequence. Default: 35")
            ("low-osl,l", po::value<int>(), "low orientational significance level. Default: 80")
            ("high-osl,i", po::value<int>(), "high orientational significance level. Default: 200")
            ("tm-osl,s", po::value<int>(), "TM-existence significance level. Default: 500")
            ("avg-osl,a", po::value<int>(), "average orientation significance level. Default: 80")
            ("gen-graph,g","if set, generates the score graphs of each sequence into output.tmpred_graph. Default: not set");
    po::variables_map vm;
    po::store(po::parse_command_line(argc, argv, desc), vm);

    if (vm.count("help") or argc < 2) {
        std::cout << "Usage: " << argv[0] << " file_name <opt>" << std::endl;
        std::cout << desc << std::endl;
        return 0;
    }
    file = argv[1];
    out_file = file;
    if (vm.count("gen-graph") ) { g = true;}
    if (vm.count("output")) { out_file = vm["output"].as<std::string>();}
    if (vm.count("print-mode")) { pr = vm["print-mode"].as<int>();}
    if (vm.count("min-len")) { m = vm["min-len"].as<int>();}
    if (vm.count("max-len")) { M = vm["max-len"].as<int>();}
    if (vm.count("low-osl")) { t1 = vm["low-osl"].as<int>();}
    if (vm.count("high-osl")) { t2 = vm["high-osl"].as<int>();}
    if (vm.count("tm-osl")) { t3 = vm["tm-osl"].as<int>();}
    if (vm.count("avg-osl")) { t4 = vm["avg-osl"].as<int>();}

    /*std::cout << "parameters: " << std::endl;
    std::cout << "out_file: " << out_file << std::endl;
    std::cout << "pr " << pr << std::endl;
    std::cout << "m " << m << std::endl;
    std::cout << "M " << M << std::endl;
    std::cout << "t1 " << t1 << std::endl;
    std::cout << "t2 " << t2 << std::endl;
    std::cout << "t3 " << t3 << std::endl;
    std::cout << "t4 " << t4 << std::endl;*/
    //std::cout << "g " << g << std::endl;

    fdgpi::tmpred tp(file, out_file, pr, g, m, M, t1, t2, t3, t4);
    tp.process_all();
    return 0;
}
