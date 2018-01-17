## TMPred C++


As presented at the TMPred web-service (http://www.ch.embnet.org/software/TMPRED_form.html): "The TMPred (Transmembrane
Predictor) makes a prediction of membrane-spanning regions and their orientation. The algorithm is based on the statistical 
analysis of TMbase, a database of naturally occurring transmembrane proteins. The prediction is made using a combination of 
several weight-matrices for scoring''. This project offers a C++11 updated off-line version of TMPred for fasta files. Also, 
our version allows to modified all the internal parameters use by TMPred to detect transmembranes. 
For more information of how TMPred works please refer to http://www.ch.embnet.org/software/tmbase/TMBASE_doc.html


## Requirements

In order to successfully compile TMPred in Ubuntu-Linux machine you need: 
- An updated g++ version compiler
- An installed updated version of the cmake build system (version ≥ 3.5)
- An installed updated version of the boost library

##Compile

To be able to compile the codes: 
- Check that the boost library is correctly pointed in the CMakeLists.txt file
- Create the build folder if it does not exists
	- mkdir build
- Access the build folder and tun
	- cmake ..
	- make

## Methods

-[TMPred] :
	Use: ./TMPred file_name <opt>
		<file_name>: Name of the fasta file to be analyzed
	Options:
  		-h [ --help ]           print usage message
  		-o [ --output ] arg     pathname for output. Default: file_name.tmpred
  		-p [ --print-mode ] arg data to be printed: 
			0 - empty output
			1 - all TMPred output
			2 - only predicted models information
			Default: 1
  		-m [ --min-len ] arg    minimal length of transmembrane sequence. Default: 17
  		-M [ --max-len ] arg    maximal length of transmembrane sequence. Default: 35
  		-l [ --low-osl ] arg    low orientational significance level. Default: 80
  		-i [ --high-osl ] arg   high orientational significance level. Default: 200
  		-s [ --tm-osl ] arg     TM-existence significance level. Default: 500
  		-a [ --avg-osl ] arg    average orientation significance level. Default: 80
  		-g [ --gen-graph ]      if set, generates the score graphs of each sequence into output.tmpred_graph. Default: not set

	
	Example: ./TMPred example -o example -p 1 -g -m 15 -M 33 -l 79 -s 501 -a 81 -i 220
	output:	example.tmpred and example.tmpred_graphs
        
This projects include an example file in the example folder. The -g options generates the scores graphs using the gnu-iostream.h code which 
is the C++ interface to gnuplot implemented by Daniel Stahlke (http://www.stahlke.org/dan/gnuplot-iostream).

			
Note: While our tool was implemented in a Linux machine, it had been successfully tested using Cygwin (https://www.cygwin.com/) in a Window machine 
(there could be problems with the gnuplot-iostream.h code, which we recommned to comment from the tmpred.h code, including any call to its methods).

