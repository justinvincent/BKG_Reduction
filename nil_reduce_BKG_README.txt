The tar file nil_reduce_BKG.tar contains all necessary source code and make files needed to generate
the executables, EKG_phase and reduce_BKG. The code compiles and runs on Sun Solaris running either
on Sun (SPARC) hardware or Solaris 10 installed under VMware running on another operating
system. In principal, this code can be compiled on other Unix (linux) operating systems.
However, several required FORTRAN modules require system subroutines malloc() and free(), which currently are not
supported by gcc. The essential algorithmic outlines are reasonably well documented. Hence, if Solaris
Unix is not available an alternative strategy would be to substitute external signal conditioning
modules.

Questions and comments may be addressed to
A. Z. Snyder	avi@npg.wustl.edu
J. L. Vincent	vincent@nmr.mgh.harvard.edu

EKG_phase
Usage:  EKG_phase <EKG_datfile>
e.g.,   EKG_phase 04_1004_ascii_EKG[.dat]
        option
        -e      dump (expanded) EKG record to enable manual selection of exemplar beat
        -q      disable automatic beat deletion/insertion
        -r      recompile IBIhist and bps after each beat deletion/insertion
        -t<flt> enter manually identified exemplar beat latency in sec
        -k<flt> specify template start to lagged product peak interval in sec (default=0.20)
        -l<flt> specify lagged product interval in msec (default=14.00)
        -Y<int> suppress use of specified component of y to compute z[0]
        -d<flt> manually delete beat within 100 msec of specified time
        -i<flt> manually insert beat at specified time in sec
        -z0<flt>        specify beat detector z[0] criterion in s.d. units (default=2.50)
        -z1<flt>        specify beat detector z[1] criterion in s.d. units (default=2.50)
        -Z      report discrepancy between z[0] vs. z[1] meeting criterion
        -b<flt> specify false beat detector bps tolerance as fraction of local bps (default=0.22)
        -h<flt> specify false beat detector IBI histogram tolerance in s.d. units (default=3.50)
        -j<flt> specify beat insertion IBI tolerance as fraction of IBI mode (default=0.24)
        -s<flt> specify spike detector criterion in s.d. units (default=3.00)
        -a<flt> specify number of seconds used to estimate initial IBI (default=6.00)
        -R      revise beat placements by reading the diag file
                use -R in conjunction with -i or -d options, -q is implied
N.B.:   option -k has no effect if exemplar beat start time (-t) given

EKG_datfile	ASCII text file listing of EKG trace. To be parsed into cardiac beats according to the algorithm
		outlined in main text section 2.1 and Appendix A.
		example file:
		04_1101_BOLD12_S_ldr_BP_1-30_FIR_48_db_zero_-0.3_2.7130_500_13_25_red_ascii_dump_EKG.dat
		EKG_datfile format:
		The first two columns give time (in sec) past the start of the record and EKG
		trace in microvolts. An optional third column can be used to code for the
		presence (1) or absence (0) of artifact. Currently this information us not used algorithmically
		but is included some EKG_phase output listings.
		In our implementation we generate EKG_datfile using an executable called NSepoch2power.
		The command is NSepoch2power -t<EKG_channel>, e.g., NSepoch2power -tEKG.
		Source code and a make file for NSepoch2powe are included in this distribution.

		EKG_phase always generates the following outputs:
*_template.dat	exemplar beat expanded into 3 dimensions
*_yz.dat	y(t) (Eqn 2.1) z(t) and second derivative (see main text Section 2.3) plus chi(t) (Appendix A)
		These signals are illustrated in Figure 3.
*_bps.dat	cardiac rate time series (beats per second) as illustrated in Figure 4A
*_IBIhist.dat	IBI histogram (as illustrated in Figure 4B
*_beatavg.dat	average beat time course (3 dimensions). 6 columns. The last 3 columns give the s.d. time course
		evaluated over all beats contributing to the averages. (Not illustrated in the paper)
*_ISE.tcl	tcl script that can be used by Neuroscan software to insert markers corresponding to cardiac beats
*_diag.txt	diagnostic listing of all detected beats. The last column codes manner of beat insertion.
			*	automatically inserted using stringent criteria (Appendix A)
			**	automatically inserted using second tier criteria
			***	manually inserted
		The *_diag.txt file is passed on to reduce_BKG.
		Manual override is achieved by specifying on the command line times of beats to delete
		and/or insert in concert with option -R.
		*_diag.txt quotes the EKG_phase command line.
N.B.:		Exemplar beats are identified via the -t option specifying the exemplar beat latency (in seconds)
		relative to the start of the record. This latency should be stated as 150 to 200 msec before the
		QRS complex. It an exemplar beat is not given, this version of EKG_phase will automatically generate
		an exemplar beat using an OBS-like strategy. Auto-exemplar finding generally works well if the cardiac
		trace is recorded using only chest leads directly connected to non-Neuroscan equipment (e.g., InVivo
		clinical monitoring equipment). 

reduce_BKG	
Usage:  reduce_BKG <NS_file.dat> <EKG_diagfile>
e.g.,   reduce_BKG 04_1004_ascii[.dat] 04_1004_EKG_diag[.txt]
        option
        -s<int> specify mGLM averaging interval in seconds (default=10.00)
        -h<int> specify mGLM basis function frequency limit in Hz (default=30.00)
        -A      reduce BKG using AAS (default use mGLM)
        -a<int> specify AAS number of averaged beats (default=11)
        -l<flt> specify model duration in seconds (default=1.0240)
        -g<flt> specify Hanning taper duration in sec (default=0.0320)
        -z<str> exclude artifact using specified layfile
        -k<flt> specify beat averaging timebase in seconds (use with -C) (default=1.5360)
        -f<int> specify FFT points (use with -P) (default=8192)
        -c<str> select (comma separated) channels for diagnostic output
        -N      output BKG reduced record in NeuroScan ASCII format
        -E      output estimated BKG      in NeuroScan ASCII format
        -t      test mode (create artificial square wave BKG in first channel)
        -v      verbose mode
        -d      debug mode
        -x      fast execution mode (suppress condition number computation)
        diagnostic output options for selected channels
        -C      EEG averaged in phase with beats pre/post reduction
        -T      EEG complete timecourse pre/post reduction
        -P      EEG spectral power density pre/post reduction
        -B      mean BKG model psd
        -M<int> save (at least) the first specified number of BKG models in groups of 25
        -L<int> average earliest and latest BKG models (and corresponding post-reduction EEG)

NS_file.dat	Neuroscan data file exported as ASCII text. reduce_BKG accepts files exported as epochs
		and as continuous format.
		example file:
		04_1101_BOLD12_S_ldr_BP_1-30_FIR_48_db_zero_-0.3_2.7130_500_13_25_red_epoch_ascii.dat
EKG_diagfile	This is the *_diag.txt file generated by EKG_phase.
		example file:
		04_1101_BOLD12_S_ldr_BP_1-30_FIR_48_db_zero_-0.3_2.7130_500_13_25_red_ascii_dump_EKG_diag.txt
layfile		ASCII text file containing a listing of epochs to be excluded from the BKG computations.
		specified via the -z option
		example file:
		04_1101_BOLD12_S_ldr_BP_1-30_FIR_48_db_zero_-0.3_2.7130_500_13_25_red.lay
		layfile format:
		critical data introduced by a line including the string "Comments"
		Each subsequent line gives the time (in sec) past the start of the record of artifact onset and
		artifact duration.
		Valid field separators are ',' and ' '.
		Each line must have 6 total fields but the last 4 are ignored. (This requirement
		can be easily disabled by a simple and obvious modification of the source code.)
		N.B.: The layfile mechanism is intended prevent brief epochs of corrupted EKG
 		from interfering with beat detection. Attempting to exclude extended parts of the record
		significantly slows down reduce_BKG execution. Such epochs should be excluded from subsequent
		processing using other mechanisms.
		
		main reduce_BKG outputs:
*BKGred_diag.dat 	Always generated. Lists time of each beat and associated algebraic parameters (e.g., mGLM
			condition number). *BKGred_diag.dat quotes the reduce_BKG command line.
		example file:
		04_1101_BOLD12_S_ldr_BP_1-30_FIR_48_db_zero_-0.3_2.7130_500_13_25_red_epoch_ascii_BKGred_diag.dat
*BKGred.dat		option -N	Neuroscan ASCII format record post BKG artifact reduction (e.g., main text Figure 5B).
*BKGred_eBKG.dat	option -E	Neuroscan ASCII format record of estimated BKG artifact (e.g., main text Figure 5C).
		reduce_BKG generates diagnostic outputs on channels selected via the -c option.
		Example option usage: -cO1,O2,C3,C4
		Specific types of diagnostic output are specified by additional options (partial listing follows).
		option -M<int>	BKG models in groups of 25 (e.g., main text Figure 6).
		option -C	EEG averaged in phase with beats pre/post reduction. (e.g., main text Figure 7).
		example file:
		04_1101_BOLD12_S_ldr_BP_1-30_FIR_48_db_zero_-0.3_2.7130_500_13_25_red_epoch_ascii_BKGred_O1_avg.dat
		option -L<int>	average earliest and latest BKG models (e.g., Supplemntary Figure 1A).
