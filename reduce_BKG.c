#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <float.h>
#include <limits.h>
#include <string.h>
#include <ctype.h>
#include <assert.h>
#include <unistd.h>

#define TWO		2		/* for storing pre/post info */
#define MAXL		256		/* filename character length */
#define MAXS		1024		/* input string length */
#define MAXT		8		/* character length of channel label */
#define MAXM		25		/* maximum number of models/file written with -M option (rooted in xmgr limits) */
#define LENTEMPS	1.024		/* mGLM model duration in seconds */
#define AVGLENS		1.536		/* averaged BKG timebase in seconds */
#define HANNINGS	0.032		/* Hanning taper duration in seconds */
#define HIFREQ		30.		/* basis function high frequency */
#define BEATWINS	10.		/* number of seconds to include in mGLM integrals */
#define NAVG		11		/* number of beats to average for AAS */
#define FFTPTS		8192		/* eeg psd pts */
#define PAD		2.		/* assumed artifact duration in sec at begining and end of record */

/***********/
/* externs */
/***********/
extern void	sinbas_ (float *b, int *npts, int *nbasis);
extern int	npad_ (int *nsample, int *margin);					/* FORTRAN librms */
extern void	deigen_ (double *A, double *W, int *ndim);				/* FORTRAN librms */
extern void	dmatinv_ (double *A, int *ndim, double *det);				/* FORTRAN librms */
extern void	fft_   (float *a, float *b, int *nseg, int *n, int *nspn, int *idir);	/* FORTRAN librms */
extern void	realt_ (float *a, float *b, int *nseg, int *n, int *nspn, int *idir);	/* FORTRAN librms */

typedef struct {
	int	i;		/* time point in EKG record corresponding to start of beat */
	int	j;	 	/* time point in EKG record corresponding to start of template */
	float	signal;		/* strength of signal used to detect beat */
	int	art;		/* presence of artifact on beat model (lentempp) */
	int	exclude;	/* presence of artifact on averaging interval (avglenp) */
	float	**bkgfit; 	/* the model for each beat bkgfit[MAXC][LENTEMPP] */
	int	mstart, mend;	/* points where the beat window should begin & end (dependent on artifacts & beatwinp */
	int	rbeatwinp; 	/* revised beatwinp after considering artifact */
	double  det;		/* determinant for inverted A matrix */
	double	cndnum;		/* condition number for A matrix */
} BEAT;

typedef struct {
	char		trode[MAXT];
	float		mean;
	int		diag;
} CHANNEL;

/********************/
/* global variables */
/********************/
	char		rcsid[] = "$Id: reduce_BKG.c,v 1.29 2006/12/17 01:16:02 avi Exp $";
	char		program[MAXL];
	int		debug = 0;
	int		verbose = 0, fast = 0;

/*********************/
/* NS file variables */
/*********************/
	char		NStime[MAXL], NSdate[MAXL];	/* N.B.: time and date are system variables */
	char		rows[MAXL];
	int		channels;
	float		rate;			/* sampling rate in Hz */

/*******/
/* eeg */
/*******/
	CHANNEL		*channel;
	int 		points, sweeps, npts, nchan;
	float 	  	**eeg;			/* eeg[ichan][ipoint] */
	float 	  	**bcg;			/* the estimated ballistocardiogram */
	short		*artifact;		/* true or false for every point on record */
	float		***timecourse;

/*************/
/* BKG model */
/*************/
	BEAT		*beat;
	int		nbeat;
	float		lentemps = LENTEMPS;
	int		lentempp;
	int		lentempp_pad;		/* for BKG psd */
	float		*hanning;
	float		hannings = HANNINGS;
	int		hanningp;

/*******/
/* GBR */
/*******/
	float		**basis;
	int		nbasis, nbasis1;
	float		beatwins = BEATWINS;	/* mGLM integration interval */
	int		beatwinp, beatwinpo2;
	float		hifreq = HIFREQ;	/* frequency cutoff of BKG model */

/*******/
/* AAS */
/*******/
	int		navg = NAVG, navgo2;

/******************/
/* beat averaging */
/******************/
	float		***beatavg, ***beatvar;
	float		avglens = AVGLENS;
	int		avglenp;

/*******/
/* FFT */
/*******/
	int		fftpts = FFTPTS;
	float		hzpbin;
	float		***eegpsd;

void errm (char* program) {
	fprintf (stdout, "%s: memory allocation error\n", program);
	exit (-1);
}

void errr (char* program, char* filespc) {
	fprintf (stdout, "%s: %s read error\n", program, filespc);
	exit (-1);
}

void errf (char* program, char* filespc) {
	fprintf (stdout, "%s: %s format error\n", program, filespc);
	exit (-1);
}

void errw (char* program, char* filespc) {
	fprintf (stdout, "%s: %s write error\n", program, filespc);
	exit (-1);
}

void getroot (char *filespc, char *filroot) {
	char	*str;
	strcpy (filroot, filespc);
	while (str = strrchr (filroot, '.')) {
			if (!strcmp (str, ".rec"))	*str = '\0';
		else	if (!strcmp (str, ".dat"))	*str = '\0';
		else	if (!strcmp (str, ".txt"))	*str = '\0';
		else	if (!strcmp (str, ".lay"))	*str = '\0';
		else	break;
	}
}

int split (char *string, char *srgv[], int maxp) {
	int	i, m;
	char	*ptr;

	if (ptr = strchr (string, '#')) *ptr = '\0';
/**********************************/
/* convert '[', ']', ',' to space */
/**********************************/
	ptr = string; while (ptr = strpbrk (ptr, "[],")) *ptr = '\x20';
	i = m = 0;
	while (m < maxp) {
		while (!isgraph ((int) string[i]) && string[i]) i++;
		if (!string[i]) break;
		srgv[m++] = string + i;
		while (isgraph ((int) string[i])) i++;
		if (!string[i]) break;
		string[i++] = '\0';
	}
	return m;
}

void write_command_line (FILE *outfp, int argc, char *argv[], int gbr) {
	int		i;

	fprintf (outfp, "#%s", program);
	for (i = 1; i < argc; i++) fprintf (outfp, " %s", argv[i]);
	fprintf (outfp, "\n");
	fprintf (outfp, "#%s\n", rcsid);
	if (gbr) {
		fprintf (outfp, "#BKG reduction by mGLM integration interval %.2f sec\n", beatwins);
		fprintf (outfp, "#freq limit %.4f Hz; basis functions %d (including DC)\n", hifreq, nbasis1);
	} else {
		fprintf (outfp, "#BKG reduction by AAS %d beats in average\n", navg);
	}
	fprintf (outfp, "#BKG model duration %.4f sec; Hanning taper %.4f sec\n", lentemps, hannings);
}

void writeNS_file (float **array, char *outfile) {
	FILE		*outfp;
	int		i, ichan;

	fprintf (stdout, "#Writing: %s\n", outfile);
	if (!(outfp = fopen (outfile, "w"))) errw (program, outfile);
	fprintf (outfp, "[Subject]\t\n");
	fprintf (outfp, "[Date]\t%s\n", NSdate);
	fprintf (outfp, "[Time]\t%s\n", NStime);
	fprintf (outfp, "[Channels]\t%i\n", channels);
	fprintf (outfp, "[Rate]\t%10.6f\n", rate);
	fprintf (outfp, "[Type]\tContinuous\n");
	fprintf (outfp, "[Rows]\t%s\n", rows);
	fprintf (outfp, "[Electrode Labels]\n");
	for (i = 0; i < nchan; i++) fprintf (outfp, "[      %s]\t", channel[i].trode);
	fprintf (outfp, "\n");
	fprintf (outfp, "[Electrode XUnits]\n");
	for (i = 0; i < nchan; i++) fprintf (outfp, "[ Default]\t");
	fprintf (outfp, "\n");
	fprintf (outfp, "[Electrode YUnits]\n");
	for (i = 0; i < nchan; i++) fprintf (outfp, "[ Default]\t");
	fprintf (outfp, "\n");
	fprintf (outfp, "[Continuous Data]\n");
	for (i = 0; i < npts; i++) { 
		for (ichan = 0; ichan < nchan; ichan++) {
			fprintf (outfp, "%10.4f\t", array[ichan][i]);
		}
		fprintf (outfp, "\n");
	}
	if (fclose (outfp)) errw (program, outfile);
}

double **calloc_double2 (int n1, int n2) {
	int	i;
	double	**a;

	if (!(a = (double **) malloc (n1 * sizeof (double *)))) errm (program);
	if (!(a[0] = (double *) calloc (n1 * n2, sizeof (double)))) errm (program);
	for (i = 1; i < n1; i++) a[i] = a[0] + i*n2;
	return a;
}

void free_double2 (double **a) {
	free (a[0]);
	free (a);
}

float **calloc_float2 (int n1, int n2) {
	int	i;
	float	**a;

	if (!(a = (float **) malloc (n1 * sizeof (float *)))) errm (program);
	if (!(a[0] = (float *) calloc (n1 * n2, sizeof (float)))) errm (program);
	for (i = 1; i < n1; i++) a[i] = a[0] + i*n2;
	return a;
}

void free_float2 (float **a) {
	free (a[0]);
	free (a);
}

float ***calloc_float3 (int n1, int n2, int n3) {
	unsigned int	i, j;
	float		***a;

	if (!(a = (float ***) malloc (n1 * sizeof (float **)))) errm (program);
	if (!(a[0] = (float **) malloc (n1 * n2 * sizeof (float *)))) errm (program);
	if (!(a[0][0] = (float *) calloc (n1 * n2 * n3, sizeof (float)))) errm (program);
	for (i = 0; i < n1; i++) {
		a[i] = a[0] + n2*i;
		for (j = 0; j < n2; j++) {
			a[i][j] = a[0][0] + n3*(n2*i + j);
		}
	}
	return a;
}

void free_float3 (float ***a) {
	free (a[0][0]);
	free (a[0]);
	free (a);
}

double zeromean (float *f, int n) {
	int		i;
	double		u;

	for (u = i = 0; i < n; i++) u += f[i];
	u /= n;
	for (i = 0; i < n; i++) f[i] -= u;
	return u;
}

void mGLM (int ibeat) {
	float		**bcomp, *coeff, *B;
	double		**W, **E, **A, det;
	int 		k, j, i, ii, mpts, nright, nleft;
	int		kbeat, mstart, mend, ichan;

/*********************************************************************/
/* mend in mGLM() is first timepoint excluded on right from integral */
/*********************************************************************/
	mstart = beat[ibeat].i; if (mstart < 0) mstart = 0;
	mend = mstart;
	nright = 0; nleft = !artifact[mstart];
	do {
		if (mstart <= 0) break;
		mstart--;
		if (!artifact[mstart]) nleft++;
	} while (nleft <= beatwinpo2);
	do {
		if (mend >= npts) break;
		mend++;
		if (!artifact[mend]) nright++;
	} while (nright + nleft <= beatwinp);
	  while (nright + nleft <= beatwinp) {
		if (mstart <= 0) break;
		mstart--;
		if (!artifact[mstart]) nleft++;
	}
	beat[ibeat].mend = mend; beat[ibeat].mstart = mstart;
	beat[ibeat].rbeatwinp = mend - mstart;

	bcomp  = calloc_float2  (nbasis1, beat[ibeat].rbeatwinp);
	E      = calloc_double2 (nbasis1, nbasis1);
	W      = calloc_double2 (nbasis1, nbasis1);
	A      = calloc_double2 (nbasis1, nbasis1);
	for (kbeat = mpts = 0; kbeat < nbeat; kbeat++) {
		ii = beat[kbeat].i;
		if (ii + lentempp <  beat[ibeat].mstart) continue;
		if (ii 	          >= beat[ibeat].mend)   break;
		for (i = 0; i < lentempp; i++, ii++) {
			if (ii <  beat[ibeat].mstart) continue;
			if (ii >= beat[ibeat].mend)   break;
			if (artifact[ii]) continue;
			for (j = 0; j < nbasis1; j++) {
				bcomp[j][ii - beat[ibeat].mstart] += basis[j][i];
			}
			mpts++;
		}
	}
	for (j = 0; j <= nbasis; j++) {
		for (i = 0; i <= nbasis; i++) {
			for (A[i][j] = k = 0; k < beat[ibeat].rbeatwinp; k++) {
				A[i][j] += bcomp[i][k]*bcomp[j][k];
			}
		}
	}

/***************************/
/* normalize linear system */
/***************************/
	for (j = 0; j <= nbasis; j++) {
		for (i = 0; i <= nbasis; i++) {
			A[i][j] /= mpts;
			E[i][j] = A[i][j];
		}
	}

/************************/
/* invert linear system */
/************************/
	if (!fast) {
		deigen_ (E[0], W[0], &nbasis1);
		beat[ibeat].cndnum = E[0][0]/E[nbasis][nbasis];
	}
	dmatinv_ (A[0], &nbasis1, &det);
	beat[ibeat].det = det;
	if (verbose) printf ("#mGLM: beat %3d cndnum %.4e det %.4e rbeatwinp %d\n",
			ibeat, beat[ibeat].cndnum, det, beat[ibeat].rbeatwinp);
	if (beat[ibeat].cndnum > 400) {
		printf ("%s: beat %3d GLM condition number is high (%.4e)\n", program, ibeat, beat[ibeat].cndnum);
		printf ("Consider reducing the modeling frequency limit or increasing the integration interval\n");
	}

/*******************************/
/* solve mGLM for each channel */
/*******************************/
	if (!(coeff = (float *) calloc (nbasis1, sizeof (float)))) errm (program);
	if (!(B     = (float *) calloc (nbasis1, sizeof (float)))) errm (program);
	for (ichan = 0; ichan < nchan; ichan++) {
		for (j = 0; j <= nbasis; j++) {
			for (B[j] = k = 0; k < beat[ibeat].rbeatwinp; k++) {
/*************************************************************************/
/* bcomp[j][k] is zero where artifact[k + beat[ibeat].mstart] is nonzero */
/*************************************************************************/
				B[j] += eeg[ichan][k + beat[ibeat].mstart]*bcomp[j][k];
			}
			B[j] /= mpts;
		}
		for (j = 0; j <= nbasis; j++) {
			for (coeff[j] = i = 0; i <= nbasis; i++) coeff[j] += A[i][j]*B[i];
		}
		for (i = 0; i < lentempp; i++) {
			for (beat[ibeat].bkgfit[ichan][i] = j = 0; j <= nbasis; j++) {
				beat[ibeat].bkgfit[ichan][i] += coeff[j]*basis[j][i];
			}
		}
	}

	free (B); free (coeff);
	free_float2 (bcomp);
	free_double2 (E); free_double2 (W); free_double2 (A);
}

void AAS (int ibeat, int ichan) {
	int		i, j, k, ii;
	int		kbeat, mstart, mend, nleft, nright;

/***********************************************/
/* mend in AAS() is last beat added to average */
/***********************************************/
	mstart = mend = ibeat;
	nright = 0; nleft = !beat[mstart].art;
	do {
		if (mstart <= 0) break;
		mstart--;
		if (!beat[mstart].art) nleft++;
	} while (nleft <= navgo2);
	do {
		if (mend >= nbeat - 1) break;
		mend++;
		if (!beat[mend].art) nright++;
	} while (nright + nleft < navg);
	  while (nright + nleft < navg) {
		if (mstart <= 0) break;
		mstart--;
		if (!beat[mstart].art) nleft++;
	}

	assert (mstart >= 0 && mend < nbeat);
	beat[ibeat].mstart = mstart;
	beat[ibeat].mend   = mend;
	k = 0; for (kbeat = beat[ibeat].mstart; kbeat <= beat[ibeat].mend; kbeat++) {
		if (beat[kbeat].art) continue;
		if (beat[kbeat].i < 0) continue;
		if (beat[kbeat].i + lentempp >= npts) continue;
		for (ii = beat[kbeat].i, i = 0; i < lentempp; i++, ii++) {
			beat[ibeat].bkgfit[ichan][i] += eeg[ichan][ii];
		}
		k++;
	}
	if (k != navg) printf ("AAS: ibeat=%d\tk=%d\tnavg=%d\tnleft=%d\tnright=%d\tmstart=%d\tmend=%d\n",
		ibeat, k, navg, nleft, nright, mstart, mend);
	for (i = 0; i < lentempp; i++) beat[ibeat].bkgfit[ichan][i] /= k;
}

/******************************************************************/
/* compile beat timecourse average and sd for channel of interest */
/******************************************************************/
int avgsd (int task, int ichan, int mchan) {
	int 		i, n, ibeat;
	double 		q;

	n = 0;
	for (ibeat = 0; ibeat < nbeat; ibeat++) {
		if (beat[ibeat].exclude) continue;
		if (beat[ibeat].i + avglenp >= npts) break;
		for (i = 0; i < avglenp; i++) {
			q = eeg[ichan][beat[ibeat].i + i];
			beatavg[task][mchan][i] += q;
			beatvar[task][mchan][i] += q*q;
		}
		n++;
	}
	for (i = 0; i < avglenp; i++) {
		beatavg[task][mchan][i] /= n; beatvar[task][mchan][i] /= (n - 1);
		beatvar[task][mchan][i] -= beatavg[task][mchan][i]*beatavg[task][mchan][i];
	}
	for (i = 0; i < npts; i++) {
		timecourse[task][mchan][i] = eeg[ichan][i];
	}
	return n;
}

/***********/
/* FFT BKG */
/***********/
void bkgfft (float *bkgmean, float *bkgpsd) {
	int		i, j, no2, nbin;
	int		one = 1, negone = -1;
	float		*a, *b, *bkgpad;
	double		q;

	no2		= lentempp_pad/2;
	nbin		= no2 + 1;
	if (!(bkgpad	= (float *) calloc (lentempp_pad, sizeof (float)))) errm (program);
	if (!(a		= (float *) calloc (nbin, sizeof (float)))) errm (program);
	if (!(b		= (float *) calloc (nbin, sizeof (float)))) errm (program);

	for (i = 0; i < lentempp; i++) bkgpad[i] = bkgmean[i];
	q = (bkgpad[0] - bkgpad[lentempp - 1])/(lentempp_pad - lentempp + 1);
	for (i = lentempp; i < lentempp_pad; i++) {
		bkgpad[i] = bkgpad[lentempp - 1] + q*(i - lentempp + 1);
	}
	for (j = i = 0; j < no2; j++) {
		a[j] = bkgpad[i++];
		b[j] = bkgpad[i++];
	}
	fft_   (a, b, &one, &no2, &one, &negone);
	realt_ (a, b, &one, &no2, &one, &negone);
	q = 1./(lentempp_pad*lentempp_pad);
	bkgpsd[0] = a[0]*a[0]*q;
	for (i = 1; i < no2; i++) bkgpsd[i] = 2.*(a[i]*a[i] + b[i]*b[i])*q;
	bkgpsd[no2] = a[no2]*a[no2]*q;
	free (bkgpad); free (a); free (b);
}

/***********/
/* FTT EEG */
/***********/
#define MAXE	128		/* for parsing state of artifact[] */
int eegfft (int ichan, float *psdout) {
	int		i, j, no2, nbin;
	int		one = 1, negone = -1;
	float		*a, *b;
	int		counter, k, m, start[MAXE], stop[MAXE];

	i = m = 0;
	while (m < MAXE) {
		while ( artifact[i]) {i++; if (i >= npts) break;}
		if (i >= npts) break;
		start[m] = i;
		while (!artifact[i]) {i++; if (i >= npts) break; if (i - start[m] == fftpts) break;}
		stop[m++] = i;
		if (i >= npts) break;
	}
if (0) {
	for (k = 0; k < m; k++) {
		printf ("%10d start%10d stop%10d %10d\n", k, start[k], stop[k], stop[k] - start[k]);
	}
	exit (0);
}
	no2 = fftpts/2;
	nbin = no2 + 1;
	if (!(a = (float *) calloc (nbin, sizeof (float)))) errm (program);
	if (!(b = (float *) calloc (nbin, sizeof (float)))) errm (program);

	for (counter = k = 0; k < m; k++) if (stop[k] - start[k] == fftpts) {
		i = start[k];
		for (j = 0; j < no2; j++) {
			a[j] = eeg[ichan][i++];
			b[j] = eeg[ichan][i++];
		}
		fft_   (a, b, &one, &no2, &one, &negone);
		realt_ (a, b, &one, &no2, &one, &negone);
		psdout[0] += a[0]*a[0];
		for (i = 1; i < no2; i++) psdout[i] += 2.*(a[i]*a[i] + b[i]*b[i]);
		psdout[no2] += a[no2]*a[no2];
		counter++;
	}

	for (i = 0; i < nbin; i++) psdout[i] /= (fftpts*fftpts*counter);
	free (a); free (b);
	return counter;
}

void usage (char *program) {
	fprintf (stdout, "Usage:\t%s <NS_file.dat> <EKG_diagfile>\n", program);
	fprintf (stdout, "e.g.,\t%s 04_1004_ascii[.dat] 04_1004_EKG_diag[.txt]\n", program);
	fprintf (stdout, "\toption\n");
	fprintf (stdout, "\t-s<int>\tspecify mGLM averaging interval in seconds (default=%.2f)\n", BEATWINS);
	fprintf (stdout, "\t-h<int>\tspecify mGLM basis function frequency limit in Hz (default=%.2f)\n", HIFREQ);
	fprintf (stdout, "\t-A\treduce BKG using AAS (default use mGLM)\n");
	fprintf (stdout, "\t-a<int>\tspecify AAS number of averaged beats (default=%i)\n", NAVG);
	fprintf (stdout, "\t-l<flt>\tspecify model duration in seconds (default=%.4f)\n", LENTEMPS);
	fprintf (stdout, "\t-g<flt>\tspecify Hanning taper duration in sec (default=%.4f)\n", HANNINGS);
	fprintf (stdout, "\t-z<str>\texclude artifact using specified layfile\n");
	fprintf (stdout, "\t-k<flt>\tspecify beat averaging timebase in seconds (use with -C) (default=%.4f)\n", AVGLENS);
	fprintf (stdout, "\t-f<int>\tspecify FFT points (use with -P) (default=%d)\n", FFTPTS);
	fprintf (stdout, "\t-c<str>\tselect (comma separated) channels for diagnostic output\n");
	fprintf (stdout, "\t-N\toutput BKG reduced record in NeuroScan ASCII format\n");	
	fprintf (stdout, "\t-E\toutput estimated BKG      in NeuroScan ASCII format\n");
	fprintf (stdout, "\t-t\ttest mode (create artificial square wave BKG in first channel)\n");
	fprintf (stdout, "\t-v\tverbose mode\n");
	fprintf (stdout, "\t-d\tdebug mode\n");
	fprintf (stdout, "\t-x\tfast execution mode (supress condition number computation)\n");
	fprintf (stdout, "\tdiagnostic output options for selected channels\n");
	fprintf (stdout, "\t-C\tEEG averaged in phase with beats pre/post reduction\n");
	fprintf (stdout, "\t-T\tEEG complete timecourse pre/post reduction\n");
	fprintf (stdout, "\t-P\tEEG spectral power density pre/post reduction\n");
	fprintf (stdout, "\t-B\tmean BKG model psd\n");
	fprintf (stdout, "\t-M<int>\tsave (at least) the first specified number of BKG models in groups of %d\n", MAXM);
	fprintf (stdout, "\t-L<int>\taverage earliest and latest BKG models (and corresponding post-reduction EEG)\n");
	exit (1);
}

int main (int argc, char *argv[]) {
	char		*ptr, command[MAXL], string[MAXS], *srgv[MAXL];
	int		c, i, j, k, l, m, n, ii;

/************/
/* file i/o */
/************/
	FILE		*datfp, *outfp, *diagfp, *layfp;
	char		datroot[MAXL], datfile[MAXL], outroot[MAXL], outfile[MAXL], diagfile[MAXL], diagroot[MAXL];
	char		layroot[MAXL], layfile[MAXL];
	char		dchanstr[MAXL] = "";
	char		trailer[] = "BKGred";
	fpos_t		filepos;	/* for parsing datfile */

/*****************/
/* EEG variables */
/*****************/
	int		ichan, mchan, isweep;
	float		delta;

/************/
/* artifact */
/************/
	int 		*zempart, *zempdur;

/*******/
/* EKG */
/*******/
	int		ibeat, ibeat0, ibeat1, mbeat, kbeat, nbeat_in_avg;

/*******/
/* BKG */
/*******/
	float		**bkgmean, *bkgpsd;

/*******/
/* FFT */
/*******/
	int		zero = 0;
	int		nepoch;			/* number of epochs included in pre/post EEG psd */

/*********/
/* flags */
/*********/
	int		continuous = 0;
	int		timedomain = 0;
	int		headers = 0;
	int		zempel = 0;
	int		test = 0;
	int		gbr = 1, aas = 0;
	int		outputN = 0;
	int		outputP = 0;
	int		outputB = 0;
	int		outputE = 0;
	int		outputC = 0;
	int		outputT = 0;
	int		outputL = 0;
	int		outputM = 0;

	printf ("#%s\n", rcsid);
	if (!(ptr = strrchr (argv[0], '/'))) ptr = argv[0]; else ptr++;
	strcpy (program, ptr);
	printf ("#%s", program); for (i = 1; i < argc; i++) printf (" %s", argv[i]); printf ("\n");

/************************/
/* process command line */
/************************/
	for (k = 0, i = 1; i < argc; i++) {
		if (*argv[i] == '-') {
			strcpy (command, argv[i]); ptr = command;
			while (c = *ptr++) switch (c) {
				case 'x': fast++; 					break;
				case 't': test++; 					break;
				case 'd': debug++;					break;
				case 'v': verbose++;					break;
				case 'A': aas++; gbr = 0;				break;
				case 'G': gbr++; aas = 0;				break;
				case 'B': outputB++;					break;
				case 'C': outputC++;					break;
				case 'E': outputE++;					break;
				case 'T': outputT++; 					break;
				case 'P': outputP++; 					break;
				case 'N': outputN++;					break;
				case 'L': outputL = atoi (ptr);				*ptr = '\0'; break;
				case 'M': outputM = atoi (ptr);				*ptr = '\0'; break;
				case 's': beatwins = atof (ptr);			*ptr = '\0'; break;
				case 'h': hifreq = atof (ptr);				*ptr = '\0'; break;
				case 'l': lentemps = atof (ptr);			*ptr = '\0'; break;
				case 'k': avglens = atof (ptr);				*ptr = '\0'; break;
				case 'a': navg = atoi (ptr);				*ptr = '\0'; break;
				case 'c': strcpy (dchanstr, ptr);	    		*ptr = '\0'; break;
				case 'z': zempel++; getroot (ptr, layroot); 		*ptr = '\0'; break;
				case 'f': fftpts = atoi (ptr);				*ptr = '\0'; break;
			}
		}
		else switch (k) {
			case 0: getroot (argv[i], datroot);	k++; break;
			case 1:	getroot (argv[i], diagroot);	k++; break;
		}	
	}
	if (k < 2) usage (program);

/***********************/
/* parse input NS file */
/***********************/
	sprintf (datfile, "%s.dat", datroot);
	fprintf (stdout, "#Reading: %s\n", datfile);
	if (!(datfp = fopen (datfile, "r"))) errr (program, datfile);
	while (fgets (string, MAXS, datfp)) {
		m = split (string, srgv, MAXL);
		if (!strcmp (srgv[0], "Date"))		strcpy (NSdate, srgv[1]);
		if (!strcmp (srgv[0], "Time"))		strcpy (NStime, srgv[1]);
		if (!strcmp (srgv[0], "Channels"))	channels = atoi (srgv[1]);
		if (!strcmp (srgv[0], "Rate"))		rate = atof (srgv[1]);
		if (!strcmp (srgv[0], "Headers")	&& !strcmp (srgv[1], "Yes")) headers++;
		if (!strcmp (srgv[0], "Points"))	points = atoi (srgv[1]);
		if (!strcmp (srgv[0], "Sweeps"))	sweeps = atoi (srgv[1]);
		if (!strcmp (srgv[0], "Rows"))		strcpy (rows, srgv[1]);		
		if (!strcmp (srgv[0], "Type")		&& !strcmp (srgv[1], "Continuous")) continuous++;		
		if (!strcmp (srgv[0], "Domain")		&& !strcmp (srgv[1], "Time")) timedomain++;
		if (!strcmp (srgv[0], "Electrode")	&& !strcmp (srgv[1], "Labels")) {
			fgets (string, MAXS, datfp);
			nchan = split (string, srgv, MAXL);
			if (nchan != channels) errf (program, datfile);
			if (!(channel = (CHANNEL *) calloc (nchan, sizeof (CHANNEL)))) errm (program);
			for (i = 0; i < nchan; i++) strcpy (channel[i].trode, srgv[i]);
		}
		if (!strcmp (srgv[0], "Epoch")		&& !strcmp (srgv[1], "Header")) break;
		if (!strcmp (srgv[0], "Continuous")	&& !strcmp (srgv[1], "Data")) break;
	}
/************************************************/
/* determine total number of time points (npts) */
/************************************************/
	if (continuous) {
		if (fgetpos (datfp, &filepos)) errr (program, datfile);
		npts = 0; while (fgets (string, MAXS, datfp)) npts++;
		if (fsetpos (datfp, &filepos)) errr (program, datfile);
	} else {
		if (!headers) {
			fprintf (stdout, "%s: %s does not contain headers\n", program, datfile);
			fprintf (stdout, "*** select Header Option when exporting to ASCII ***\n");
			exit (-1);
		}
		if (!timedomain) {
			fprintf (stdout, "%s: %s not time domain data\n", program, datfile);
			exit (-1);
		}
		npts = points*sweeps;
	}

/*************************************************/
/* count channels selected for diagnostic output */
/*************************************************/
	m = split (dchanstr, srgv, MAXL);
	for (ichan = 0; ichan < nchan; ichan++) for (k = 0; k < m; k++) {
		if (!strcmp (channel[ichan].trode, srgv[k])) channel[ichan].diag++;
	}
	for (mchan = ichan = 0; ichan < nchan; ichan++) if (channel[ichan].diag) mchan++;

/***********************/
/* allocate EEG memory */
/***********************/
	delta		= 1./rate;
	lentempp	= rate*lentemps + 0.5;
	hanningp	= rate*hannings + 0.5;
	beatwinp	= rate*beatwins + 0.5; 
	avglenp		= rate*avglens  + 0.5;
	lentempp_pad	= npad_ (&lentempp, &zero);	/* for BKG psd */
	beatwinpo2	= beatwinp/2;
	navgo2		= navg/2;

	eeg		= calloc_float2 (nchan, npts);
	bcg		= calloc_float2 (nchan, npts);
	beatavg		= calloc_float3 (TWO, mchan, avglenp);
	beatvar		= calloc_float3 (TWO, mchan, avglenp);
	timecourse	= calloc_float3 (TWO, mchan, npts);
	eegpsd		= calloc_float3 (TWO, mchan, fftpts/2 + 1);
	bkgmean		= calloc_float2 (4, lentempp);
	if (!(bkgpsd	= (float *) calloc (lentempp_pad/2 + 1,	sizeof (float))))	errm (program);
	if (!(zempart	= (int *)   calloc (MAXL,		sizeof (int))))		errm (program);
	if (!(zempdur	= (int *)   calloc (MAXL,		sizeof (int))))		errm (program);
	if (!(hanning	= (float *) calloc (hanningp,		sizeof (float))))	errm (program);
	if (!(artifact	= (short *) calloc (npts,		sizeof (short))))	errm (program);

/****************************/
/* read data into eeg array */
/****************************/
	if (continuous) {
		for (i = 0; i < npts; i++) {
			fgets (string, MAXS, datfp);
			if (split (string, srgv, MAXL) != nchan) errf (program, datfile);
			for (ichan = 0; ichan < nchan; ichan++) eeg[ichan][i] = atof (srgv[ichan]);
		}
	} else {
		rewind (datfp);
		for (isweep = k = 0; isweep < sweeps; isweep++) {
			do {
				fgets (string, MAXS, datfp);
			} while (!strstr (string, "Epoch Header"));
			do {
				fgets (string, MAXS, datfp);
			} while (!strstr (string, "Epoch Data"));
			for (i = 0; i < points; i++) {
				fgets (string, MAXS, datfp);
				if (split (string, srgv, MAXL) != nchan) errf (program, datfile);
				for (ichan = 0; ichan < nchan; ichan++) eeg[ichan][k] = atof (srgv[ichan]);
				k++;
			}
		}
		if (k != npts) errf (program, datfile);
	}
	if (fclose (datfp)) errr (program, datfile);
	printf ("#npts=%d\n", npts);

/****************************/
/* read EKG_phase diag file */
/****************************/
	nbeat = ibeat = 0;
	sprintf (diagfile, "%s.txt", diagroot);
	if (!(diagfp = fopen (diagfile, "r"))) errr (program, diagfile);
	while (fgets (string, MAXS, diagfp)) {
		m = split (string, srgv, MAXL);
		if (!m || (!strcmp(srgv[0], "delta"))) continue;
		nbeat++;
	}
	rewind (diagfp);
	if (!(beat = (BEAT *) calloc (nbeat, sizeof (BEAT)))) errm (program);
	for (ibeat = 0; ibeat < nbeat; ibeat++) {
		beat[ibeat].bkgfit  = calloc_float2 (nchan, lentempp);
	}
	ibeat = 0;
	while (fgets (string, MAXS, diagfp)) {
		m = split (string, srgv, MAXL);
		if (!m) continue;
		ibeat = (atoi(srgv[0]) - 1);
		beat[ibeat].i = (atof(srgv[1])*rate + 0.01) - hanningp;	/* shift all beats to the left by hanning taper */
		if (beat[ibeat].i < 0) beat[ibeat].i = 0;		/* keep beat in memory domain */
		beat[ibeat].signal = atof(srgv[3]);
		if (m == 13) {
			for (k = i = 0; i < strlen(srgv[12]); i++) if (srgv[12][i] == '*') k--;
			beat[ibeat].j = k;
		}
	}
	if (fclose (diagfp)) errr (program, diagfile);

if (test) {
/*************************/
/* set up artificial BKG */
/*************************/
	strcpy (channel[0].trode, "test");	
	strcpy (dchanstr, channel[0].trode);
	for (i = 0; i < npts; i++) eeg[0][i] = 0.;
	for (k = 0; k < nbeat; k++) {
		ii = beat[k].i + 50;
		for (i = 0; i < 50; i++, ii++) {
			if (ii < 0) continue;
			if (ii >= npts) break;
			eeg[0][ii] += 50.;
		}
	}
}

if (zempel) {
/************************************************/
/* parse John Zempel's artifact assessment file */
/************************************************/
	sprintf (layfile, "%s.lay", layroot);
	fprintf (stdout, "#Reading: %s\n", layfile);
	if (!(layfp = fopen (layfile, "r"))) errr (program, layfile);
	i = 0;
	do {
		fgets (string, MAXS, layfp);
	} while (!strstr (string, "Comments"));
	fgets (string, MAXS, layfp);
	do {
		m = split (string, srgv, MAXL);
		if (m != 6) {
			fprintf (stdout, "%s: %s format error (number of fields in line != 6)\n", program, layfile);
			exit (-1);
		}
		zempart[i] = (rate * atof (srgv[0]));
		zempdur[i] = (rate * atof (srgv[1])) + 0.5;	
		i++;
	} while (fgets (string, MAXS, layfp));
	for (k = 0; k < i; k++) {
		for (l = zempart[k]; l < (zempart[k] + zempdur[k]); l++) {
			if (l < 0) continue;
			if (l >= npts) break;
			artifact[l] = 1;
		}
	}
	if (fclose (layfp)) errr (program, layfile);
}

/**********************************************************************/
/* always set artifact for PAD seconds at beginning and end of record */
/**********************************************************************/
	for (k = 0; k < PAD*rate; k++) artifact[k] = 1;
	for (k = (npts - 1); k > ((npts - 1) - PAD*rate); k--) artifact[k] = 1;

/****************************************/
/* search for beats containing artifact */
/****************************************/
	for (ibeat = 0; ibeat < nbeat; ibeat++) {
		for (i = beat[ibeat].i; i < beat[ibeat].i + lentempp; i++) {
			if (i < 0) continue; if (i >= npts) break;
			if (artifact[i]) beat[ibeat].art++;
		}
		for (i = beat[ibeat].i; i < beat[ibeat].i + avglenp; i++) {
			if (i < 0) continue; if (i >= npts) break;
			if (artifact[i]) beat[ibeat].exclude++;
		}
	}

/*******************************/
/* make each channel zero mean */
/*******************************/
	printf ("%10s%10s%10s%10s\n", "channel", "electrode", "mean", "diagout");
	for (ichan = 0; ichan < nchan; ichan++) {
		channel[ichan].mean = zeromean (eeg[ichan], npts);
		printf ("%10i%10s%10.4f%10d\n", ichan + 1, channel[ichan].trode, channel[ichan].mean, channel[ichan].diag);
	}

/******************************************************/
/* pre BKG reduction diagnostics on selected channels */
/******************************************************/
for (mchan = ichan = 0; ichan < nchan; ichan++) if (channel[ichan].diag) {
	if (outputC) avgsd (0, ichan, mchan);
	if (outputP) eegfft (ichan, eegpsd[0][mchan]);
	mchan++;
}

/*************************/
/* set up hanning window */
/*************************/
	for (i = 0; i < hanningp; i++) {
		hanning[i] = .5*(1 + cos((M_PI*(i + 1))/(hanningp)));
	}

if (gbr) {
/**************************/
/* set up basis functions */
/**************************/
	nbasis = 2 * (int) (lentemps*hifreq + 0.5); nbasis1 = nbasis + 1;
	printf ("#number of basis functions = %d (including DC)\n", nbasis1);
	basis = calloc_float2 (nbasis1, lentempp);
	sinbas_ (basis[0], &lentempp, &nbasis);

/****************************/
/* solve mGLM for each beat */
/****************************/
	for (ibeat = 0; ibeat < nbeat; ibeat++) {
		if (beat[ibeat].i <  beatwinpo2) continue;
		if (beat[ibeat].i >= beat[nbeat - 1].i - beatwinpo2) break;
		mGLM (ibeat);
	}
	free_float2 (basis);

	for (ichan = 0; ichan < nchan; ichan++) {
/******************************************/
/* eliminate discontinuities with hanning */
/******************************************/
		for (ibeat = 0; ibeat < nbeat; ibeat++) {
			for (i = 0; i < hanningp; i++) {
				if (beat[ibeat].i <  beatwinpo2) continue;
				if (beat[ibeat].i >= beat[nbeat - 1].i - beatwinpo2) break;
				beat[ibeat].bkgfit[ichan][lentempp - hanningp + i] *= hanning[i];
				beat[ibeat].bkgfit[ichan][i] *= hanning[hanningp - 1 - i];
			}
		}

/************************/
/* reduce BKG with mGLM */
/************************/
		for (ibeat = 0; ibeat < nbeat; ibeat++) {
			if (beat[ibeat].i < beatwinpo2) {
				for (kbeat = 0; kbeat < nbeat; kbeat++) if (beat[kbeat].i > beatwinpo2) break;
			} else if (beat[ibeat].i >= beat[nbeat - 1].i - beatwinpo2) {
				for (kbeat = nbeat - 1; kbeat > 0; kbeat--) {
					if (beat[kbeat].i < beat[nbeat - 1].i - beatwinpo2) break;
				}
			} else {
				kbeat = ibeat;
			}
			for (ii = beat[ibeat].i, i = 0; i < lentempp; i++, ii++) {
				if (ii < 0) continue;
				if (ii >= npts) break;
				eeg[ichan][ii] -= beat[kbeat].bkgfit[ichan][i];
				bcg[ichan][ii] += beat[kbeat].bkgfit[ichan][i];
			}
		}
	}

/*********************************************************************/
/* find first fully modeled beat (excluding effect of artifact gaps) */
/*********************************************************************/
	for (ibeat0 = ibeat = 0; ibeat < nbeat; ibeat++) {
		if (beat[ibeat].i <  beatwinpo2) continue;
		ibeat0 = ibeat; break;
	}
/**********************************************************************/
/* find last+1 fully modeled beat (excluding effect of artifact gaps) */
/**********************************************************************/
	for (ibeat1 = ibeat = nbeat - 1; ibeat > 0; ibeat--) {
		if (beat[ibeat].i > beat[nbeat - 1].i - beatwinpo2) continue;
		ibeat1 = ibeat; break;
	}
}

if (aas) {
	for (ichan = 0; ichan < nchan; ichan++) {
		for (ibeat = 0; ibeat < nbeat; ibeat++) {
			if (ibeat < navgo2) continue;
			if (ibeat >= nbeat - navgo2) continue;
			AAS (ibeat, ichan);
		}

/******************************************/
/* eliminate discontinuities with hanning */
/******************************************/
		for (ibeat = 0; ibeat < nbeat; ibeat++) {
			for (i = 0; i < hanningp; i++) {
				if (ibeat < navgo2) continue;
				if (ibeat >= nbeat - navgo2) break;
				beat[ibeat].bkgfit[ichan][lentempp - hanningp + i] *= hanning[i];
				beat[ibeat].bkgfit[ichan][i] *= hanning[hanningp - 1 - i];
			}
		}

/***********************/
/* reduce BKG with AAS */
/************************/
		for (ibeat = 0; ibeat < nbeat; ibeat++) {
			kbeat = ibeat;
			if (ibeat < navgo2) kbeat = navgo2;
			if (ibeat >= nbeat - navgo2) kbeat = nbeat - navgo2 - 1;
			for (ii = beat[ibeat].i, i = 0; i < lentempp; i++, ii++) {
				if (ii < 0) continue;
				if (ii >= npts) break;
				eeg[ichan][ii] -= beat[kbeat].bkgfit[ichan][i];
				bcg[ichan][ii] += beat[kbeat].bkgfit[ichan][i];
			}
		}
	}

/****************************************************************/
/* first fully modeled beat (excluding effect of artifact gaps) */
/****************************************************************/
	ibeat0 = navgo2;
/*****************************************************************/
/* last+1 fully modeled beat (excluding effect of artifact gaps) */
/*****************************************************************/
	ibeat1 = nbeat - navgo2;
}

/**********/
/* output */
/**********/
	if (!(ptr = strrchr (datroot, '/'))) ptr = datroot; else ptr++;
	strcpy (outroot, ptr);

	if (outputN) {
/*********************/
/* write NS dat file */
/*********************/
		sprintf (outfile, "%s_%s.dat", outroot, trailer);
		writeNS_file (eeg, outfile);
	}

	if (outputE) {
/**********************************/
/* write BKG estimate NS dat file */
/**********************************/
		sprintf (outfile, "%s_%s_eBKG.dat", outroot, trailer);
		writeNS_file (bcg, outfile);
	}

for (mchan = ichan = 0; ichan < nchan; ichan++) if (channel[ichan].diag) {
	if (outputC) {
/************************************************************************************/
/* output beat timecourse average and s.d. for channels of interest after reduction */
/************************************************************************************/
		nbeat_in_avg = avgsd (1, ichan, mchan);
		sprintf (outfile, "%s_%s_%s_avg.dat", outroot, trailer, channel[ichan].trode);
		fprintf (stdout, "#Writing: %s\n", outfile);
		if (!(outfp = fopen (outfile, "w"))) errw (program, outfile);
		write_command_line (outfp, argc, argv, gbr);
		fprintf (outfp, "#average and s.d. of %d beats\n", nbeat_in_avg);
		fprintf (outfp, "#%9s%10s%10s%10s%10s\n", "latency", "mean_pre", "sd_pre", "mean_post", "sd_post");
		for (i = 0; i < avglenp; i++) {
			fprintf (outfp, "%10.4f", i*delta);
			fprintf (outfp, "%10.2f%10.2f", beatavg[0][mchan][i], sqrt (beatvar[0][mchan][i]));
			fprintf (outfp, "%10.2f%10.2f", beatavg[1][mchan][i], sqrt (beatvar[1][mchan][i]));
			fprintf (outfp, "\n");
		}
		if (fclose (outfp)) errw (program, outfile);
	}

	if (outputT) {
/**************************/
/* write timecourse files */
/**************************/
		sprintf (outfile, "%s_%s_%s.dat", outroot, trailer, channel[ichan].trode);
		fprintf (stdout, "#Writing: %s\n", outfile);
		if (!(outfp = fopen (outfile, "w"))) errw (program, outfile);
		write_command_line (outfp, argc, argv, gbr);
		fprintf (outfp, "#%9s%10s%10s%10s\n", "time(sec)", "EEG_pre", "EEG_post", "artifact");
		for (i = 0; i < npts; i++) {
			fprintf (outfp, "%10.4f%10.2f%10.2f%10d\n",
			i*delta, timecourse[0][mchan][i], timecourse[1][mchan][i], artifact[i]);
		}
		if (fclose (outfp)) errw (program, outfile);
	}

	if (outputM > 0) {
/*********************/
/* write model files */
/*********************/
		mbeat = outputM;
		kbeat = ibeat0;
		m = 0; do {
			sprintf (outfile, "%s_%s_%s_BKGmodels%d.dat", outroot, trailer, channel[ichan].trode, m);
			fprintf (stdout, "#Writing: %s\n", outfile);
			if (!(outfp = fopen (outfile, "w"))) errw (program, outfile);
			write_command_line (outfp, argc, argv, gbr);
			sprintf (command, "#beats (counting from 0)");
			for (i = 0; i < lentempp; i++) {
				fprintf (outfp, "%10.4f", i*delta);
				k = 0; for (ibeat = kbeat; ibeat < ibeat1; ibeat++) {
					if (beat[ibeat].art) continue;
					fprintf (outfp, "%10.2f", beat[ibeat].bkgfit[ichan][i]);
					if (!i) {
						mbeat--;
						l = strlen (command); sprintf (command + l, " %d", ibeat);
						if (0) printf ("k=%d mbeat=%d %s\n", k, mbeat, command);
					}
					if (++k == MAXM) break;
				}
				fprintf (outfp, "\n");
			}
			fprintf (outfp, "%s\n", command);
			if (fclose (outfp)) errw (program, outfile);
			m++;
			kbeat = ibeat + 1;
			if (0) printf ("mbeat=%d kbeat=%d\n", mbeat, kbeat);
		} while (mbeat > 0 && kbeat < ibeat1);
	}

	if (outputL > 0) {
/********************************************************************************************/
/* write mean of first and last specified number of beats and corresponding BKG-reduced EEG */
/********************************************************************************************/
		for (k = 0; k < 4; k++) for (i = 0; i < lentempp; i++) bkgmean[k][i] = 0.0;
		m = 0; for (ibeat = ibeat0; ibeat < ibeat1; ibeat++) {
			if (beat[ibeat].art) continue;
			for (i = 0; i < lentempp; i++) {
				bkgmean[0][i] += beat[ibeat].bkgfit[ichan][i];
				bkgmean[2][i] += eeg[ichan][beat[ibeat].i + i];
			}
			if (++m == outputL) break;
		}
		for (i = 0; i < lentempp; i++) {
			bkgmean[0][i] /= m;
			bkgmean[2][i] /= m;
		}
		l = 0; for (ibeat = ibeat1 - 1; ibeat > ibeat0; ibeat--) {
			if (beat[ibeat].art) continue;
			for (i = 0; i < lentempp; i++) {
				bkgmean[1][i] += beat[ibeat].bkgfit[ichan][i];
				bkgmean[3][i] += eeg[ichan][beat[ibeat].i + i];
			}
			if (++l == outputL) break;
		}
		for (i = 0; i < lentempp; i++) {
			bkgmean[1][i] /= l;
			bkgmean[3][i] /= l;
		}

		sprintf (outfile, "%s_%s_%s_BKGmean.dat", outroot, trailer, channel[ichan].trode);
		fprintf (stdout, "#Writing: %s\n", outfile);
		if (!(outfp = fopen (outfile, "w"))) errw (program, outfile);
		write_command_line (outfp, argc, argv, gbr);
		fprintf (outfp, "#number of waveforms averaged = %d (first) and %d (last)\n", m, l);
		fprintf (outfp, "#%9s%10s%10s%10s%10s\n", "time(sec)", "firstBKG", "lastBKG", "firstEEG", "lastEEG");
		for (i = 0; i < lentempp; i++) {
			fprintf (outfp, "%10.4f", i*delta);
			for (k = 0; k < 4; k++) fprintf (outfp, "%10.2f", bkgmean[k][i]);
			fprintf (outfp, "\n");
		}
		if (fclose (outfp)) errw (program, outfile);
	}

	if (outputB) {
/***************************************/
/* compute mean BKG and write psd file */
/***************************************/
		for (i = 0; i < lentempp; i++) bkgmean[0][i] = 0.0;
		k = 0; for (ibeat = ibeat0; ibeat < ibeat1; ibeat++) {
			if (beat[ibeat].art) continue;
			for (i = 0; i < lentempp; i++) bkgmean[0][i] += beat[ibeat].bkgfit[ichan][i];
			k++;
		}
		for (i = 0; i < lentempp; i++) bkgmean[0][i] /= k;
		bkgfft (bkgmean[0], bkgpsd);
		hzpbin = rate/lentempp_pad;
		sprintf (outfile, "%s_%s_%s_BKGpsd.dat", outroot, trailer, channel[ichan].trode);
		fprintf (stdout, "#Writing: %s\n", outfile);
		if (!(outfp = fopen (outfile, "w"))) errw (program, outfile);
		write_command_line (outfp, argc, argv, gbr);
		fprintf (outfp, "#%9s%10s\n", "freq (Hz)", channel[ichan].trode);
		fprintf (outfp, "%10.4f%10.4f\n", 0., 0.);
		for (i = 0; i <= lentempp_pad/2; i++) {
			fprintf (outfp, "%10.4f%10.4f\n", hzpbin*i,       bkgpsd[i]);	
			fprintf (outfp, "%10.4f%10.4f\n", hzpbin*(i + 1), bkgpsd[i]);	
		}
		fprintf (outfp, "%10.4f%10.4f\n", hzpbin*lentempp_pad/2, 0.);
		if (fclose (outfp)) errw (program, outfile);
	}

	if (outputP) {
/*******************************************/
/* FFT channel of interest after reduction */
/*******************************************/
		hzpbin = rate/fftpts;
		nepoch = eegfft (ichan, eegpsd[1][mchan]);
		sprintf (outfile, "%s_%s_%s_EEGpsd.dat", outroot, trailer, channel[ichan].trode);
		fprintf (stdout, "#Writing: %s\n", outfile);
		if (!(outfp = fopen (outfile, "w"))) errw (program, outfile);
		write_command_line (outfp, argc, argv, gbr);
		fprintf (outfp, "#fftpts = %d hz per bin = %.4f number of epochs = %d\n", fftpts, hzpbin, nepoch);
		fprintf (outfp, "#%9s%10s%10s\n", "freq (Hz)", "pre", "post");
		fprintf (outfp, "%10.4f%10.4f%10.4f\n", 0., 0., 0.);
		for (i = 0; i <= fftpts/2; i++) {
			fprintf (outfp, "%10.4f%10.4f%10.4f\n", i*hzpbin,       eegpsd[0][mchan][i], eegpsd[1][mchan][i]);
			fprintf (outfp, "%10.4f%10.4f%10.4f\n", (i + 1)*hzpbin, eegpsd[0][mchan][i], eegpsd[1][mchan][i]);
		}
		fprintf (outfp, "%10.4f%10.4f%10.4f\n", i*hzpbin, 0., 0.);
		if (fclose (outfp)) errw (program, outfile);
	}

	mchan++;
}	/* ichan loop for diagnostic output */

/***********************************/
/* write reduce_BKG diagnosis file */
/***********************************/
	sprintf (outfile, "%s_%s_diag.dat", outroot, trailer);
	fprintf (stdout, "#Writing: %s\n", outfile);
	if (!(outfp = fopen (outfile, "w"))) errw (program, outfile);
	write_command_line (outfp, argc, argv, gbr);
	fprintf (outfp, "#%9s%10s%12s%12s%10s%10s%10s\n",
		"beat", "latency", "det", "cndnum", "mstart", "mend", "artifact");
	for (ibeat = 0; ibeat < nbeat; ibeat++) {
		fprintf (outfp, "%10i%10.4f%12.4e%12.4e", ibeat, delta*beat[ibeat].i, beat[ibeat].det, beat[ibeat].cndnum); 
		fprintf (outfp, "%10i%10i%10i\n", beat[ibeat].mstart, beat[ibeat].mend, beat[ibeat].art);
	}
	if (fclose (outfp)) errw (program, outfile);

/***************/
/* free memory */
/***************/
FREE:	free_float2 (eeg);
	free_float2 (bcg);
	free_float2 (bkgmean); free (bkgpsd);
	free (artifact);
	free (hanning);
	free (zempdur); free (zempart);
	free_float3 (beatavg);
	free_float3 (beatvar);
	free_float3 (timecourse);
	free_float3 (eegpsd);
	if (debug) printf ("deallocating beat");
	for (ibeat = 0; ibeat < nbeat; ibeat++) {if (debug) printf (" %d", ibeat); fflush (stdout);
		free_float2 (beat[ibeat].bkgfit);
	}
	if (debug) printf ("\n"); fflush (stdout);
	free (beat); free (channel);
	exit (0);
}
