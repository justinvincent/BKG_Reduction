#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <float.h>
#include <string.h>
#include <ctype.h>
#include <assert.h>

#define MAXL		256	/* filename character length */
#define SDCRIT		3.0	/* spike detector criterion (deviation from mean in units of sd) */
#define NZ		4	/* number of z signals */
#define ZCRIT0		2.5	/* z[0] (undifferentiated signal) beat detector criterion in sd units (1.24% level) */
#define ZCRIT1		2.5	/* z[1] (signal curvature)        beat detector criterion in sd units (1.24% level) */
#define BPSTOL		0.22	/* bps tolerance */
#define IBIFRAC		0.24	/* IBI histogram tolerance as fraction of IBI mode */
#define IBITOL		3.5	/* IBI histogram tolerance sd deviation from mode in units of sd */
#define NDIM		3	/* dimensionality of timeseries */
#define TAU		14.0	/* lagged product lag in msec */
#define RSEC		6.0	/* epoch length in sec used in meanIBI () when exemplar beat is given */
#define LSHIFT		0.2	/* left shift of template start relative to lagged product peak in sec */

typedef struct {
	int	i;		/* time point in EKG record corresponding to start of beat */
	int	j;		/* time point in EKG record corresponding to start of exemplar beat */
	float	signal;		/* strength of signal used to detect beat */
	short	longR, shortR, longL, shortL;	/* IBI tests */
	short	fastR, slowR,  fastL, slowL;	/* bps tests */
} BEAT;

static int beatcomp (const void *i, const void *j) {
	BEAT  *pi, *pj;
	pi = (BEAT *) i; pj = (BEAT *) j;
	return (pi->signal - pj->signal > 0) ? -1 : 1;
}

void errm (char* program) {
	fprintf (stderr, "%s: memory allocation error\n", program);
	exit (-1);
}

void errr (char* program, char* filespc) {
	fprintf (stderr, "%s: %s read error\n", program, filespc);
	exit (-1);
}

void errw (char* program, char* filespc) {
	fprintf (stderr, "%s: %s write error\n", program, filespc);
	exit (-1);
}

void getroot (char *filespc, char *filroot) {
	char	*str;
	strcpy (filroot, filespc);
	while (str = strrchr (filroot, '.')) {
			if (!strcmp (str, ".rec"))	*str = '\0';
		else	if (!strcmp (str, ".dat"))	*str = '\0';
		else	break;
	}
}

int split (char *string, char *srgv[], int maxp) {
	int	i, m;
	char	*ptr;

	if (ptr = strchr (string, '#')) *ptr = '\0';
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

double zeromean (float *f, int n) {
	int		i;
	double		u;

	for (u = i = 0; i < n; i++) u += f[i];
	u /= n;
	for (i = 0; i < n; i++) f[i] -= u;
	return u;
}

double unitvar (float *f, int npts) {
	int		i;
	double		v;

	for (v = i = 0; i < npts; i++) v += f[i]*f[i];
	v /= npts - 1; v = sqrt (v);
	for (i = 0; i < npts; i++) f[i] /= v;
	return v;
}

double trendout (float *f, int n) {
	int	i;
	double	x, sxx, sy, sxy, a[2];

	sxx = (double) n*(n+1) / (3.*(n - 1));
	for (sy = sxy = i = 0; i < n; i++) {
		x = -1. + 2.*i/(n - 1);
		sy  += f[i];
		sxy += f[i]*x;
	}
	a[0] = sy/n;
	a[1] = sxy/sxx;
	for (i = 0; i < n; i++) {
		x = -1. + 2.*i/(n - 1);
		f[i] -= a[1]*x;
	}
	return a[1];
}

/********************/
/* global variables */
/********************/
char	rcsid[] = "$Id: EKG_phase.c,v 1.30 2006/12/16 23:07:10 avi Exp $";
char	program[MAXL];
int	debug = 0;

void errb () {
	fprintf (stderr, "%s: beat buffer allocation exceeded\n", program);
	exit (-1);
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

void smooth (float *f, int n, int nsmooth) {
	float	*h, *ptr1, *ptr2, *ptr;
	int	k, j;

	if (!(h = (float *) calloc (n, sizeof (float)))) errm (program);
	ptr1 = f; ptr2 = h;
	for (k = 0; k < nsmooth; k++) {
		ptr2[0]   = ptr1[0]; 
		ptr2[n-1] = ptr1[n-1];
		for (j = 1; j < n - 1; j++) {
			ptr2[j] = 0.25*ptr1[j - 1] + 0.5*ptr1[j] + 0.25*ptr1[j + 1];
		}
		ptr = ptr1; ptr1 = ptr2; ptr2 = ptr;
	}
	if (nsmooth % 2) for (j = 1; j < n - 1; j++) f[j] = h[j];
	free (h);
}

void smoothc (float *f, int n, int nsmooth) {
	float	*h, *ptr1, *ptr2, *ptr;
	int	k, j;

	if (!(h = (float *) calloc (n, sizeof (float)))) errm (program);
	ptr1 = f; ptr2 = h;
	for (k = 0; k < nsmooth; k++) {
		for (j = 0; j < n; j++) {
			ptr2[j] = 0.25*ptr1[(n + j - 1) % n] + 0.5*ptr1[j] + 0.25*ptr1[(j + 1) % n];
		}
		ptr = ptr1; ptr1 = ptr2; ptr2 = ptr;
	}
	if (nsmooth % 2) for (j = 0; j < n; j++) f[j] = h[j];
	free (h);
}

int meanIBI0 (float rate, int maxu, float *ekg[NDIM], int npts) {
	int		l, k, i, kmax;
	float	        delta;
	double		*acfsum;
	
	delta = 1./rate;
	if (!(acfsum  = (double *) calloc (maxu, sizeof (double)))) errm (program);
	for (k = 0; k < maxu; k++) {
		for (i = 0; i < npts - maxu; i++) {
			for (l = 1; l < NDIM; l++) {	/* not using ekg[0] probably is better */
				acfsum[k] += ekg[l][i]*ekg[l][i + k];
			}
		}
	}
	for (k = 0; k < maxu; k++) {
		acfsum[k] /= npts - maxu;
		if (0) printf ("%10.4f%10.4f%10.4f\n", delta*k, ekg[0][k], acfsum[k]*10);
	}
	for (kmax = k = rate/2; k < maxu; k++) {
		if (acfsum[k] > acfsum[k - 1] && acfsum[k] > acfsum[k + 1] && acfsum[k] > acfsum[kmax]) kmax = k;
	}
	free (acfsum);
	return kmax;
}

int meanIBI (float rate, int maxu, float *ekg[NDIM], int npts, int istart) {
	int		l, k, i, kmax, mpts;
	float	        delta;
	double		*acfsum;

	delta = 1./rate;
	if (!(acfsum  = (double *) calloc (maxu, sizeof (double)))) errm (program);

	mpts = (istart) ? npts : npts - maxu;
	for (k = 0; k < maxu; k++) {
		for (i = 0; i < npts; i++) {
			for (l = 1; l < NDIM; l++) {	/* not using ekg[0] probably is better */
				acfsum[k] += ekg[l][istart + i]*ekg[l][istart + i + k];
			}
		}
	}
	for (k = 0; k < maxu; k++) {
		acfsum[k] /= mpts;
		if (0) printf ("%10.4f%10.4f%10.4f\n", delta*k, ekg[0][k], acfsum[k]*10);
	}
	for (kmax = k = rate/2; k < maxu; k++) {
		if (acfsum[k] > acfsum[k - 1] && acfsum[k] > acfsum[k + 1] && acfsum[k] > acfsum[kmax]) kmax = k;
	}
	free (acfsum);
	return kmax;
}

int beat_delete (int ibeat, BEAT *beat, int nbeat0) {
	int	i;

	for (i = ibeat; i < nbeat0 - 1; i++) beat[i] = beat[i + 1];
	return nbeat0 - 1;
}

int beat_insert (BEAT beatnew, int ibefore, BEAT *beat, int nbeat0) {
	int	i;

	for (i = nbeat0; i > ibefore; i--) beat[i] = beat[i - 1];
	beat[i] = beatnew;
	return nbeat0 + 1;
}

int spikeout (float *f, int n, float sdcrit) {
	double	sum, var, sd, q, q1, q2, sign;
	int	i, nspike, m;

	if (n < 3) return -1;
	nspike = q = sum = var = 0.0;
	for (i = 1; i < (n - 1); i++) {
		q = .5*fabs (2.*f[i] - f[i + 1] - f[i - 1]);
		sum += q; var += q*q;
	}
	m = n - 2;
	var -= sum*sum/m; sum /= m; var /= (m - 1); sd = sqrt (var);
	printf ("#spike measure mean and s.d. %10.6f % 10.6f\n", sum, sd);

	for (i = 1; i < n - 1; i++) {
		sign = (2.*f[i] - f[i + 1] - f[i - 1] > 0.) ? +1. : -1.;
		q1 = sign*(f[i] - f[i + 1]);
		q2 = sign*(f[i] - f[i - 1]);
		q = sdcrit*sd + sum;
		if (q1 > q && q2 > q) {
			f[i] = 0.5*(f[i + 1] + f[i - 1]);
			nspike++;
			if (0) printf ("#spike removed at point %10d (counting from 0)\n", i);
		}
	}
	return nspike;
}

void write_command_line (FILE *outfp, int argc, char *argv[]) {
	int		i;

	fprintf (outfp, "#%s", program);
	for (i = 1; i < argc; i++) fprintf (outfp, " %s", argv[i]);
	fprintf (outfp, "\n");
	fprintf (outfp, "#%s\n", rcsid);
}

extern void butt1dba_ (float *data, int *n, float *delta, float *fhalf_lo, int *iorder_lo,
							  float *fhalf_hi, int *iorder_hi);	/* butt1d.f */
extern void eigen_ (float *A, float *W, int *ndim);						/* eigen.f */

void usage (char *program) {
	fprintf (stderr, "Usage:\t%s <EKG_datfile>\n", program);
	fprintf (stderr, "e.g.,\t%s 04_1004_ascii_EKG[.dat]\n", program);
	fprintf (stderr, "\toption\n");
	fprintf (stderr, "\t-e\tdump (expanded) EKG record to enable manual selection of exemplar beat\n");
	fprintf (stderr, "\t-q\tdisable automatic beat deletion/insertion\n");
	fprintf (stderr, "\t-r\trecompile IBIhist and bps after each beat deletion/insertion\n");
	fprintf (stderr, "\t-t<flt>\tenter manually identified exemplar beat latency in sec\n");
	fprintf (stderr, "\t-k<flt>\tspecify template start to lagged product peak interval in sec (default=%.2f)\n", LSHIFT);
	fprintf (stderr, "\t-l<flt>\tspecify lagged product interval in msec (default=%.2f)\n", TAU);
	fprintf (stderr, "\t-Y<int>\tsuppress use of specified component of y to compute z[0]\n");
	fprintf (stderr, "\t-d<flt>\tmanually delete beat within 100 msec of specified time\n");
	fprintf (stderr, "\t-i<flt>\tmanually insert beat at specified time in sec\n");
	fprintf (stderr, "\t-z0<flt>\tspecify beat detector z[0] criterion in s.d. units (default=%.2f)\n", ZCRIT0);
	fprintf (stderr, "\t-z1<flt>\tspecify beat detector z[1] criterion in s.d. units (default=%.2f)\n", ZCRIT1);
	fprintf (stderr, "\t-Z\treport discrepancy between z[0] vs. z[1] meeting criterion\n");
	fprintf (stderr, "\t-b<flt>\tspecify false beat detector bps tolerance as fraction of local bps (default=%.2f)\n",
				BPSTOL);
	fprintf (stderr, "\t-h<flt>\tspecify false beat detector IBI histogram tolerance in s.d. units (default=%.2f)\n",
				IBITOL);
	fprintf (stderr, "\t-j<flt>\tspecify beat insertion IBI tolerance as fraction of IBI mode (default=%.2f)\n",
				IBIFRAC);
	fprintf (stderr, "\t-s<flt>\tspecify spike detector criterion in s.d. units (default=%.2f)\n", SDCRIT);
	fprintf (stderr, "\t-a<flt>\tspecify number of seconds used to estimate initial IBI (default=%.2f)\n", RSEC);
	fprintf (stderr, "\t-R\trevise beat placements by reading the diag file\n");
	fprintf (stderr, "\t\tuse -R in conjunction with -i or -d options, -q is implied\n");
	fprintf (stderr, "N.B.:\toption -k has no effect if exemplar beat start time (-t) given\n");
	exit (1);
}

int main (int argc, char *argv[]) {
	char		*ptr, command[MAXL], string[MAXL], *srgv[MAXL];
	char		stars[8];
	int		c, i, j, k, l, m, n, r;
	double		q;

/************/
/* file i/o */
/************/
	FILE		*datfp, *outfp, *diagfp;
	char		datroot[MAXL], datfile[MAXL], outfile[MAXL], diagfile[MAXL];

/*******/
/* EKG */
/*******/
	int		npts, nspike;
	int		maxu, lentemp, mlentemp, ndim = NDIM;
	int		ii, itau;
	float		rate, delta, tau = TAU, lshift = LSHIFT;
	float		*ekg[NDIM];
	float		*beatavg[NDIM], *beatvar[NDIM];
	float		time0, time1;
	float		mean, sd, HR;
	BEAT		*beat, beatnew;
	int		ibeat, jbeat, nbeat, mbeat;
	int		lbeat = 100;	/* number beat used to compute beat covariance */
	float		sdcrit = SDCRIT;

/****************/
/* correlations */
/****************/
	int		uIBI;		/* mean interbeat interval returned by meanIBI () */
	int		kmax, istart;
	float		rsec = RSEC;	/* used to get initial esitmate of cardiac period */

/**********/
/* signal */
/**********/
	double		zsum[2], zvar[2], zcrit[2] = {ZCRIT0, ZCRIT1};
	double		x2sum, x2var, x2thresh, wsum;
	float		*y[NDIM], *z[NZ], *x2, **beatcov, **beatwei, **template;
	float		*accept;
	double		signal, lnS;

/**********************/
/* manual delete list */
/**********************/
	float		dtime[64];
	int		nmdel = 0, imdel;

/**********************/
/* manual insert list */
/**********************/
	float		itime[64];
	int		nmins = 0, imins;

/*****************/
/* IBI histogram */
/*****************/
	float		*g, sigma, IBItol = IBITOL, IBIfrac = IBIFRAC;
	int		IBImode, IBI, lIBI, rIBI, FWHM, nsmooth, IBIklo, IBIkhi;
	int		IBItest, IBImin, IBImax;

/********************/
/* instantaneous HR */
/********************/
	float		bpstol = BPSTOL;
	float		*bps[2], bpsleft, bpsright;
	float		bpsrate = 10.0, bpsdelta = 0.1;
	int		bpsn, bpsnsmooth = 40;
	int		bpstest, kleft, krigt;

/**********************/
/* butterworth filter */
/**********************/
	float		fhalf_lo = 2.0, fhalf_hi = 100.0;
	int		iorder_lo = 2, iorder_hi = 0;	/* "lo" refers to low freq limit (e.g., high pass) */

/*********/
/* flags */
/*********/
	int		report_z = 0;
	int		usey[NDIM] = {1, 1, 1};
	int		iter;
	int		have_accept;		/* set when EKG_datfile input has three columns */
	int		recompile1 = 0;		/* recompile IBIhist and bps after each modification */
	int		dump_y = 0;		/* dump y to a file */
	int		modify = 1;		/* enable beat deletion and insertion */
	int		revise = 0;		/* revise beat placements by reading the diag file */
	float		ttemp = 0.0;		/* time of manually identified exemplar beat */

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
				case 'e': dump_y++;			break;
				case 'q': modify = 0;			break;
				case 'r': recompile1++;			break;
				case 'Z': report_z++;			break;
				case 'R': revise = 1; modify = 0;	break;
				case 'a': rsec = atof (ptr);		*ptr = '\0'; break;
				case 'd': dtime[nmdel++] = atof (ptr);	*ptr = '\0'; break;
				case 'i': itime[nmins++] = atof (ptr);	*ptr = '\0'; break;
				case 'k': lshift = atof (ptr);		*ptr = '\0'; break;
				case 'l': tau = atof (ptr);		*ptr = '\0'; break;
				case 'h': IBItol = atof (ptr);		*ptr = '\0'; break;
				case 'j': IBIfrac = atof (ptr);		*ptr = '\0'; break;
				case 's': sdcrit = atof (ptr);		*ptr = '\0'; break;
				case 't': ttemp = atof (ptr);		*ptr = '\0'; break;
				case 'b': bpstol = atof (ptr);		*ptr = '\0'; break;
				case 'z': switch (*ptr++) {
						case '0': zcrit[0] = atof (ptr); break;
						case '1': zcrit[1] = atof (ptr); break;
					  }				*ptr = '\0'; break;
				case 'Y': l = atoi (ptr); if (l >= 0 && l <= 2) usey[l] = 0;
									*ptr = '\0'; break;
			}
		}
		else switch (k) {
			case 0:	getroot (argv[i], datroot);	k++; break;
		}	
	}
	if (k < 1) usage (program);

	sprintf (datfile, "%s.dat", datroot);
	if (!(datfp = fopen (datfile, "r"))) errr (program, datfile);
	npts = 0; have_accept = 1; while (fgets (string, MAXL, datfp)) {
		m = split (string, srgv, MAXL);
		if (!m) continue;
		if (m < 2) {
			fprintf (stderr, "%s: %s format error (number of fields in line < 2)\n", program, datfile);
			exit (-1);
		}
		if (m < 3) have_accept = 0;
		npts++;
	}
	if (!have_accept) fprintf (stderr, "WARNING: %s has no accept column\n", datfile);
	rewind (datfp);
	n = 1; while (n < npts + 1024) n *= 2;

/***********************/
/* allocate EKG memory */
/***********************/
	for (l = 0; l < ndim; l++) {
		if (!(ekg[l] = (float *) calloc (n, sizeof (float)))) errm (program);
		if (!(y[l]   = (float *) calloc (n, sizeof (float)))) errm (program);
	}
	for (l = 0; l < NZ; l++) {
		if (!(z[l]   = (float *) calloc (n, sizeof (float)))) errm (program);
	}
	if (!(accept = (float *) calloc (n, sizeof (float)))) errm (program);
	i = 0; while (fgets (string, MAXL, datfp)) {
		m = split (string, srgv, MAXL);
		if (!m) continue;
		if (!i) time0 = atof (srgv[0]);
		ekg[0][i] = atof (srgv[1]);
		accept[i] = (m >= 3) ? atof (srgv[2]) : 0;
		i++;
	}
	time1 = atof (srgv[0]);
	if (fclose (datfp)) errr (program, datfile);
	rate = (npts - 1) / (time1 - time0);
	mlentemp = 2*rate;
	delta = 1./rate;
	maxu = 1.5*rate;
	itau = 0.5 + rate*tau/1000.;
	printf ("#npts=%d n=%d rate=%.4f itau=%d\n", npts, n, rate, itau);

/***********************/
/* allocate bps memory */
/***********************/
	bpsn = npts*bpsrate/rate;
	printf ("#bpsn=%d\n", bpsn);
	for (l = 0; l < 2; l++) {
		if (!(bps[l] = (float *) malloc (bpsn * sizeof (float)))) errm (program);
	}
/*****************************************/
/* allocate beat average and s.d. arrays */
/*****************************************/
	for (l = 0; l < ndim; l++) {
		if (!(beatavg[l] = (float *) calloc (mlentemp, sizeof (float)))) errm (program);
		if (!(beatvar[l] = (float *) calloc (mlentemp, sizeof (float)))) errm (program);
	}
/******************/
/* pad and filter */
/******************/
	q = (ekg[0][0] - ekg[0][npts - 1]) / (n - npts + 1);
	for (i = 1; i <= n - npts; i++) {
		ekg[0][npts - 1 + i] = ekg[0][npts - 1] + q*i;
	}
	butt1dba_ (ekg[0], &n, &delta, &fhalf_lo, &iorder_lo, &fhalf_hi, &iorder_hi);

/************************************/
/* remove spikes in EKG time series */
/************************************/
	nspike = spikeout (ekg[0], npts, sdcrit);
	printf ("#%d spikes removed\n", nspike);

/********************************/
/* differentiate the timeseries */
/********************************/
	for (i = 0; i < npts + 1024; i++) {
		ekg[1][i] = 5.*(ekg[0][i + 1] - ekg[0][i]);
		ekg[2][i] = 0.001*ekg[1][i]*ekg[0][i + itau];
	}

/***************************/
/* zero mean unit variance */
/***************************/
	for (l = 0; l < ndim; l++) {
		mean = zeromean (ekg[l], npts + 1024); sd = unitvar (ekg[l], npts + 1024);
		printf ("#ekg[%d] mean and s.d. %10.2f%10.2f\n", l, mean, sd);
	}

/***************************/
/* estimate cardiac period */
/***************************/
	istart = rate*ttemp;
	if (istart < 0 || istart > npts - maxu) {
		fprintf (stderr, "%s: invalid exemplar beat time (%.2f)\n", program, ttemp);
		usage (program);
	}
	k = (istart) ? rate*rsec + 0.5: npts;
	uIBI = meanIBI (rate, maxu, ekg, k, istart);
	printf ("#estimated IBI = %.2f sec; heart rate = %3.1f bpm\n", uIBI*delta, (rate/uIBI)*60.);
	printf ("#estimated number of beats = %d\n", npts/uIBI);
	mbeat = 2.0*npts/uIBI;
	lentemp = .7*uIBI;
	if (!(beat = (BEAT *) calloc (mbeat, sizeof (BEAT))))  errm (program);
	if (!(g =   (float *) calloc (maxu,  sizeof (float)))) errm (program);
	template = calloc_float2 (ndim, lentemp);

	if (dump_y) {
/********/
/* dump */
/********/
		sprintf (outfile, "%s_dump.dat", datroot);
		printf ("#Writing: %s\n", outfile);
		if (!(outfp = fopen (outfile, "w"))) errw (program, outfile);
		write_command_line (outfp, argc, argv);
		for (i = 0; i < npts; i++) {
			fprintf (outfp, "%10.4f", i*delta);
			for (l = 0; l < ndim; l++) fprintf (outfp, "%10.4f", ekg[l][i]);
			fprintf (outfp, "\n");
		}
		if (fclose (outfp)) errw (program, outfile);
	}

if (ttemp <= 0.0) {	/* no exemplar beat specified */
/*****************************************************/
/* create tentative beat list by thresholding ekg[2] */
/*****************************************************/
	if (!(x2 = (float *) calloc (npts, sizeof (float)))) errm (program);
	for (i = 0; i < npts; i++) x2[i] = ekg[2][i];
	smoothc (x2, npts, 2);
	for (x2sum = x2var = i = 0; i < npts; i++) {q = x2[i]; x2sum += q; x2var += q*q;}
	x2var -= x2sum*x2sum/npts; x2sum /= npts; x2var /= npts - 1; sd = sqrt (x2var);
	x2thresh = x2sum + 2.5*sd;
	printf ("#x2 mean and s.d. %10.6f%10.6f\tthreshold = %.6f\n", x2sum, sd, x2thresh);
	nbeat = 0; for (i = 1; i < npts - 1; i++) {
		if (x2[i] > x2thresh && x2[i] > x2[i - 1] && x2[i] > x2[i + 1]) {
			ii = i - lshift*rate + 0.5; if (ii < 0) continue;
			if (nbeat && ii - beat[nbeat - 1].i < (int) (0.8*uIBI)) continue;
			if (ii + lentemp >= npts) break;
			if (nbeat >= mbeat) errb ();
			beat[nbeat].signal = x2[i];
			beat[nbeat++].i = ii;
		}
	}
	free (x2);
	qsort ((void *) beat, nbeat, sizeof (BEAT), beatcomp);
	printf ("#tentative nbeat = %d\n", nbeat);
	for (ibeat = 0; ibeat < nbeat; ibeat++) {
		if (debug) printf ("%10d%10.4f%10.4f\n", ibeat, delta*beat[ibeat].i, beat[ibeat].signal);
	}

	if (nbeat < lbeat) lbeat = nbeat;
/***************************************************************************/
/* compute covariance of top lbeat beats sorted by smoothed lagged product */
/***************************************************************************/
	beatcov = calloc_float2 (lbeat, lbeat);
	beatwei = calloc_float2 (lbeat, lbeat);
	for (ibeat = 0; ibeat < lbeat; ibeat++) {
	for (jbeat = ibeat; jbeat < lbeat; jbeat++) {
		for (q = k = 0; k < lentemp; k++) for (l = 0; l < ndim; l++) {
			q += ekg[l][beat[ibeat].i + k]*ekg[l][beat[jbeat].i + k];
		}
		beatcov[ibeat][jbeat] = beatcov[jbeat][ibeat] = q/lentemp;
	}}
	if (debug) {
		k = (lbeat < 10) ? lbeat : 10;
		printf ("#first %d cols of beatcov\n", k);
		for (ibeat = 0; ibeat < lbeat; ibeat++) {
			for (jbeat = 0; jbeat < k; jbeat++) printf ("%10.4f", beatcov[ibeat][jbeat]);
			printf ("\n");
		}
	}
	printf ("#beatcov eigenvalues\n");
	eigen_ (beatcov[0], beatwei[0], &lbeat);
	if (debug) for (ibeat = 0; ibeat < lbeat; ibeat++) {
		printf ("%10d%12.4e\n", ibeat, beatcov[ibeat][ibeat]);
	}
	printf ("#weights of principal beatcov eigenvector\n");
	if (debug) for (ibeat = 0; ibeat < lbeat; ibeat++) {
		printf ("%10d%10.4f\n", ibeat, beatwei[0][ibeat]);
	}
/***********************************************************/
/* create template as weighted sum using first eigenvector */
/***********************************************************/
	j = 0;		/* signals no exemplar beat */
	for (wsum = ibeat = 0; ibeat < lbeat; ibeat++) {
		for (l = 0; l < ndim; l++) for (i = 0; i < lentemp; i++) {
			template[l][i] += beatwei[0][ibeat]*ekg[l][beat[ibeat].i + i];
		}
		wsum += beatwei[0][ibeat];
	}
	for (l = 0; l < ndim; l++) for (i = 0; i < lentemp; i++) template[l][i] /= wsum;
	free_float2 (beatcov); free_float2 (beatwei);
} else {	/* (ttemp > 0.0) */
/*************************************************/
/* create template using specified exemplar beat */
/*************************************************/
	j = ttemp*rate + 0.5;
	for (l = 0; l < ndim; l++) for (i = 0; i < lentemp; i++) {
		template[l][i] = ekg[l][j + i];
	}
}

/******************/
/* write template */
/******************/
	sprintf (outfile, "%s_template.dat", datroot);
	printf ("#Writing: %s\n", outfile);
	if (!(outfp = fopen (outfile, "w"))) errw (program, outfile);
	write_command_line (outfp, argc, argv);
	for (i = 0; i < lentemp; i++) {
		fprintf (outfp, "%10.4f", i*delta);
		for (l = 0; l < ndim; l++) fprintf (outfp, "%10.4f", template[l][i]);
		fprintf (outfp, "\n");
	}
	if (fclose (outfp)) errw (program, outfile);

/*******************************/
/* compute y and z timecourses */
/*******************************/
	printf ("#computing z[0]; usey[] = {"); for (l = 0; l < ndim; l++) printf (" %d", usey[l]); printf ("}\n");
	for (i = 0; i < npts; i++) {
		for (q = l = 0; l < ndim; l++) {
			for (y[l][i] = k = 0; k < lentemp; k++) y[l][i] += template[l][k]*ekg[l][i + k];
			y[l][i] /= lentemp;
			if (usey[l]) q += y[l][i];
		}
		z[0][i] = q;
	}
	for (i = 1; i < npts - 1; i++) z[1][i] = 20.0*(2.*z[0][i] - z[0][i - 1] - z[0][i + 1]);

/**************************************/
/* locate beats by finding peaks in z */
/**************************************/
	for (l = 0; l < 2; l++) {
		zsum[l] = zvar[l] = 0.0;
		for (i = 0; i < npts; i++) {zsum[l] += z[l][i]; zvar[l] +=z[l][i]*z[l][i];}
		zvar[l] -= zsum[l]*zsum[l]/npts; zsum[l] /= npts; zvar[l] /= npts - 1; sd = sqrt (zvar[l]);
		printf ("#z[%d] mean and s.d. %10.6f%10.6f\t threshold = %.6f\n", l, zsum[l], sd, zsum[l] + zcrit[l]*sd);
	}
	nbeat = 0; for (i = 1; i < npts; i++) {
		for (z[2][i] = l = 0; l < 2; l++) {
			q = z[l][i] - zsum[l];
			if (q > 0.0) z[2][i] += q*q/zvar[l];
		}
		z[2][i] = sqrt (z[2][i]);	/* z[2] is essentially a chi */
		if (z[0][i] < z[0][i - 1] || z[0][i] < z[0][i + 1]) continue;
		m = (z[0][i] - zsum[0] > zcrit[0]*sqrt(zvar[0]));
		if  (z[1][i] - zsum[1] > zcrit[1]*sqrt(zvar[1])) m += 2;
		if (report_z && m == 1) {
			printf ("#z[0] pos z[1] neg at %10.4f sec\n", i*delta);
		}
		if (report_z && m == 2) {
			printf ("#z[1] pos z[0] neg at %10.4f sec\n", i*delta);
		}
		if (m < 3) continue;
		if (!revise) {
			if (nbeat >= mbeat) errb ();
			beat[nbeat].i = i;
			beat[nbeat].j = j;
			beat[nbeat].signal = z[0][i];
			nbeat++;
		}
	}

/***************************************************/
/* revise beat placements by reading the diag file */
/***************************************************/
	if (revise) {
		nbeat = ibeat = 0;
		sprintf (diagfile, "%s_diag.txt", datroot);
		if (!(diagfp = fopen (diagfile, "r"))) errr (program, diagfile);
		while (fgets (string, MAXL, diagfp)) {
			m = split (string, srgv, 4*MAXL);
			if (!m) continue;
			ibeat = atoi(srgv[0]) - 1;
			beat[ibeat].i = atof(srgv[1])/delta + 0.01;
			beat[ibeat].signal = atof(srgv[3]);
			if (m == 13) {
				for (k = i = 0; i < strlen(srgv[12]); i++) if (srgv[12][i] == '*') k--;
				beat[ibeat].j = k;
			} else {
				beat[ibeat].j = ttemp*rate + 0.5;
			}
			nbeat++;
		}
		if (fclose (diagfp)) errr (program, diagfile);
	}

/**********************/
/* manual beat insert */
/**********************/
	memset (&beatnew, '\0', sizeof (BEAT));
	for (imins = 0; imins < nmins; imins++) {
		beatnew.i = itime[imins]*rate;
		beatnew.j = -3;
		for (ibeat = 0; ibeat < nbeat; ibeat++) if (beat[ibeat].i > beatnew.i) break;
		printf ("#manual insert beat at %8.2f sec", itime[imins]);
		nbeat = beat_insert (beatnew, ibeat, beat, nbeat);
		printf ("\tnbeat = %d\n", nbeat);
	}

/**********************/
/* manual beat delete */
/**********************/
	for (imdel = 0; imdel < nmdel; imdel++) {
		for (ibeat = 0; ibeat < nbeat; ibeat++) {
			if (fabs (beat[ibeat].i*delta - dtime[imdel]) < 0.1) {
				printf ("#manual delete beat at %8.2f sec z[0] = %.4f\n",
					beat[ibeat].i*delta, beat[ibeat].signal);
				nbeat = beat_delete (ibeat, beat, nbeat);
			}
		}
	}

	do {	iter = 0;
/*************************/
/* compile IBI histogram */
/*************************/
		for (i = 0; i < maxu; i++) g[i] = 0.0;
		IBImin = maxu; IBImax = 0;
		for (ibeat = 1; ibeat < nbeat; ibeat++) {
			IBI = beat[ibeat].i - beat[ibeat - 1].i;
			if (IBI >= 0 && IBI < maxu) g[IBI] += 1;
			if (IBI > IBImax) IBImax = IBI; if (IBI < IBImin) IBImin = IBI;
		}
		nsmooth = 30*rate/nbeat;
		smooth (g, maxu, nsmooth);
		for (IBImode = k = rate/2; k < maxu; k++) {
			if (g[k] > g[IBImode]) IBImode = k;
		}
		IBIklo = IBImode; while (g[IBIklo] > g[IBImode]/2) IBIklo--;
		IBIkhi = IBImode; while (g[IBIkhi] > g[IBImode]/2) IBIkhi++; IBIkhi--;
		FWHM = IBIkhi - IBIklo;
		sigma = FWHM/(2.*sqrt(2.*log(2.)));
		HR = (rate/IBImode)*60.;
		printf ("#mode IBI = %.4f +/- %.4f sec heart rate = %3.1f bpm\n", IBImode*delta, sigma*delta, HR);
		if (IBImode > 1.2*uIBI || IBImode < .8*uIBI) {
			fprintf (stderr, "WARNING: IBImode differs from initial IBI estimate by more than 20 percent\n");
		}

/*****************************/
/* compute instantanous rate */
/*****************************/
		for (ibeat = 1; ibeat < nbeat; ibeat++) {
			q = rate/(beat[ibeat].i - beat[ibeat - 1].i);
			kleft = (ibeat == 1) ? 0 : beat[ibeat - 1].i*delta*bpsrate;
			krigt = (ibeat == nbeat - 1) ? bpsn : beat[ibeat].i*delta*bpsrate;
			for (k = kleft; k < krigt; k++) bps[0][k] = bps[1][k] = q;
		}
		smooth (bps[1], bpsn, bpsnsmooth);

/*****************************/
/* compile IBI and bps tests */
/*****************************/
		for (ibeat = 0; ibeat < nbeat; ibeat++) {
			l = (ibeat == nbeat - 1) ? npts : beat[ibeat + 1].i;
			IBI = l - beat[ibeat].i;
			beat[ibeat].longR = (IBI - IBImode > IBItol*sigma);
			bpsright = rate/IBI;
			krigt = (beat[ibeat].i*delta + 2.0)*bpsrate; if (krigt >= bpsn) krigt = bpsn - 1;
			beat[ibeat].slowR = (bpsright < (1.0 - bpstol)*bps[1][krigt]);
			if (ibeat < nbeat - 1) {
				beat[ibeat].shortR = (IBI - IBImode < -IBItol*sigma);
				beat[ibeat].fastR  = (bpsright > (1.0 + bpstol)*bps[1][krigt]);
			}
			l = (ibeat) ? beat[ibeat - 1].i : 0;
			IBI = beat[ibeat].i - l;
			beat[ibeat].longL = (IBI - IBImode > IBItol*sigma);
			bpsleft = rate/IBI;
			kleft = (beat[ibeat].i*delta - 2.0)*bpsrate; if (kleft < 0) kleft = 0;
			beat[ibeat].slowL = (bpsleft < (1.0 - bpstol)*bps[1][kleft]);
			if (ibeat > 0) {
				beat[ibeat].shortL = beat[ibeat - 1].shortR;
				beat[ibeat].fastL  = (bpsleft > (1.0 + bpstol)*bps[1][kleft]);
			}
		}

/**************************/
/* delete too close beats */
/**************************/
		if (modify) for (ibeat = 0; ibeat < nbeat - 1; ibeat++) {
			IBItest = beat[ibeat].shortL && beat[ibeat].shortR;
			bpstest = beat[ibeat].fastL  && beat[ibeat].fastR;
			if (IBItest && bpstest && beat[ibeat].j >= 0) {
				printf ("#delete beat at %8.2f sec z[0] = %8.4f", beat[ibeat].i*delta, beat[ibeat].signal);
				nbeat = beat_delete (ibeat, beat, nbeat);
				printf ("  nbeat = %d\n", nbeat);
				printf ("#\tbeats before and after at %.2f and %.2f sec\n",
					beat[ibeat - 1].i*delta, beat[ibeat].i*delta);
				iter++; if (recompile1) goto ONE;
			}
		}

/**********************************************************************************/
/* delete beats with too long IBI on one side and too short IBI on the other side */
/**********************************************************************************/
		if (modify) for (ibeat = 0; ibeat < nbeat; ibeat++) {
			IBItest = (beat[ibeat].longL && beat[ibeat].shortR) || (beat[ibeat].shortL && beat[ibeat].longR);
			bpstest = (beat[ibeat].slowL && beat[ibeat].fastR)  || (beat[ibeat].fastL  && beat[ibeat].slowR); 
			if (IBItest && bpstest && beat[ibeat].j >= 0) {
				printf ("#delete beat at %8.2f sec z[0] = %8.4f", beat[ibeat].i*delta, beat[ibeat].signal);
				nbeat = beat_delete (ibeat, beat, nbeat);
				printf ("  nbeat = %d\n", nbeat);
				printf ("#\tbeats before and after at %.2f and %.2f sec\n",
					beat[ibeat - 1].i*delta, beat[ibeat].i*delta);
				iter++; if (recompile1) goto ONE;
			}
		}

/******************************************************/
/* delete beats that do not satisfy the above         */
/* criteria but nevertheless are in extreme proximity */
/******************************************************/
		if (modify) for (ibeat = 0; ibeat < nbeat - 1; ibeat++) {
			if ((beat[ibeat + 1].i - beat[ibeat].i) < IBImode/3) {
/* printf ("#beats %d and %d too close at %.4f\n", ibeat, ibeat+1, beat[ibeat].i*delta); */
				m = -1;
				if (z[2][beat[ibeat].i] < z[2][beat[ibeat + 1].i] && beat[ibeat].j >= 0) {
					m = ibeat;
				} else if (z[2][beat[ibeat + 1].i] < z[2][beat[ibeat].i] && beat[ibeat + 1].j >= 0) {
					m = ibeat + 1;
				}
				if (m > -1) {
					printf ("#delete beat at %8.2f sec z[0] = %8.4f", beat[m].i*delta, beat[m].signal);
					nbeat = beat_delete (m, beat, nbeat);
					printf ("  nbeat = %d\n", nbeat);
					printf ("#\tbeats before and after at %.2f and %.2f sec\n",
						beat[m - 1].i*delta, beat[m].i*delta);
					iter++; if (recompile1) goto ONE;
				}
			}
		}

/************************/
/* fill in missed beats */
/************************/
 		if (modify) for (ibeat = 0; ibeat < nbeat; ibeat++) {
			memset (&beatnew, '\0', sizeof (BEAT));
 			IBI = (ibeat) ? beat[ibeat].i - beat[ibeat - 1].i : beat[0].i;
			if (IBI > IBImode*(1. + IBIfrac)) {
				ii = beat[ibeat].i;
/***********************************************************************************/
/* look for z[0] and z[1] signal peaks to left of beat having at least a 32% level */
/***********************************************************************************/
				m = 0;
				for (kmax = k = IBImode*(1. - IBIfrac); k < IBImode*(1. + IBIfrac); k++) {
					if (z[0][ii - k] - zsum[0] > sqrt(zvar[0])
					&&  z[1][ii - k] - zsum[1] > sqrt(zvar[1])
 					&&  z[2][ii - k] > z[2][ii - kmax]) {kmax = k; m++;}
				}
 				beatnew.signal = z[0][ii - kmax];
 				if (m) {
					beatnew.j = -1;
					goto FOUND;
				}
/********************************************************************/
/* failing to find a 32% level signal locate new beat using IBIhist */
/********************************************************************/
				for (kmax = k = IBImode*(1. - IBIfrac); k < IBImode*(1. + IBIfrac); k++) {
 					if (g[k] > g[kmax]) kmax = k;
				}
 				beatnew.signal = z[0][ii - kmax];
				beatnew.j = -2;
FOUND: 				beatnew.i = ii - kmax;
/********************************************************************/
/* ensure that inserted beat would not create too short an interval */
/********************************************************************/
				l = (ibeat) ? beat[ibeat - 1].i : 0; r = beat[ibeat].i;
				lIBI = beatnew.i - l; rIBI = r - beatnew.i;
				bpsleft = rate/lIBI; bpsright = rate/rIBI;
				kleft = (l*delta - 0.1)*bpsrate; if (kleft <     0) kleft = 0;
				krigt = (r*delta + 0.1)*bpsrate; if (krigt >= bpsn) krigt = bpsn - 1;
				beatnew.shortL = (ibeat && lIBI < IBImode*(1. - IBIfrac));
				beatnew.fastL  = (ibeat && bpsleft > (1.0 + bpstol)*bps[0][kleft]);
				beatnew.shortR = (ibeat && rIBI < IBImode*(1. - IBIfrac));
				beatnew.fastR  = (ibeat && bpsright > (1.0 + bpstol)*bps[0][krigt]);
				if (beatnew.shortL && beatnew.fastL) {
					printf ("#beat insertion at %8.2f sec prevented by ", beatnew.i*delta);
					printf ("existing beat to the left at %8.2f\n", beat[ibeat - 1].i*delta);
				} else if (beatnew.shortR && beatnew.fastR) {
					printf ("#beat insertion at %8.2f sec prevented by ", beatnew.i*delta);
					printf ("existing beat to the right at %8.2f\n", beat[ibeat].i*delta);
				} else {
					if (nbeat == mbeat) errb ();
					nbeat = beat_insert (beatnew, ibeat, beat, nbeat);
					printf ("#insert beat at %8.2f sec z[0] = %8.4f  nbeat = %d\n",
						beatnew.i*delta, beatnew.signal, nbeat);
					iter++; if (recompile1) goto ONE;
				}
			}
		}
ONE:	;
	} while (iter);
	printf ("#final number of beats = %d\n", nbeat);

/****************************************************/
/* compute y S/N using unamibguously detected beats */
/****************************************************/
	for (lnS = m = ibeat = 0; ibeat < nbeat; ibeat++) {
		if (beat[ibeat].longR || beat[ibeat].slowR || beat[ibeat].shortR || beat[ibeat].fastR
		||  beat[ibeat].longL || beat[ibeat].slowL || beat[ibeat].shortL || beat[ibeat].fastL
		|| beat[ibeat].j < -1) continue;
		lnS += log (beat[ibeat].signal);
		m++;
	}
	lnS /= m; signal = exp (lnS);

/*******************/
/* write diag file */
/*******************/
	sprintf (outfile, "%s_diag.txt", datroot);
	if (!(outfp = fopen (outfile, "w"))) errw (program, outfile);
	printf ("#Writing: %s\n", outfile);
	write_command_line (outfp, argc, argv);
	fprintf (outfp, "#y[] mask = {");
	for (l = 0; l < ndim; l++) fprintf (outfp, " %d", usey[l]); fprintf (outfp, "}\n");
	sd = sqrt (zvar[0]);
	fprintf (outfp, "#z[0] mean and s.d. %10.6f%10.6f  threshold = %.4f  2nd threshold = %.4f\n",
			zsum[0], sd, zsum[0] + zcrit[0]*sd, zsum[0] + sd);
	fprintf (outfp, "#z[0] S/N (dB) = %.4f\n", 10.*M_LOG10E*(2.0*lnS - log (zvar[0])));
	sd = sqrt (zvar[1]);
	fprintf (outfp, "#z[1] mean and s.d. %10.6f%10.6f  threshold = %.4f  2nd threshold = %.4f\n",
			zsum[1], sd, zsum[1] + zcrit[1]*sd, zsum[1] + sd);
	fprintf (outfp, "#npts=%d n=%d rate=%.4f\n", npts, n, rate);
	fprintf (outfp, "#initial IBI estimate = %.2f sec; heart rate = %3.1f bpm\n", uIBI*delta, (rate/uIBI)*60.);
	fprintf (outfp, "#estimated number of beats = %d  final number = %d\n", npts/uIBI, nbeat);
	fprintf (outfp, "#mode IBI = %.4f +/- %.4f sec\n", IBImode*delta, sigma*delta);
	fprintf (outfp, "#heart rate by IBImode = %3.1f bpm\n", HR);
	fprintf (outfp, "#ibeat     time     LIBI   signal longR slowR shortR fastR longL slowL shortL fastL insrt\n");
	for (ibeat = 0; ibeat < nbeat; ibeat++) {
		fprintf (outfp, "%6d%9.3f", ibeat + 1, beat[ibeat].i*delta);
		fprintf (outfp, "%9.3f", ((ibeat) ? beat[ibeat].i - beat[ibeat - 1].i : 0)*delta);
		fprintf (outfp, "%9.4f", beat[ibeat].signal);
		fprintf (outfp, "%6d%6d%6d%6d%6d%6d%6d%6d",
			beat[ibeat].longR, beat[ibeat].slowR, beat[ibeat].shortR, beat[ibeat].fastR,
			beat[ibeat].longL, beat[ibeat].slowL, beat[ibeat].shortL, beat[ibeat].fastL);
		strcpy (stars, ""); for (k = 0; k < -beat[ibeat].j; k++) strcat (stars, "*");
		fprintf (outfp, " %-5s\n", stars);
	}
	if (fclose (outfp)) errw (program, outfile);

/***********************/
/* write IBI histogram */
/***********************/
	sprintf (outfile, "%s_IBIhist.dat", datroot);
	if (!(outfp = fopen (outfile, "w"))) errw (program, outfile);
	printf ("#Writing: %s\n", outfile);
	write_command_line (outfp, argc, argv);
	fprintf (outfp, "#IBI min  = %10.4f sec\tIBI max = %10.4f sec\n", IBImin*delta, IBImax*delta);
	fprintf (outfp, "#IBI mode = %10.4f sec\t(histogram value=%.4f)\n", IBImode*delta, g[IBImode]);
	fprintf (outfp, "#IBI histogram nmooth = %d\n", nsmooth);
	fprintf (outfp, "#FWHM = %.4f sec (%.4fto%.4f sec)\n", FWHM*delta, IBIklo*delta, IBIkhi*delta);
	for (j = 0; j < maxu; j++) fprintf (outfp, "%10.2f%10.4f\n", j*delta, g[j]);
	if (fclose (outfp)) errw (program, outfile);

/*************************************/
/* write instantaneous HR timecourse */
/*************************************/
	sprintf (outfile, "%s_bps.dat", datroot);
	if (!(outfp = fopen (outfile, "w"))) errw (program, outfile);
	printf ("#Writing: %s\n", outfile);
	write_command_line (outfp, argc, argv);
	fprintf (outfp, "#bpsnsmooth = %d\n", bpsnsmooth);
	for (k = 0; k < bpsn; k++) fprintf (outfp, "%10.2f%10.4f%10.4f\n", k*bpsdelta, bps[0][k], bps[1][k]);
	if (fclose (outfp)) errw (program, outfile);

/******************************************/
/* compile beat timecourse average and sd */
/******************************************/
	ii = 0.3*uIBI;
	for (i = 0; i < uIBI; i++) for (l = 0; l < ndim; l++) beatavg[l][i] = beatvar[l][i] = 0.0;
	for (n = ibeat = 0; ibeat < nbeat; ibeat++) {
		if (beat[ibeat].i - ii < 0 || beat[ibeat].j < 0) continue; /* exclude inserted beats */
		for (i = 0; i < uIBI; i++) for (l = 0; l < ndim; l++) {
			q = ekg[l][beat[ibeat].i - ii + i];
			beatavg[l][i] += q;
			beatvar[l][i] += q*q;
		}
		n++;
	}
	for (i = 0; i < uIBI; i++) for (l = 0; l < ndim; l++) {
		beatvar[l][i] -= beatavg[l][i]*beatavg[l][i]/n;
		beatavg[l][i] /= n; beatvar[l][i] /= n - 1;
	}
	sprintf (outfile, "%s_beatavg.dat", datroot);
	printf ("#Writing: %s\n", outfile);
	if (!(outfp = fopen (outfile, "w"))) errw (program, outfile);
	write_command_line (outfp, argc, argv);
	fprintf (outfp, "#average and s.d. of %d beats\n", n);
	for (i = 0; i < uIBI; i++) {
		fprintf (outfp, "%10.4f", (i - ii)*delta);
		for (l = 0; l < ndim; l++) fprintf (outfp, "%10.4f", beatavg[l][i]);
		for (l = 0; l < ndim; l++) fprintf (outfp, "%10.4f", sqrt (beatvar[l][i]));
		fprintf (outfp, "\n");
	}
	if (fclose (outfp)) errw (program, outfile);

/******************************/
/* output y and z timecourses */
/******************************/
	for (i = 0; i < npts; i++) z[3][i] = 0.0;
	for (ibeat = 0; ibeat < nbeat; ibeat++) z[3][beat[ibeat].i] = -2.0;
	sprintf (outfile, "%s_yz.dat", datroot);
	printf ("#Writing: %s\n", outfile);
	if (!(outfp = fopen (outfile, "w"))) errw (program, outfile);
	write_command_line (outfp, argc, argv);
	fprintf (outfp, "#%9s%10s%10s%10s%10s%10s%10s%10s%10s\n",
		"time", "y[0]", "y[1]", "y[2]", "z[0]", "z[1]", "chi", "beat", "accept");
	for (i = 0; i < npts; i++) {
		fprintf (outfp, "%10.4f", i*delta);
		for (l = 0; l < ndim; l++) fprintf (outfp, "%10.4f", y[l][i]);
		for (l = 0; l < NZ;   l++) fprintf (outfp, "%10.4f", z[l][i]);
		fprintf (outfp, "%10.4f", accept[i]);
		fprintf (outfp, "\n");
	}
	if (fclose (outfp)) errw (program, outfile);

/*********************/
/* output tcl script */
/*********************/
	sprintf (outfile, "%s_ISE.tcl", datroot);
	printf ("#Writing: %s\n", outfile);
	if (!(outfp = fopen (outfile, "w"))) errw (program, outfile);
	write_command_line (outfp, argc, argv);
	fprintf (outfp, "#nbeat = %d\n", nbeat);
	fprintf (outfp, "#IBImode (msec) = %.2f\n", IBImode*delta*1000.);
	for (ibeat = 0; ibeat < nbeat; ibeat++) {
		fprintf (outfp, "INSERTSTIMEVENT %d 3 \"\" \"\" NORESPONSE\n", beat[ibeat].i);
	}
	if (fclose (outfp)) errw (program, outfile);

/***************/
/* free memory */
/***************/
FREE:	free (beat); free (g);
	free_float2 (template);
	for (l = 0; l < 2; l++) free (bps[l]);
	for (l = 0; l < ndim; l++) {
		free (ekg[l]);
		free (y[l]);
		free (beatavg[l]); free (beatvar[l]);
	}
	for (l = 0; l < NZ; l++) free (z[l]);
	exit (0);
}
