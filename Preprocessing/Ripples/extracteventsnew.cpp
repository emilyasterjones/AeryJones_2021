/*
    extractevents.cpp 

    mex function for extracting events where a continuous signal goes above and below a given threshold value.  

*/

#include "mex.h"
#include <unistd.h>
#include <string.h>
#include <math.h>


#define NEVENT_MEASURES    10     	/* ten fields per event:
				   start_index
				   end_index
				   peak_value
				   peak_index
				   total_area
				   mid_area
				   total_energy
				   mid_energy 
				   max_thresh */
#define MAXEVENTS    100000     /* maximum 100000 events */

void calcmeasures(double *eptr, double *startptr, double *endptr, double *dataptr, int mindur, double thresh, double* firstaboveptr);
void Usage();

/******************************************************************************
  INTERFACE FUNCTION
  */

void mexFunction(
    int           nlhs,           /* number of expected outputs */
    mxArray       *plhs[],        /* array of pointers to output arguments */
    int           nrhs,           /* number of inputs */
    const mxArray *prhs[])         /* array of pointers to input arguments */
{
    double     	*data;
    double	*dataptr, *dtmp, *dataendptr;
    double	*lastbaseptr;
    double	*nextbaseptr;
    double	*firstaboveptr;
    double	*lastaboveptr;
    char        tmpstring[200];
    double      *tmpdata;
    double	event[MAXEVENTS * NEVENT_MEASURES];
    double	*eptr;
    double 	thresh;
    double 	baseline;
    double	tmpd;
    int		minsep;
    int		mindur;
    bool	allowstartend;
    bool	inevent;
    bool	validevent;
    int 	datalen;
    int 	nevents;
    int         dims[2]= {1,1};
    int         i, j, k;
    

    /* Check numbers of arguments */
    if ((nrhs != 6) || (nlhs > 1)) {
	Usage();
    }

    /* get a pointer to the data */
    data = mxGetPr(prhs[0]);
    
    datalen = mxGetN(prhs[0]) * mxGetM(prhs[0]);
    dataendptr = data + datalen;

    /* get the values for the other inputs */
    thresh = *mxGetPr(prhs[1]);
    baseline = *mxGetPr(prhs[2]);
    minsep = round(*mxGetPr(prhs[3]));
    mindur = round(*mxGetPr(prhs[4]));
    allowstartend = round(*mxGetPr(prhs[5]));

    dataptr = data;
    lastbaseptr = data;
    /* if we don't want events that started over baseline at the beginning we
     * advanced the pointers past the suprathreshold event */
    if (!allowstartend) {
	while ((*lastbaseptr > baseline) && (lastbaseptr <= dataendptr)) {
	    lastbaseptr++;
	    dataptr++;
	}
	if (dataptr == dataendptr) {
	    mexErrMsgTxt("Error, no below baseline values in data");
	}
    }
    else {
	/* to account for the possibility that we have an event that begins at
	 * the start of the data we need to decrement lastbaseptr by 1.  This
	 * decrement is compensated for  when we call calcmeasures below.
	 * Note that this memory address should NOT be accessed */
	lastbaseptr--;
    }
    
    /* go through the data and look for threshold crossings */
    nevents = 0;
    eptr = event;
    while (dataptr <= dataendptr) {
	/* check to see if we are below baseline, and if so, set the
	 * lastbaseptr */
	if (*dataptr <= baseline) {
	    lastbaseptr = dataptr;
	}
	/* at this point in the loop we have not yet found an event, so we
	 * check the data for a threshold crossing */
	if (*dataptr >= thresh) {
	    firstaboveptr = dataptr;
	    /* this is an event that starts at lastbaseoffset */
	    /* find the time when the signal goes below threshold and then the
	     * time when it returns to baseline. */
	    dtmp = dataptr;
	    inevent = 1;
	    validevent = 1;
	    while (inevent) {
		while ((*dtmp >= thresh) && (dtmp <= dataendptr)) {
		    dtmp++;
		}
		/* the next point is the first point below threshold, so we
		 * check the length of this threshold crossing*/
		if ((dtmp - firstaboveptr) < mindur) {
		    /* this was a short threshold crossing, so we advance the
		     * data pointer the current value and continue as though no
		     * threshold crossing was detected */
		    dataptr = dtmp;
		    validevent = 0;
		    break;
		}
		else {
		    lastaboveptr = dtmp;
		}

		/* continue until we go below baseline */
		while ((*dtmp > baseline) && (dtmp <= dataendptr)) {
		    dtmp++;
		}
		nextbaseptr = dtmp - 1;
		/* we have found the end of this event if we are either at the
		 * end of the file or if we can go minsep samples forward
		 * without getting another threshold crossing */
		if (dtmp + minsep <= dataendptr) {
		    inevent = 0;
		    for (i = 0; ((i < minsep) && (dtmp <= dataendptr)); 
			    i++, dtmp++) {
			if (*dtmp > baseline) {
			    inevent = 1;
			    break;
			}
		    }
		}
		if (dtmp == dataendptr) {
		    /* we're at the end of the data so we're done */
		    inevent = 0;
		}
	    }
	    /* calculate the various measures for this event if it is valid*/
	    if (validevent) {
		if ((dtmp != dataendptr) || allowstartend) {
		    //mexPrintf("event end at %d, dataend %d, duration %d\n", 
		//	    dtmp - data, dataendptr == dtmp, 
		//	    lastaboveptr - firstaboveptr);
		    calcmeasures(eptr, ++lastbaseptr, nextbaseptr, data, 
			    	mindur, thresh, firstaboveptr);
		    eptr += NEVENT_MEASURES;
		    nevents++;
		}
		/* if we had a valid event, move dataptr and lastbaseptr to the 
		 * current sample */
		lastbaseptr = nextbaseptr;
		dataptr = nextbaseptr;
	    }
	}
	dataptr++;
    }

    plhs[0] = mxCreateDoubleMatrix(NEVENT_MEASURES, nevents, mxREAL);
    tmpdata = mxGetPr(plhs[0]);
    /* copy the events to the output structure */
    memcpy(tmpdata, event, nevents * NEVENT_MEASURES * sizeof(double));
}

void calcmeasures(double *eptr, double *startptr, double *endptr, double *dataptr, int mindur, double thresh, double *firstaboveptr)
{
    /* calculate the eight measures for each event: start_index end_index 
     * peak_value peak_index total_area mid_area total_energy mid_energy
     * max_thresh*/
    double *tmpptr, *tmpptr2;
    double *peakptr;
    double peak;
    int	   peakind;
    double area;
    double halfarea;
    double areatmp;
    int midarea;
    double energy;
    double halfenergy;
    double energytmp;
    double maxthresh;
    int midenergy;



    peak = -1e10;
    area = 0;
    energy = 0;
    tmpptr = startptr;
    /* calculate peak, total area, and energy */
    while (tmpptr <= endptr) {
	/* peak */
	if (*tmpptr > peak) {
	    peak = *tmpptr;
	    peakind = tmpptr - dataptr;
	    peakptr = tmpptr;
	}
	/* area */
	area += *tmpptr;
	/* energy */
	energy += *tmpptr * *tmpptr;
	tmpptr++;
    }

    /* work down from the peak until the two points are separated by mindur */
    tmpptr2 = peakptr - 1;
    tmpptr = peakptr + 1;
    while ((tmpptr - tmpptr2 < mindur) && (*tmpptr >= thresh) && 
	    (*tmpptr2 >= thresh)) {
	/* the maximum threshold for this event is the smaller of the two 
	 * values.  We compute this at each step so that when we step out of
	 * the loop we are done */
	maxthresh = (*tmpptr < *tmpptr2) ? *tmpptr : *tmpptr2;
	// increment or decrement to which ever point has the largest value 
	if (*(tmpptr+1) > *(tmpptr2-1)) {
	    tmpptr++;
	}
	else {
	    tmpptr2--;
	}
    }
    
    tmpptr = startptr;
    halfarea = area / 2.0;
    midarea = -1;
    halfenergy = energy / 2.0;
    midenergy = -1;
    areatmp = 0;
    energytmp = 0;
    /* find the indeces for the middle of the area and the energy */
    while (tmpptr <= endptr) {
	areatmp += *tmpptr;
	energytmp += *tmpptr * *tmpptr;
	if ((midarea < 0) && (areatmp >= halfarea)) {
	    midarea = tmpptr - dataptr;
	}
	if ((midenergy < 0) && (energytmp >= halfenergy)) {
	    midenergy = tmpptr - dataptr;
	}
	if ((midarea > 0) && (midenergy > 0)) {
	    /* break out of loop when both are assigned */
	    break;
	}
	tmpptr++;
    }
    
    /* assign the values */
    /* set the start and end index. Add one to all indeces to change from 0
     * base C to 1 base matlab */
    eptr[0] = startptr - dataptr + 1;
    eptr[1] = endptr - dataptr + 1 ;
    eptr[2] = peak;
    eptr[3] = peakind + 1;
    eptr[4] = area;
    eptr[5] = midarea + 1;
    eptr[6] = energy;
    eptr[7] = midenergy + 1;
    eptr[8] = maxthresh;
    eptr[9] = firstaboveptr- dataptr + 1;;
}


void Usage()
{
    mexErrMsgTxt("Usage: events = extractevents(signal, threshold, baseline, minspacing, allowstartend)\n\t see helpfile for more details.");
}




