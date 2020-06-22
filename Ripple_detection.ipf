#pragma rtGlobals=3	// Use modern global access method.

//Function detects ripple events (3xSD amplitude, 3 oscillations), sorts them according to mean frequency (<250Hz for SWRs and >250Hz for FRs), calculates their number, amlitude, oscillation frequency, etc. 

Function Ripple_detection()

SVAR dataname // filename for LFP recording
SVAR path // symbolic path to save HFB times to server
NVAR wave_number
SVAR list

string hfb_copy


Variable main_threshold
Variable local_threshold

Variable sampletime = 0.0002 //sample time, sec - 5KHz
Variable baseline
Variable eventstart, eventstartf
string/g status

Variable eventend, eventendf
Variable ripple, ripple250, ripple600
variable entropy
Variable pnumber
Variable rippleend,ripplenumber, ripplewidth
Variable localmax, peak, local_peak
Variable eventwidth, eventmax
Variable tracewidth
Variable n, k, j, i, c
Variable traceduration
String int="intermediate"
String int1="intermediate1"
String int2="intermediate2"
String rfreq="ripple frequency"
String rmode="ripple mode"
String rampl="ripple amplitude"
String rdr="ripple duration"
Variable avrfr, avramp, avrdur

string temp, temp2, temp3

string eventmode = "HFBfullmode_"+ num2str(wave_number)

make/O/N = 100000 $eventmode = NaN

string eventfreq =dataname[0,strlen(dataname)-3] + "_rippleMeanFreq" + dataname[strlen(dataname)-2,strlen(dataname)-1]//= "HFBfullfreq_"+ num2str(wave_number)

make/O/N = 100000 $eventfreq = NaN

string evententropy = "HFBfullentropy_"+ num2str(wave_number)

make/O/N = 100000 $evententropy = NaN

string eventindex = "HFBfullindex_"+ num2str(wave_number)

make/O/N = 100000 $eventindex = NaN


wave emode = $eventmode
wave eentr = $evententropy
wave eindex = $eventindex
wave efreq = $eventfreq

wave stdev = $"StDev"

wave HFB250f = $"HFB_250_freq"
wave HFB250dur  = $"HFB_250_duration"
wave HFB250occ = $"HFB_250_occupation"
wave HFB600f = $"HFB_600_freq"
wave HFB600dur  = $"HFB_600_duration"
wave HFB600occ = $"HFB_600_occupation"
wave HFB250ind = $"HFB_250_index"
wave HFB600ind = $"HFB_600_index"

wave HFB250pf = $"HFB_250_peakF"
wave HFB600pf = $"HFB_600_peakF"

wave HFB250mf = $"HFB_250_meanF"
wave HFB600mf = $ "HFB_600_meanF"

wave HFB250pA = $"HFB_250_peakAmp"
wave HFB600pA = $"HFB_600_peakAmp"

wave HFB250entropy = $"HFB_250_entropy"
wave HFB600entropy = $"HFB_600_entropy"


string t250times = dataname[0,strlen(dataname)-3] + "_125to250times" + dataname[strlen(dataname)-2,strlen(dataname)-1]
make/O/N=(10000,2) $t250times = NaN
wave times250 = $t250times

string t600times = dataname[0,strlen(dataname)-3] + "_250to600times" + dataname[strlen(dataname)-2,strlen(dataname)-1]
make/O/N=(10000,2) $t600times = NaN
wave times600 = $t600times




Make/O/n=100000 $"rfreq250" = NaN
Make/O/n=100000 $"rampl250" = NaN
Make/O/n=100000 $"rmode250" = NaN
Make/O/n=100000 $"rdur250" = NaN
Make/O/n=100000 $"rentr250" = NaN
Make/O/n=100000 $"rindex250" = NaN
//
wave rfr250 = $"rfreq250"
wave ramp250 = $"rampl250"
wave rmode250 = $"rmode250"
wave rentr250 = $"rentr250"
wave rdur250 = $"rdur250"
wave rind250 = $"rindex250"


Make/O/n=100000 $"rfreq600" = NaN
Make/O/n=100000 $"rampl600" = NaN
Make/O/n=100000 $"rmode600" = NaN
Make/O/n=100000 $"rdur600" = NaN
Make/O/n=100000 $"rentr600" = NaN
Make/O/n=100000 $"rindex600" = NaN
//
wave rfr600 = $"rfreq600"
wave ramp600 = $"rampl600"
wave rmode600 = $"rmode600"
wave rentr600 = $"rentr600"
wave rdur600 = $"rdur600"
wave rind600 = $"rindex600"


Wave ripplefreq = $"Ripple_frequency"

Wave ripplemode = $"Ripple_mode"

Wave rippleindex = $"Ripple_index"

Wave rippleentropy = $"Ripple_entropy"




//filter the signal for 125-600Hz band

wave main = $dataname

duplicate/o $dataname, $"filtered"

wave main = $"filtered"

FilterFIR/LO={0.12,0.122,500}/HI={0.023,0.025,500} main

//end of filtering



variable trim = 0

wavestats/q main
pnumber=V_npnts


trim = ceil(log(v_npnts)/log(2))   //find nearest power of 2



if(2^trim > v_npnts)  //find max power of 2 to pad to
	insertpoints v_npnts,(2^trim-v_npnts), main   //pad to nearest power of 2 with 0s
endif


hilbertTransform/DEST=ht main

duplicate/o ht, $"envelope"

wave hw = $"envelope"
hw = real(r2polar(cmplx(main[p], ht[p])))


smooth 20,hw  //4ms gaussian smooth

wavestats/q hw

deletepoints pnumber, (v_npnts-pnumber), hw    //get rid of padded zeros
deletepoints pnumber, (v_npnts-pnumber), main


variable ripplewidth_sum = 0

string execstring
string psdname
variable SD
variable comp

variable crossings = 0
variable freq

Wave hfbcopy
wave stdev = $"StDev"

Silent 1


// set thresholds from HT wave

WaveStats/Q/R=[0,pnumber] hw
SD = V_sdev
main_threshold = v_avg + SD*3  //all events exceeding 3xSDs
baseline = V_avg
stdev[wave_number] = SD  // save SD value for comparison

WaveStats/Q/R=[0,pnumber] main  
local_threshold=V_rms

traceduration=pnumber*sampletime


Differentiate main/D=differ 

n=0
k=0
j=0
ripple=0
ripple250 = 0
ripple600 = 0

i = 0


findlevels/P/Q/D=up/EDGE=1 hw,	main_threshold    //get all up crossings
crossings = V_LevelsFound
findlevels/P/Q/D=down/EDGE=2 hw,	main_threshold   // get all down crossings
variable crossingsDown = V_LevelsFound

eventendf = 0

do 
	
	do
	
		if((up[i] > eventendf)||(i>crossings))  //make sure the new event start is past the previous event
				break
		endif
		
		i+=1
	
	while(1)  
	
	
	eventstart = up[i]

	k = eventstart
	
	do
		k = k - 1
		
		if(k<0)   
			break
		endif
		
	while(hw[k]>baseline)   //go backwards until find 1xSD up crossing as start of event
	
	eventstartf = k
	
	k = 0
	
	
	do    //find stop of event	
	
		k = k+1
		
		if(k>=crossingsDown)   //breakpoint to detect last event
			break
		endif
		
	while(down[k] < eventstart)  //until next down crossing following event start
	
	eventend = down[k]
		

	k = down[k]
	

	do    // find 1xSD crossing point
	
		k = k + 1
		
		if(k>=pnumber)
			break
		endif
		
	while(hw[k] > baseline)  //go forward until find 1xSD crossing as end of event

	eventendf = k
	

			ripplewidth=(eventendf-eventstartf)*sampletime*1000  //ripple width, in ms
//		
			 	findlevels/M=0.002/P/D=dif_cross/R=[eventstart,eventend]/Q/EDGE=2 differ,0   //count positive peaks within putative event
			 	peak = V_LevelsFound
			
			 	k = 0
			 	
			 	
		 		do   //secondary check to make sure all peaks are real
		 	
				local_peak = peak
		 	
			 	if(main[dif_cross[k]]	< local_threshold )   //if peak is not positive, delete/substract
			 	 	peak = peak-1
			 	endif 	
			 	
			 	k+=1
			 	
			 	while(k < local_peak)
			 		

			 		
			 		if ((peak>=3)&(eventendf>eventstartf))  // if event has minimum 3 oscillations 
	 	 
	 	 					ripple+=1  //count ripple
			
				 			wavestats/q/R=[eventstartf,eventendf] hw
				 				
				 			localmax = v_max   //get max amplitude
	
			
					 			hfb_copy = dataname + "_HFB" + num2str(ripple)
								duplicate/O/R=[eventstartf-625,eventendf+625] main, $hfb_copy  
								
								wave hfbw = $hfb_copy
																
								psdname = fPowerSpectralDensity(hfbw,256,"hamming",1)   //get PSD
								
								wave psdw = $psdname
					
								wavestats/q psdw
													
																				
								eindex[j] = area(psdw,250,600)/area(psdw)   //calculate ripple index ( % of FR/total)
																
								//calculate mean frequency
								duplicate/o psdw,psdx
								psdx = x
								
								duplicate/o psdw,$"meanf"
								wave mnf = $"meanf"
								
								execstring ="psdx*" + psdname
								setformula mnf, execstring
								doupdate
//								
								 freq = sum(mnf)/sum(psdw)
								 if(freq>0)
									efreq[j] = freq
								 endif
								
								setformula mnf,""   //delete dependency
								doupdate
								
								
								//calculate spectral entropy
								
								temp = "psd_n"
								
								duplicate/o psdw, $temp
								
								execstring =  temp + ":= " + psdname + "/" + "area(" + psdname + ")"  //build command: normalize PSD to PDF
							
								execute execstring  //execute command
								
								
								wave psdw = $temp
					
								temp2 = "entr"
								
								duplicate/o psdw, $temp2
								 
								wave psentr = $temp2
								
								execstring = temp2 + ":=" + temp + "*log(" + temp + ")/Log(2)"
								
								execute execstring
								
								eentr[j] = -sum(psentr)*10  //record entropy value
								
								smooth 20, $temp
								
								wavestats/q $temp
								
									
								emode[j] = v_maxloc   //get peak power location
							
								
						
						
						
						
						if(efreq[j]<250)   //if slow ripple
								
								times250[ripple250][0] = eventstartf*sampletime
								times250[ripple250][1] = eventendf*sampletime
								
								rdur250[ripple250] = (eventendf-eventstartf)*sampletime*1000
								ramp250[ripple250] = localmax
								rind250[ripple250] = eindex[j]
								rfr250[ripple250] = efreq[j]
								rmode250[ripple250] = emode[j]
								rentr250[ripple250] = eentr[j]
						
								ripple250 = ripple250 + 1
							
								
						else     //if fast ripple
						
								times600[ripple600][0] = eventstartf*sampletime
								times600[ripple600][1] = eventendf*sampletime
								
								rdur600[ripple600] = (eventendf-eventstartf)*sampletime*1000
								ramp600[ripple600] = localmax
								rind600[ripple600] = eindex[j]
								rfr600[ripple600] = efreq[j]
								rmode600[ripple600] = emode[j]
								rentr600[ripple600] = eentr[j]
						
								ripple600 = ripple600 + 1
								
	
						endif  //end of ripple sorting and saving
						
						
								killwaves psentr 
								
								killwaves psdw
								
								killwaves psdx
								
							      killwaves $psdname  
								
							      killwaves hfbw
							   

					 	
		 					j = j+1
		 					
			 		endif  //end of peak > 3
				
			
	//end of ripple event
							 	
		i = i+ 1
		
while(i<crossings)

Redimension/N=(ripple) emode, eindex, eentr, efreq
Redimension/N=(ripple250) rdur250,ramp250,rind250,rfr250,rmode250,rentr250
Redimension/N=(ripple600) rdur600,ramp600,rind600,rfr600,rmode600,rentr600


if (ripple>0)  //if any ripples detected in trace
	
	wavestats/q efreq
	ripplefreq[wave_number] = V_avg
	print "mean freq: ", v_avg
	wavestats/q emode
	ripplemode[wave_number] = V_avg
	print "peak freq: ", v_avg
	wavestats/q eindex
	rippleindex[wave_number] = V_avg
	wavestats/q eentr
	rippleentropy[wave_number] = V_avg
	
	
	//250HFB averaging and saving
	
	HFB250f[wave_number] = (ripple250/(traceduration/60)) // ripple freq per minute
	
	wavestats/q rdur250
	HFB250dur[wave_number] = V_avg // average ripple duration, ms
	print "average slow duration: ", v_avg
	HFB250occ[wave_number] = (V_sum/(1000*traceduration))*100 //% occupation
	wavestats/q rmode250
	HFB250pf[wave_number] = V_avg
	wavestats/q rfr250
	HFB250mf[wave_number] = V_avg
	wavestats/q ramp250
	HFB250pA[wave_number] = V_avg
	wavestats/q rentr250
	HFB250entropy[wave_number]= V_avg
	wavestats/q rind250
	HFB250ind[wave_number]= V_avg
	
	
	//600HFB averaging and saving
	
	HFB600f[wave_number] = (ripple600/(traceduration/60)) // ripple freq per minute
	
	wavestats/q rdur600
	HFB600dur[wave_number] = V_avg  // average ripple duration, ms
	HFB600occ[wave_number] = (V_sum/(1000*traceduration))*100 //% occupation
	print "average fast duration: ", v_avg
	wavestats/q rmode600
	HFB600pf[wave_number] = V_avg
	wavestats/q rfr600
	HFB600mf[wave_number] = V_avg
	wavestats/q ramp600
	HFB600pA[wave_number] = V_avg
	wavestats/q rentr600
	HFB600entropy[wave_number]= V_avg
	wavestats/q rind600
	HFB600ind[wave_number]= V_avg
	
		
endif

print "Counted HFBs:", ripple
print ripple250, " slow and ", ripple600, " fast"


save/O/J/P=HFBtimesSave times250     //save event timestamps to server
save/O/J/P=HFBtimesSave times600   


end

				

		
		
		
		
		






