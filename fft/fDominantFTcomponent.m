function [vDominantAmplitude, vDominantFrequency, vDominantPhase, vDominantPhaseShift, vDominantPhaseStep, vDominantPeriodStep, vDominantFrequencyTsRange] ...
= fDominantFTcomponent(x, Ts, vFrequencyTsRange)
%function [vDominantAmplitude, vDominantFrequency, vDominantPhase, vDominantPhaseShift, vDominantPhaseStep, vDominantPeriodStep, vDominantFrequencyTsRange] ...
%= fDominantFTcomponent(x, Ts, vFrequencyTsRange)
%%USES
%	iFrequencyStepRefinement = 1;
%	vTSFrequencyStep = vFrequencyTsRange(2)-vFrequencyTsRange(1);
%	vTsFrequencyies = [(0:vTSFrequencyStep:(L*vFrequencyTsRange(2))), ((floor(L*vFrequencyTsRange(1))+vTSFrequencyStep):vTSFrequencyStep:(-vTSFrequencyStep))].'/L;
%	[X, vFreqTs, vFreqTs_max, vFreqTs_max_isNyquist] = fFourierTransform(x, iFrequencyStepRefinement, vTsFrequencyies)
%
	iFrequencyStepRefinement=1;

	%L = length(x);
	%n = (0:(L-1));
	%k = (0:(L-1))'/L;
	%M = exp(-i*2*pi*k*n/L);
	%	vFreqTs = [(0:ceil(L/2-1)),(ceil(-L/2):(-1))].'/L; %this is not frequency, but f*Ts; in "fftshift" order
	%	i2pivFreqTs = -i*2*pi*vFreqTs; %performance optimisation
	%	M = exp(i2pivFreqTs*n);
	%X = M*x;

	L = size(x,1); %number of samples

	%% input
	if(~exist('Ts') || any(~isfinite(Ts)) || isempty(Ts)), Ts=1; end
	if(~exist('vFrequencyTsRange') || any(~isfinite(vFrequencyTsRange)) || isempty(vFrequencyTsRange)), vFrequencyTsRange=[-1/2+1/L,1/2]; end

	%% output
	vDominantAmplitude = NaN;
	vDominantFrequency = NaN;
	vDominantPhaseShift = NaN;
	vDominantFrequencyTsRange = NaN;
	vDominantPeriod = NaN;

	%% analysed (positive) frequency spectrum
	iPositiveSide_vFreqTs = (1:(L/2)+1);
	%F = vFreqTs(iPositiveSide_vFreqTs)/Ts;

	%% search spectrum in (frequency * Ts)
	if(~all(vFrequencyTsRange==[-1/2+1/L,1/2]))
		%L=10;vTsFrequencies=[(0:(L/2)),(floor(-L/2+1):(-1))].'/L, vFrequencyTsRange=[-1/2+1/L,1/2];Ts=1;
		vTsFrequencies = [vFrequencyTsRange(1):(vFrequencyTsRange(2)-vFrequencyTsRange(1))/(L-1):vFrequencyTsRange(2)].'; %vector length L
		iTsFreqMean = find(vTsFrequencies==0);
		if(~isempty(iTsFreqMean))
			iPositiveSide_vFreqTs = iTsFreqMean:L;
		else
			iPositiveSide_vFreqTs = 1:L;
		end
	else
		%%default vTsFrequencies=[(0:(L/2)),(floor(-L/2+1):(-1))].'/L
		vTsFrequencies = [];
	end %if(~all(vFrequencyTsRange==[-1/2+1/L,1/2]))
	[x_FT, vFreqTs, ~, vFreqTs_max_isNyquist, iTransformedSampleN] = fNonuniformFourierTransform(x, iFrequencyStepRefinement, vTsFrequencies);

	%% calculate amplitude spectrum
	vAmplitudeSpectrum = abs(x_FT)/iTransformedSampleN;
	% double all amplitude values that have a negative frequency pair
	%	double all except the mean component iPositiveSide_vFreqTs(1) and the Nyquist's frequency component iPositiveSide_vFreqTs(end) for vFreqTs_max_isNyquist
	vAmplitudeSpectrum(iPositiveSide_vFreqTs(2:end-1)) = 2 * vAmplitudeSpectrum(iPositiveSide_vFreqTs(2:end-1));
	if(~vFreqTs_max_isNyquist) %if last is not Nyquist's frequency
		vAmplitudeSpectrum(iPositiveSide_vFreqTs(end)) = 2 * vAmplitudeSpectrum(iPositiveSide_vFreqTs(end));
	end
	if(vFreqTs(1)~=0) %if first is not the mean component (for a special vFrequencyTsRange case)
		vAmplitudeSpectrum(iPositiveSide_vFreqTs(1)) = 2 * vAmplitudeSpectrum(iPositiveSide_vFreqTs(1));
	end

%	% find the one sided (positive) frequency position of harmonics with peak power peaks
%	vSensitivity = [1 1/3]; % peaks > mean+3*std
%	if(~isempty(vSensitivity) && (vSensitivity(1)>0) && (vSensitivity(2)>0))
%		vSignificantAmplitude = mean(vAmplitudeSpectrum(iPositiveSide_vFreqTs))/vSensitivity(2) +std(vAmplitudeSpectrum(iPositiveSide_vFreqTs))/vSensitivity(1);
%	else
%		vSignificantAmplitude = max(vAmplitudeSpectrum(iPositiveSide_vFreqTs))-eps;
%	end
%	[vDominantAmplitude, iSignificant] = findpeaks(vAmplitudeSpectrum(iPositiveSide_vFreqTs),'MinPeakHeight',vSignificantAmplitude,'SortStr','descend'); 

	%% find frequency index of the dominant harmonic, of maximal amplitude
	[vDominantAmplitude, iDominant] = max(vAmplitudeSpectrum(iPositiveSide_vFreqTs));

	%% frequency corresponding to dominant harmonic
	vDominantFrequency = vFreqTs(iDominant)/Ts; %==F(iDominant)

	%% period corresponding to dominant harmonic
	vDominantPeriodStep = 1./(vDominantFrequency*Ts);

	%% phase 'phi' corresponding to dominant harmonic in Amp*cos(2*pi*Freq*t +phi) % sin(t) = cos(t-pi/2)
	vDominantPhase = angle(x_FT(iDominant));

	%% phase time shift 'tshift' corresponding to dominant harmonic Amp*cos(2*pi*Freq*(t-tshift))
	%rounded to the nearest (2*pi*vDominantFrequency) multiplicand
	vDominantPhaseShift = mod(round(vDominantPhase/(2*pi*vDominantFrequency))*(2*pi*vDominantFrequency),2*pi)./(2*pi*vDominantFrequency);

	%% phase time shift 'tshift' in number of sampling times corresponding to dominant harmonic Amp*cos(2*pi*Freq*(t-tshift))
	%vDominantPhaseStep = 2*pi*vDominantFrequency.*vDominantPhaseShift/Ts
	%rounded to the nearest (Ts) multiplicand
	vDominantPhaseStep = mod(round(vDominantPhase/Ts)*Ts,2*pi)/Ts;

	%% set vDominantFrequencyTsRange the FrequencyTs range of certainty
	if(iDominant>1)
		iPrev = max(find(x(iPositiveSide_vFreqTs(1:(iDominant-1))))); % prev valid sample
		vDominantFrequencyTsRange(1) = vFreqTs(iPrev);
	else
		vDominantFrequencyTsRange(1) = vFreqTs(iDominant);
	end
	if(iDominant<length(iPositiveSide_vFreqTs))
		iNext = iDominant + min(find(x(iPositiveSide_vFreqTs(iDominant+1):end))); %  next valid sample
		vDominantFrequencyTsRange(2) = vFreqTs(iNext);
	else
		vDominantFrequencyTsRange(2) = vFreqTs(iDominant);
	end

return

%%%TST & learn !!!
%%
%% !!! % the mean part of the noise is immediately compensated for
%% !!! % seemingly the 0.98838_rms(vNoise) adds only a ~0.4191 (=1.4014-0.98234) points inrease to the signal rms
%%
%% !!! % the mean part of the noise can be compensated for with harmonics decomposition
%%mean(vNoise) :: 0.033133
%% !!! % the randomness caused rms part of the noise can NOT be compensated for with harmonics decomposition
%%rms(vNoise) :: 0.98838 =~~= rms(x_2_residual-mean(x_2_residual)) :: 0.98133
%%
%% !!! % the random noise degrades the harmonics identification by some rms_0.0719
%%
%% !!! % trying to identify the non existing 3rd harmonics does not benefit the identification quality
%
%
%clr
%iFigureN=1;
%
%Fs = 1000;			% Sampling frequency                    
%Ts = 1/Fs;			% Sampling period == t(2)-t(1) == min(diff(t)) != mean(diff(t))
%
%Ts = 1/1000;		% milliseconds Sampling period
%Fs = 1/Ts;			% Sampling frequency                    
%L = 1500;			% Number of signal samples "Length of signal"
%t = (0:L-1)'*Ts;	% Time vector
%
%vAmp = [1.2; 0.7]
%vFreq = [120; 50]
%vPeriodStep = 1./(vFreq*Ts) %[8.3333, 20.0000]
%vPhaseSteps = [3; 13] %= vPhase/Ts
%vPhase = vPhaseSteps *Ts %[0.0030, 0.0130]
%vPhaseShift = vPhase./(2*pi*vFreq) %1.0e-04*[0.0398, 0.4138]
%%vAmp = [1.2000, 0.7000]
%%vFreq = [120, 50]
%%vPeriodStep = [8.3333, 20.0000]
%%vPhaseSteps = [3, 13]
%%vPhase = [0.0030, 0.0130]
%%vPhaseShift = 1.0e-04*[0.0398, 0.4138]
%
%% a signal containing a 50 Hz sinusoid of amplitude 0.7 and a 120 Hz sinusoid of amplitude 1.
%vSignal = vAmp(1)*cos(2*pi*vFreq(1)*t +vPhase(1)) + vAmp(2)*sin(2*pi*vFreq(2)*t +vPhase(2));
%
%% base 2 freq signal
%x = vSignal;
%fMyFigure(iFigureN);iFigureN=iFigureN+1; plot(t,[vSignal]); title('Signal')
%disp(['mean(x) :: ',num2str(mean(x))])
%%mean(x) :: -2.0236e-16
%x_0 = x-mean(x);
%disp(['rms(x-mean(x)) :: ',num2str(rms(x_0))])
%%rms(x-mean(x)) :: 0.98234
%
%
%vFrequencyTsRange = [-1/2+1/L, 1/2] %[-0.5+1/L, 0.5] :: vector length is L
%%vFrequencyTsRange = [-0.4993    0.5000]
%vFrequencyTsPrecision = 1e-4;
%
%%vAmp = [1.2000, 0.7000]
%%vFreq = [120, 50]
%%vPhase = [0.0030, 0.0130]
%%vPhaseShift = 1.0e-04*[0.0398, 0.4138] = [3.9789e-06, 4.1380e-05]
%%vPhaseSteps = [3, 13]
%%vPeriodStep = [8.3333, 20]
%[vDominantAmplitude1, vDominantFrequency1, vDominantPhase1, vDominantPhaseShift1, vDominantPhaseStep1, vDominantPeriodStep1, vDominantFrequencyTsRange1] = fDominantFTcomponent(x_0, Ts, vFrequencyTsRange)
%%vDominantAmplitude1 = 1.2000 %OK
%%vDominantFrequency1 = 120 %OK
%%vDominantPhase1 = 0.0030 %OK
%%vDominantPhaseShift1 = 3.9789e-06 %OK
%%vDominantPhaseStep1 = 3 %OK
%%vDominantPeriodStep1 = 8.3333 %OK
%%vDominantFrequencyTsRange1 = 0.1193    0.1207
%
%x_0_restored = fRestoreHarmonicSignal(t, vDominantAmplitude1, vDominantFrequency1, vDominantPhase1);
%fMyFigure(iFigureN);iFigureN=iFigureN+1; plot(t,[x_0_restored]); title('x_0_restored')
%x_0_residual = x_0 - x_0_restored;
%fMyFigure(iFigureN);iFigureN=iFigureN+1; plot(t,[x_0_residual]); title('x_0_residual')
%disp(['mean(x_0_residual) :: ',num2str(mean(x_0_residual))])
%%mean(x_0_residual) :: -9.5849e-17
%x_1 = x_0_residual-mean(x_0_residual);
%disp(['rms(x_0_residual-mean(x_0_residual)) :: ',num2str(rms(x_1))])
%%rms(x_0_residual-mean(x_0_residual)) :: 0.49497
%
%%vAmp = [1.2000, 0.7000]
%%vFreq = [120, 50]
%%vPhase = [0.0030, 0.0130]
%%vPhaseShift = 1.0e-04*[0.0398, 0.4138] = [3.9789e-06, 4.1380e-05]
%%vPhaseSteps = [3, 13]
%%vPeriodStep = [8.3333, 20]
%[vDominantAmplitude2, vDominantFrequency2, vDominantPhase2, vDominantPhaseShift2, vDominantPhaseStep2, vDominantPeriodStep2, vDominantFrequencyTsRange2] = fDominantFTcomponent(x_1, Ts, vFrequencyTsRange)
%%vDominantAmplitude2 = 0.7000 %OK
%%vDominantFrequency2 = 50 %OK
%%sin(x)==cos(x-pi/2),cos(x)==sin(x+pi/2) :: vPhaseD = +pi/2; vPhaseShiftD=vPhaseD/(2*pi*vFreq)
%%vDominantPhase2 = -1.5578 %sin(x) :: vDominantPhase2+(pi/2) == 0.0130 %OK
%%vDominantPhaseShift2 = 0.0150 %sin(x) :: #TOTO(vDominantPhaseShift2+(pi/2)/(2*pi*vDominantFrequency2) == 0.0200)
%%vDominantPhaseStep2 = 4.7254e+03 %sin(x) :: mod(vDominantPhaseStep2+(pi/2)/Ts,2*pi/Ts) == 13 %OK
%%vDominantPeriodStep2 = 20 %OK
%%vDominantFrequencyTsRange2 = 0.0493    0.0507
%	
%x_1_restored = fRestoreHarmonicSignal(t, vDominantAmplitude2, vDominantFrequency2, vDominantPhase2);
%fMyFigure(iFigureN);iFigureN=iFigureN+1; plot(t,[x_1_restored]); title('x_1_restored')
%x_1_residual = x_1 - x_1_restored;
%fMyFigure(iFigureN);iFigureN=iFigureN+1; plot(t,[x_1_residual]); title('x_1_residual')
%disp(['mean(x_1_residual) :: ',num2str(mean(x_1_residual))])
%%mean(x_1_residual) :: 2.9134e-16
%x_2 = x_1_residual-mean(x_1_residual);
%disp(['rms(x_1_residual-mean(x_1_residual)) :: ',num2str(rms(x_2))])
%%rms(x_1_residual-mean(x_1_residual)) :: 9.2476e-15
%
%[vDominantAmplitude3, vDominantFrequency3, vDominantPhase3, vDominantPhaseShift3, vDominantPhaseStep3, vDominantPeriodStep3, vDominantFrequencyTsRange3] = fDominantFTcomponent(x_2, Ts, vFrequencyTsRange)
%%vDominantAmplitude3 = 5.4349e-15 % pointless...
%%vDominantFrequency3 = 120
%%vDominantPhase3 = -1.7356
%%vDominantPhaseShift3 = 0.0060
%%vDominantPhaseStep3 = 4.5476e+03
%%vDominantPeriodStep3 = 8.3333
%%vDominantFrequencyTsRange3 = 0.1193    0.1207
%
%x_2_restored = fRestoreHarmonicSignal(t, vDominantAmplitude3, vDominantFrequency3, vDominantPhase3);
%fMyFigure(iFigureN);iFigureN=iFigureN+1; plot(t,[x_2_restored]); title('x_2_restored')
%x_2_residual = x_2 - x_2_restored;
%fMyFigure(iFigureN);iFigureN=iFigureN+1; plot(t,[x_2_residual]); title('x_2_residual')
%disp(['mean(x_2_residual) :: ',num2str(mean(x_2_residual))])
%%mean(x_2_residual) :: 6.6422e-31
%x_3 = x_2_residual-mean(x_2_residual);
%disp(['rms(x_2_residual-mean(x_2_residual)) :: ',num2str(rms(x_3))])
%%rms(x_2_residual-mean(x_2_residual)) :: 8.4113e-15
%fMyFigure(iFigureN);iFigureN=iFigureN+1; plot(t,[x_3]); title('3stepError')
%
%
%x_restored = ...
%	+ mean(x)+fRestoreHarmonicSignal(t, vDominantAmplitude1, vDominantFrequency1, vDominantPhase1)...
%	+ mean(x_0_residual)+fRestoreHarmonicSignal(t, vDominantAmplitude2, vDominantFrequency2, vDominantPhase2)...
%	+ mean(x_1_residual)+fRestoreHarmonicSignal(t, vDominantAmplitude3, vDominantFrequency3, vDominantPhase3);
%
%fMyFigure(iFigureN);iFigureN=iFigureN+1; plot(t,[x, x_restored]); title('x, x_restored')
%fMyFigure(iFigureN);iFigureN=iFigureN+1; plot(t,[x-x_restored]); title('x - x_restored')
%
%disp(['rms(x - x_restored) :: ',num2str(rms(x - x_restored))])
%%rms(x - x_restored) :: 8.4122e-15
%rms(x_2_residual) % 8.4113e-15
%rms(x_3) % 8.4113e-15
%
%
%%%Corrupt the signal with zero-mean white noise with a variance of 4.
%vNoise = randn(size(t));
%disp(['mean(vNoise) :: ',num2str(mean(vNoise))])
%%mean(vNoise) :: 0.033133
%disp(['rms(vNoise) :: ',num2str(rms(vNoise))])
%%rms(vNoise) :: 0.98838
%
%x = vSignal + vNoise;
%fMyFigure(iFigureN);iFigureN=iFigureN+1; plot(t,[vSignal,x]); title('Signal, noiseCorupted')
%%mean(x) :: -2.0236e-16
%disp(['mean(x) :: ',num2str(mean(x))])
%%mean(x) :: 0.033133
%x_0 = x-mean(x);
%%rms(x-mean(x)) :: 0.98234
%disp(['rms(x-mean(x)) :: ',num2str(rms(x_0))])
%%rms(x-mean(x)) :: 1.4014
%
%
%% !!! % the mean part of the noise is immediately compensated for
%% !!! % seemingly the 0.98838_rms(vNoise) adds only a ~0.4191 (=1.4014-0.98234) points inrease to the signal rms
%
%
%%vDominantAmplitude1 = 1.2000
%%vDominantFrequency1 = 120
%%vDominantPhase1 = 0.0030
%%vDominantPhaseShift1 = 3.9789e-06
%%vDominantPhaseStep1 = 3.0000
%%vDominantPeriodStep1 = 8.3333
%%vDominantFrequencyTsRange1 = 0.1193    0.1207
%[vDominantAmplitude1, vDominantFrequency1, vDominantPhase1, vDominantPhaseShift1, vDominantPhaseStep1, vDominantPeriodStep1, vDominantFrequencyTsRange1] = fDominantFTcomponent(x_0, Ts, vFrequencyTsRange)
%%vDominantAmplitude1 = 1.1965
%%vDominantFrequency1 = 120
%%vDominantPhase1 = 2.1118e-04
%%vDominantPhaseShift1 = 2.8009e-07
%%vDominantPhaseStep1 = 0.2112
%%vDominantPeriodStep1 = 8.3333
%%vDominantFrequencyTsRange1 = 0.1193    0.1207
%
%x_0_restored = fRestoreHarmonicSignal(t, vDominantAmplitude1, vDominantFrequency1, vDominantPhase1);
%fMyFigure(iFigureN);iFigureN=iFigureN+1; plot(t,[x_0_restored]); title('x_0_restored')
%x_0_residual = x_0 - x_0_restored;
%fMyFigure(iFigureN);iFigureN=iFigureN+1; plot(t,[x_0_residual]); title('x_0_residual')
%%mean(x_0_residual) :: -9.5849e-17
%disp(['mean(x_0_residual) :: ',num2str(mean(x_0_residual))])
%%mean(x_0_residual) :: -8.941e-17
%x_1 = x_0_residual-mean(x_0_residual);
%%rms(x_0_residual-mean(x_0_residual)) :: 0.49497
%disp(['rms(x_0_residual-mean(x_0_residual)) :: ',num2str(rms(x_1))])
%%rms(x_0_residual-mean(x_0_residual)) :: 1.1172
%
%%vDominantAmplitude2 = 0.7000 %OK
%%vDominantFrequency2 = 50 %OK
%%sin(x)==cos(x-pi/2),cos(x)==sin(x+pi/2) :: vPhaseD = +pi/2; vPhaseShiftD=vPhaseD/(2*pi*vFreq)
%%vDominantPhase2 = -1.5578 %sin(x) :: vDominantPhase2+(pi/2) == 0.0130 %OK
%%vDominantPhaseShift2 = 0.0150 %sin(x) :: #TOTO(vDominantPhaseShift2+(pi/2)/(2*pi*vDominantFrequency2) == 0.0200)
%%vDominantPhaseStep2 = 4.7254e+03 %sin(x) :: mod(vDominantPhaseStep2+(pi/2)/Ts,2*pi/Ts) == 13 %OK
%%vDominantPeriodStep2 = 20 %OK
%%vDominantFrequencyTsRange2 = 0.0493    0.0507
%[vDominantAmplitude2, vDominantFrequency2, vDominantPhase2, vDominantPhaseShift2, vDominantPhaseStep2, vDominantPeriodStep2, vDominantFrequencyTsRange2] = fDominantFTcomponent(x_1, Ts, vFrequencyTsRange)
%%vDominantAmplitude2 = 0.7391
%%vDominantFrequency2 = 50
%%vDominantPhase2 = -1.5549
%%vDominantPhaseShift2 = 0.0151
%%vDominantPhaseStep2 = 4.7283e+03
%%vDominantPeriodStep2 = 20
%%vDominantFrequencyTsRange2 = 0.0493    0.0507
%
%x_1_restored = fRestoreHarmonicSignal(t, vDominantAmplitude2, vDominantFrequency2, vDominantPhase2);
%fMyFigure(iFigureN);iFigureN=iFigureN+1; plot(t,[x_1_restored]); title('x_1_restored')
%x_1_residual = x_1 - x_1_restored;
%fMyFigure(iFigureN);iFigureN=iFigureN+1; plot(t,[x_1_residual]); title('x_1_residual')
%%mean(x_1_residual) :: 2.9134e-16
%disp(['mean(x_1_residual) :: ',num2str(mean(x_1_residual))])
%%mean(x_1_residual) :: 2.9958e-16
%x_2 = x_1_residual-mean(x_1_residual);
%%rms(x_1_residual-mean(x_1_residual)) :: 9.2476e-15
%disp(['rms(x_1_residual-mean(x_1_residual)) :: ',num2str(rms(x_2))])
%%rms(x_1_residual-mean(x_1_residual)) :: 0.98743
%
%%vDominantAmplitude3 = 5.4349e-15 % pointless...
%%vDominantFrequency3 = 120
%%vDominantPhase3 = -1.7356
%%vDominantPhaseShift3 = 0.0060
%%vDominantPhaseStep3 = 4.5476e+03
%%vDominantPeriodStep3 = 8.3333
%%vDominantFrequencyTsRange3 = 0.1193    0.1207
%[vDominantAmplitude3, vDominantFrequency3, vDominantPhase3, vDominantPhaseShift3, vDominantPhaseStep3, vDominantPeriodStep3, vDominantFrequencyTsRange3] = fDominantFTcomponent(x_2, Ts, vFrequencyTsRange)
%%vDominantAmplitude3 = 0.1550 % this has relevance
%%vDominantFrequency3 = 450.6667
%%vDominantPhase3 = 2.1791
%%vDominantPhaseShift3 = 7.6954e-04
%%vDominantPhaseStep3 = 2.1791e+03
%%vDominantPeriodStep3 = 2.2189
%%vDominantFrequencyTsRange3 = 0.4500    0.4513
%
%x_2_restored = fRestoreHarmonicSignal(t, vDominantAmplitude3, vDominantFrequency3, vDominantPhase3);
%fMyFigure(iFigureN);iFigureN=iFigureN+1; plot(t,[x_2_restored]); title('x_2_restored')
%x_2_residual = x_2 - x_2_restored;
%fMyFigure(iFigureN);iFigureN=iFigureN+1; plot(t,[x_2_residual]); title('x_2_residual')
%%mean(x_2_residual) :: 6.6422e-31
%disp(['mean(x_2_residual) :: ',num2str(mean(x_2_residual))])
%%mean(x_2_residual) :: 3.2761e-17
%x_3 = x_2_residual-mean(x_2_residual);
%%rms(x_2_residual-mean(x_2_residual)) :: 8.4113e-15
%disp(['rms(x_2_residual-mean(x_2_residual)) :: ',num2str(rms(x_3))])
%%rms(x_2_residual-mean(x_2_residual)) :: 0.98133
%
%% !!! % the mean part of the noise can be compensated for with harmonics decomposition
%%mean(vNoise) :: 0.033133
%% !!! % the randomness caused rms part of the noise can NOT be compensated for with harmonics decomposition
%%rms(vNoise) :: 0.98838 =~~= rms(x_2_residual-mean(x_2_residual)) :: 0.98133
%
%
%fMyFigure(iFigureN);iFigureN=iFigureN+1; plot(t,[x_3]); title('3stepError')
%
%x_restored = ...
%	+ mean(x)+fRestoreHarmonicSignal(t, vDominantAmplitude1, vDominantFrequency1, vDominantPhase1)...
%	+ mean(x_0_residual)+fRestoreHarmonicSignal(t, vDominantAmplitude2, vDominantFrequency2, vDominantPhase2)...
%	+ mean(x_1_residual)+fRestoreHarmonicSignal(t, vDominantAmplitude3, vDominantFrequency3, vDominantPhase3);
%
%fMyFigure(iFigureN);iFigureN=iFigureN+1; plot(t,[x, x_restored]); title('x, x_restored')
%fMyFigure(iFigureN);iFigureN=iFigureN+1; plot(t,[x-x_restored]); title('x - x_restored')
%
%%rms(x - x_restored) :: 8.4122e-15
%disp(['rms(x - x_restored) :: ',num2str(rms(x - x_restored))])
%%rms(x - x_restored) :: 0.98133
%
%rms(x_2_residual) % 8.4113e-15
%rms(x_3) % 8.4113e-15
%
%%% harmonic component identification bias introduced by random noise
%
%x_biased = ...
%	+ mean(x)+fRestoreHarmonicSignal(t, vDominantAmplitude1, vDominantFrequency1, vDominantPhase1)...
%	+ mean(x_0_residual)+fRestoreHarmonicSignal(t, vDominantAmplitude2, vDominantFrequency2, vDominantPhase2);
%
%rms(x - x_biased - vNoise - mean(vNoise)) %0.0719
%
%%rms(x - x_biased) :: 8.4122e-15
%disp(['rms(x - x_biased) :: ',num2str(rms(x - x_biased))])
%%rms(x - x_biased) :: 0.98743
%
%disp(['rms(x - x_biased - vNoise - mean(vNoise)) :: ',num2str(rms(x - x_biased))])
%%rms(x - x_biased - vNoise - mean(vNoise)) %0.0719
%
%% !!! % the random noise degrades the harmonics identification by some rms_0.0719
%
%% !!! % trying to identify the non existing 3rd harmonics does not benefit the identification quality
%