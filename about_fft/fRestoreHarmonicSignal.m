function signal_sumCos = fRestoreHarmonicSignal(t, vAmplitude, vFrequency, vPhase)
%function signal_sumCos = fRestoreHarmonicSignal(t, vAmplitude, vFrequency, vPhase)
% signal_sumCos = SUM vAmplitude(iHarmonic)*cospi(2*vFrequency(iHarmonic)*t +vPhase(iHarmonic));
% 	signal = Amp*cos(2*pi*Freq*(t) + Phase))

	signal_sumCos = zeros(size(t));
	for iHarmonic = 1 : length(vAmplitude)
		signal_sumCos = signal_sumCos + vAmplitude(iHarmonic)*cos(2*pi*vFrequency(iHarmonic)*t +vPhase(iHarmonic));
	end %for iHarmonic

return

%%TST
%
%clr
%
%Fs = 1000;            % Sampling frequency                    
%Ts = 1/Fs;            % Sampling period == t(2)-t(1) == min(diff(t)) != mean(diff(t))
%
%Ts = 1/1000;          % milli seconds Sampling period
%Fs = 1/Ts;            % Sampling frequency                    
%N = 1500;             % Number of signal samples "Length of signal"
%t = (0:N-1)'*Ts;       % Time vector
%
%vAmp = [1.2; 0.7];
%vFreq = [120; 50];
%vPeriodStep = 1./(vFreq*Ts);
%vPhaseShiftStep = [3; 13];
%vPhase = vPhaseShiftStep *Ts;
%vPhaseShift = (Ts/(2*pi))*vPhaseShiftStep./vFreq;
%
%% a signal containing a 50 Hz sinusoid of amplitude 0.7 and a 120 Hz sinusoid of amplitude 1.
%S = vAmp(1)*cos(2*pi*vFreq(1)*t +vPhase(1)) + vAmp(2)*sin(2*pi*vFreq(2)*t +vPhase(2));
%
%y = S;
%
%iLogTo = 1;
%iLogFileIDs = [1];
%iFigureN = 1;
%
%vMinRMSrefinement=[]; vFrequencyTsPrecision=[]; vMinAmplitude=[]; vMaxHarmonicsN=[];
%[vDominantAmplitudes, vDominantFrequencies, vDominantPhases, vDominantPhaseShifts, vDominantPhaseSteps, vDominantPeriodSteps, vDominantFrequencyTsRanges, y_means, y_harmonics, y_restored, y_restored_errors, y_restored_RMSs, iFigureN] ...
%= fDominantHarmonics(y, Ts, vMinRMSrefinement, vFrequencyTsPrecision, vMinAmplitude, vMaxHarmonicsN, iLogTo,iLogFileIDs,iFigureN);
%y_restored_RMSs
%%0.4951    1.0e-14*0.9255    1.0e-14*0.8415
%
%fHorizontalMerge(vDominantAmplitudes', vAmp)
%fHorizontalMerge(vDominantFrequencies', vFreq)
%fHorizontalMerge(vDominantPhases', vPhase)
%fHorizontalMerge(vDominantPhases'+pi/2, vPhase) %sin(x) = cos(x-pi/2)
%fHorizontalMerge(vDominantPhaseShifts', vPhaseShift)
%fHorizontalMerge(vDominantPhaseSteps', vPhaseShiftStep)
%fHorizontalMerge(vDominantPeriodSteps', vPeriodStep)
%
%signal_sumCos = fRestoreHarmonicSignal(t, vDominantAmplitudes, vDominantFrequencies, vDominantPhases);
%
%max(max(abs(S-signal_sumCos)))
%max(max(abs(y-signal_sumCos)))
%
%fMyFigure(iFigureN);iFigureN=iFigureN+1; plot(t,[S-signal_sumCos]);
%fMyFigure(iFigureN);iFigureN=iFigureN+1; plot(t,[y-signal_sumCos]);
%
%
%%Corrupt the signal with zero-mean white noise with a variance of 4.
%vNoise = randn(size(t));
%disp(['rms(vNoise)::',num2str(rms(vNoise))];
%%rms(vNoise)::0.97511
%y_noisy = S + vNoise; %signal range is ~2, noise range is ~1, so 50%...
%disp('signal range is ~2, noise range is ~1, so 50%...')
%
%[vDominantAmplitudes_noisy, vDominantFrequencies_noisy, vDominantPhases_noisy, vDominantPhaseShifts_noisy, vDominantPhaseSteps_noisy, vDominantPeriodSteps_noisy, vDominantFrequencyTsRanges_noisy, y_means_noisy, y_harmonics_noisy, y_restored_noisy, y_restored_errors_noisy, y_restored_RMSs_noisy, iFigureN] ...
%= fDominantHarmonics(y_noisy, Ts, vMinRMSrefinement, vFrequencyTsPrecision, vMinAmplitude, vMaxHarmonicsN, iLogTo,iLogFileIDs,iFigureN);
%y_restored_RMSs_noisy
%%1.0821    0.9747    0.9710
%%% already the second iteration reaches < rms(vNoise)::0.97511
%
%fHorizontalMerge(vDominantAmplitudes_noisy', vAmp)
%fHorizontalMerge(vDominantFrequencies_noisy', vFreq)
%fHorizontalMerge(vDominantPhases_noisy', vPhase)
%fHorizontalMerge(vDominantPhases_noisy'+pi/2, vPhase) %sin(x) = cos(x-pi/2)
%fHorizontalMerge(vDominantPhaseShifts_noisy', vPhaseShift)
%fHorizontalMerge(vDominantPhaseSteps_noisy', vPhaseShiftStep)
%fHorizontalMerge(vDominantPeriodSteps_noisy', vPeriodStep)
%
%signal_sumCos_noisy = fRestoreHarmonicSignal(t, vDominantAmplitudes, vDominantFrequencies, vDominantPhases);
%
%max(max(abs(S-signal_sumCos_noisy)))
%max(max(abs(y_noisy-signal_sumCos_noisy)))
%
%fMyFigure(iFigureN);iFigureN=iFigureN+1; plot(t,[S-signal_sumCos_noisy]);
%fMyFigure(iFigureN);iFigureN=iFigureN+1; plot(t,[y_noisy-signal_sumCos_noisy]);
%