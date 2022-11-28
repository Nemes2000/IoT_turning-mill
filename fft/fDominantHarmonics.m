function [vDominantAmplitudes, vDominantFrequencies, vDominantPhases, vDominantPhaseShifts, vDominantPhaseSteps, vDominantPeriodSteps, vDominantFrequencyTsRanges, y_means, y_harmonics, y_restored, y_restored_errors, y_restored_RMSs, iFigureN] ...
= fDominantHarmonics(y, Ts, vMinRMSrefinement, vFrequencyTsPrecision, vMinAmplitude, vMaxHarmonicsN, iLogTo,iLogFileIDs,iFigureN)
%function [vDominantAmplitudes, vMinRMSrefinement, vDominantFrequencies, vDominantPhaseShifts, vDominantFrequencyTsRanges, vDominantPeriods, y_means, y_harmonics, y_restored, y_restored_e, y_restored_RMS, iFigureN] ...
%= fDominantHarmonics(y, Ts, vFrequencyTsPrecision, vMinAmplitude, vMaxHarmonicsN, iLogTo,iLogFileIDs,iFigureN)
%
%%USES 
%	[vDominantAmplitude, vDominantFrequency, vDominantPhase, vDominantPhaseShift, vDominantPhaseStep, vDominantPeriodStep, vDominantFrequencyTsRange, iFigureN] ...
%	= fDominantFTcomponent(x, Ts, vFrequencyTsRange, iFigureN)
%	(IF iFigureN) x_DominantHarmonic = fRestoreHarmonicSignal(t, vDominantAmplitude, vDominantFrequency, vDominantPhase);
%	y_restored_DominantHarmonic = fRestoreHarmonicSignal(t, vDominantAmplitude, vDominantFrequency, vDominantPhaseShift);
%	y_means = [y_means, fFiniteRowMean(y_residual)];
%	y_harmonics = [y_harmonics, y_restored_DominantHarmonic];
%   y_restored = sum(y_means)+sum(y_harmonics,2);
%
%% default
%	Ts = 1;
%	vFrequencyTsPrecision = 0.05; end %0.1*(Nyquist freq=1/2)/(Ts=1)
%	 to refine each dominant frequency while(((vDominantFrequencyTsRange(2)-vDominantFrequencyTsRange(1)) > vFrequencyTsPrecision))
%	vMinAmplitude = 0.067; %0.1*rms(whiteNoise=0.67)
%	 to look for additional dominant frequencies while(vDominantAmplitude>vMinAmplitude)
%	vMinRMSrefinement = 0.1; % each iteration has to reduce the RMS at lest by 10%
%	vMaxHarmonicsN = 100; % maximum number of harmonics to evaluate

	L = size(y,1); %number of samples

	%% input
	if(~exist('Ts') || any(~isfinite(Ts)) || isempty(Ts)), Ts=1; end
	if(~exist('vMinRMSrefinement') || any(~isfinite(vMinRMSrefinement)) || isempty(vMinRMSrefinement)), vMinRMSrefinement=0.1; end %min 10% improvements required for each iteration
	if(~exist('vFrequencyTsPrecision') || any(~isfinite(vFrequencyTsPrecision)) || isempty(vFrequencyTsPrecision)), vFrequencyTsPrecision=0.05/Ts; end %0.1*(Nyquist freq=1/2)/(Ts=1)
	if(~exist('vMinAmplitude') || any(~isfinite(vMinAmplitude)) || isempty(vMinAmplitude)), vMinAmplitude=0.1*fFiniteRowRMS(y); end %0.067 = 0.1*rms(whiteNoise=0.67)
	if(~exist('vMaxHarmonicsN') || any(~isfinite(vMaxHarmonicsN)) || isempty(vMaxHarmonicsN)), vMaxHarmonicsN=100; end

	if(~exist('iFigureN') || isempty(iFigureN) || ~isfinite(iFigureN)), iFigureN=[]; end
	pInitLogging
%	%keeps existing iLogTo and iLogFileIDs
%	if(~exist('iLogTo') || isempty(iLogTo) || ~isfinite(iLogTo)), iLogTo=0; end %no logging
%	iLogToScreen = 0b0001;
%	iLogToError = 0b0010;
%	iLogToReport = 0b0100;
%	iLogToResult = 0b1000;
%	iLogToAll = bitand(iLogTo, (iLogToScreen+iLogToReport+iLogToResult)); %stdout + "_report.txt" + "_result.csv" IF not overruled by the initial iLogTo
%	iLogTo = bitand(iLogTo, (iLogToScreen+iLogToReport)); %stdout + "_report.txt" IF not overruled by the initial iLogTo
%	if(~exist('iLogFileIDs') || isempty(iLogFileIDs) || any(~isfinite(iLogFileIDs))), iLogFileIDs=[1,2]; end
%
%	fLog(sLog, iLogToScreen, iLogFileIDs); %iLogTo & stdout
%	fLog(sLog, iLogToError, iLogFileIDs); %iLogTo & stderr
%	fLog(sLog, iLogTo, iLogFileIDs); %iLogTo & (stdout + "_report.txt")
%	fLog(sLog, iLogToAll, iLogFileIDs); %iLogTo & (stdout + "_report.txt" + "_result.csv")

	%% output
	y_means = [];
	vDominantAmplitudes = [];
	vDominantFrequencies = [];
	vDominantPhases = [];
	vDominantPhaseShifts = [];
	vDominantPhaseSteps = [];
	vDominantPeriodSteps = [];
	vDominantFrequencyTsRanges = [];
	y_harmonics = [];
	y_restored_errors = []; % y - sum_i(y_means(i)+y_harmonics(:,i))
	y_restored_RMSs = [];

	%% init
	t = (0:(L-1))'*Ts; %sample time index in Ts steps; n=(t-min(t))./Ts
	y_residual = y;
	y_residual_RMS = fFiniteRowRMS(y_residual);
	y_restored = zeros(size(y), 'like',y);

	%% repeat until the desired amplitude precision is reached
	bRepeatAmp = true;
	iDominantHarmonic = 1;
	while(bRepeatAmp) %repeat

		sLog = append('. evaluating the [',num2str(iDominantHarmonic),']th dominant harmonic component');
		fLog(sLog, iLogTo, iLogFileIDs); %
		sLog = append('.. remaining signal has RMS :: [',num2str(y_residual_RMS),']');
		fLog(sLog, iLogTo, iLogFileIDs); %

		%%init new dominant harmonic search
		y_means = [y_means, NaN];
		vDominantAmplitudes 		= [vDominantAmplitudes, NaN];
		vDominantFrequencies 		= [vDominantFrequencies, NaN];
		vDominantPhases 			= [vDominantPhases, NaN];
		vDominantPhaseShifts 		= [vDominantPhaseShifts, NaN];
		vDominantPhaseSteps 		= [vDominantPhaseSteps, NaN];
		vDominantPeriodSteps 		= [vDominantPeriodSteps, NaN];
		vDominantFrequencyTsRanges 	= [vDominantFrequencyTsRanges,[NaN, NaN]'];
		y_harmonics = [y_harmonics, nan(size(y))];
		y_restored_errors = [y_restored_errors, nan(size(y))]; % y - sum_i(y_means(i)+y_harmonics(:,i))
		y_restored_RMSs = [y_restored_RMSs, NaN];

		%% init frequency refinement
		y_mean = fFiniteRowMean(y_residual);
		y_residual = y_residual -y_mean;

		%% repeat until the desired frequency precision is reached
		bRepeatFreq = true;
		iFreqRefinement = 0;
		%% for each new frequency component start with the full domain search
		vTsFrequencyies = [];
		vFrequencyTsRange=[-1/2+1/L,1/2]; %default vTsFrequencies=[(0:(L/2)),(floor(-L/2+1):(-1))].'/L
		while(bRepeatFreq) %repeat

			sLog = append('.. evaluating the [',num2str(iFreqRefinement),']th frequency refinement');
			fLog(sLog, iLogTo, iLogFileIDs); %
			sLog = append('... the analysed frequency*Ts interval is :: [',num2str(vFrequencyTsRange),']');
			fLog(sLog, iLogTo, iLogFileIDs); %

			%% the fDominantFTcomponent
			[vDominantAmplitude, vDominantFrequency, vDominantPhase, vDominantPhaseShift, vDominantPhaseStep, vDominantPeriodStep, vDominantFrequencyTsRange] = fDominantFTcomponent(y_residual, Ts, vFrequencyTsRange);

			%% iterate until desired precision
			vFrequencyTsRange = vDominantFrequencyTsRange;
			bRepeatAmp = (vDominantAmplitude > vMinAmplitude);
			bRepeatFreq = ((vDominantFrequencyTsRange(2)-vDominantFrequencyTsRange(1)) > vFrequencyTsPrecision);
			if(bRepeatAmp || iDominantHarmonic==1 || iFreqRefinement==0)
				%%update new result / store first result
				y_means(end) = y_mean;
				vDominantAmplitudes(end) = vDominantAmplitude;
				vDominantFrequencies(end) = vDominantFrequency;
				vDominantPhases(end) = vDominantPhase;
				vDominantPhaseShifts(end) = vDominantPhaseShift;
				vDominantPhaseSteps(end) = vDominantPhaseStep;
				vDominantPeriodSteps(end) = vDominantPeriodStep;
				vDominantFrequencyTsRanges(:,end) = vDominantFrequencyTsRange';
				y_harmonics(:,end)  = fRestoreHarmonicSignal(t, vDominantAmplitudes(end), vDominantFrequencies(end), vDominantPhases(end));
				y_restored_errors(:,end) = y -(sum(y_means)+sum(y_harmonics,2));
				y_restored_RMSs(end) = fFiniteRowRMS(y_restored_errors(:,end));

				sLog = append('.... has updated parameters of the dominant harmonic');
				fLog(sLog, iLogTo, iLogFileIDs); %
			else
				sLog = append('.... NO parameter update of the dominant harmonic');
				fLog(sLog, iLogTo, iLogFileIDs); %
			end %if(bRepeatAmp)

			%% log
			if(bRepeatAmp)
				if(iFreqRefinement==0)
					sLog = append('. found a [',num2str(iDominantHarmonic),']th significant dominant harmonic component');
				else
					sLog = append('. updating the [',num2str(iDominantHarmonic),']th significant dominant harmonic component');
				end
			else
				sLog = append('. the [',num2str(iDominantHarmonic),']th dominant harmonic component is NOT significant');
			end
			fLog(sLog, iLogTo, iLogFileIDs); %
			sLog = append('.. evaluated the [',num2str(iFreqRefinement),']th frequency refinement');
			fLog(sLog, iLogTo, iLogFileIDs); %
			sLog = append('... reached frequency precision :: [',num2str(vDominantFrequencyTsRange(2)-vDominantFrequencyTsRange(1)),']; precision limit[',num2str(vFrequencyTsPrecision),']');
			fLog(sLog, iLogTo, iLogFileIDs); %
			sLog = append('... dominant amplitude : [',num2str(vDominantAmplitude),']; amplitude limit[',num2str(vMinAmplitude),']');
			fLog(sLog, iLogTo, iLogFileIDs); %
			sLog = append('... dominant frequency : [',num2str(vDominantFrequency),']');
			fLog(sLog, iLogTo, iLogFileIDs); %
			sLog = append('... dominant phase : [',num2str(vDominantPhase),']');
			fLog(sLog, iLogTo, iLogFileIDs); %
			sLog = append('... dominant phase in Ts steps : [',num2str(vDominantPhaseStep),']');
			fLog(sLog, iLogTo, iLogFileIDs); %
			sLog = append('... dominant period in Ts steps : [',num2str(vDominantPeriodStep),']');
			fLog(sLog, iLogTo, iLogFileIDs); %

			iFreqRefinement = iFreqRefinement +1;
		end %while(bRepeatFreq && bRepeatAmp)

		%% restore signal by dominant harmonics
		y_restored = sum(y_means)+sum(y_harmonics,2);
		y_residual = y -y_restored;
		y_residual_RMS_prev = y_residual_RMS;
		y_residual_RMS = fFiniteRowRMS(y_residual);
		iDominantHarmonic = iDominantHarmonic +1;

		bRepeatAmp = (vDominantAmplitudes(end) > vMinAmplitude);
		bRepeatAmp = bRepeatAmp && (iDominantHarmonic <= vMaxHarmonicsN);
		bRepeatAmp = bRepeatAmp && ((y_residual_RMS_prev-y_residual_RMS) > (vMinRMSrefinement*y_residual_RMS_prev) || (y_residual_RMS-y_residual_RMS_prev) < (vMinRMSrefinement*y_residual_RMS_prev));
		%requires minimal RMS gain of vMinRMSrefinement
		% AND allows oscillations (RMS increase) of not more than vMinRMSrefinement

	end %while(bRepeatAmp)

	sLog = append('. the final residual signal has RMS :: [',num2str(y_residual_RMS),']');
	fLog(sLog, iLogTo, iLogFileIDs); %

	if(~isempty(iFigureN))
		fMyFigure(iFigureN, 'DominantHarmonics');
			sPlotType = '-';
			iLineWidth = 2;
			sLabelX = ['t'];
			sLabelY = ['y'];
		subplot(2,1,1)
			fMyMinAxisMargin(1,1);
			sTitle = ['Signal & Restored'];
			cLegend = {'signal','harmonic'};
			vMean = fFiniteRowMean(y);
			vStd = fFiniteRowStd(y);
			fMyPlot(t, [y, y_restored], sPlotType, iLineWidth, sTitle, sLabelX, sLabelY, cLegend, vMean, vStd);
		subplot(2,1,2)
			fMyMinAxisMargin(1,1);
			sTitle = ['Residual'];
			cLegend = {'residual'};
			vMean = fFiniteRowMean(y_residual);
			vStd = fFiniteRowStd(y_residual);
			fMyPlot(t, y_residual, sPlotType, iLineWidth, sTitle, sLabelX, sLabelY, cLegend, vMean, vStd);
		iFigureN = iFigureN+1;
	end

	if(~isempty(iFigureN))
		fMyFigure(iFigureN, 'Dominant Harmonics');
			sPlotType = '-';
			iLineWidth = 2;
			sLabelX = ['t'];
			sLabelY = ['signal'];
		subplot(4,1,1)
			fMyMinAxisMargin(1,1);
			sTitle = ['Signal'];
			cLegend = {'signal'};
			vMean = fFiniteRowMean(y);
			vStd = fFiniteRowStd(y);
			fMyPlot(t, [y], sPlotType, iLineWidth, sTitle, sLabelX, sLabelY, cLegend, vMean, vStd);
		subplot(4,1,2)
			fMyMinAxisMargin(1,1);
			sTitle = ['Signal-mean & Harmonics'];
			cLegend = {'centred signal','harmonics'};
			vMean = fFiniteRowMean(y-fFiniteRowMean(y));
			vStd = fFiniteRowStd(y-fFiniteRowMean(y));
			fMyPlot(t, [y-fFiniteRowMean(y), y_harmonics], sPlotType, iLineWidth, sTitle, sLabelX, sLabelY, cLegend, vMean, vStd);
		subplot(4,1,3)
			fMyMinAxisMargin(1,1);
			sTitle = ['Signal & Restored(means+DominantHarmonics)'];
			cLegend = {'signal', 'restored'};
			vMean = fFiniteRowMean(y_restored);
			vStd = fFiniteRowStd(y_restored);
			fMyPlot(t, [y, y_restored], sPlotType, iLineWidth, sTitle, sLabelX, sLabelY, cLegend, vMean, vStd);
		subplot(4,1,4)
			fMyMinAxisMargin(1,1);
			sTitle = ['Residual of restored signal'];
			cLegend = {'residual'};
			vMean = fFiniteRowMean(y_restored_errors(:,end));
			vStd = fFiniteRowStd(y_restored_errors(:,end));
			fMyPlot(t, y_restored_errors(:,end), sPlotType, iLineWidth, sTitle, sLabelX, sLabelY, cLegend, vMean, vStd);
		iFigureN = iFigureN+1;
	end

return

%%%TST
%clr
%%SampleMul = 1; %RMS_harmonics=0.0165
%%SampleMul = 10; %RMS_harmonics=0.0430
%SampleMul = 100; %RMS_harmonics=0.0041
%Ts = 0.01;
%t=(0:Ts:SampleMul-Ts)'; disp(['t=(0:',num2str(Ts),':',num2str(SampleMul),'-',num2str(Ts),')''']) %#100*SampleMul
%
%%% single harmonic
%%vSignal = cos(2*pi*t); 
%%disp(['vSignal = cos(2*pi*t)'])
%
%%% 2 harmonics
%vSignal = 0.3*sin(2*pi*0.5*(t- 0.7)) + 0.7*cos(2*pi*0.3*(t- 0.5)); 
%disp(['vSignal = 0.3*sin(2*pi*0.5*(t- 0.7)) + 0.7*cos(2*pi*0.3*(t- 0.5))'])
%disp(['cos.vPhaseShift : vPhase :: [0.3*sin(2*pi*0.5*(t- 0.7))] : -2*pi*0.5*0.7+pi/2 =',num2str(-2*pi*0.5*0.7+pi/2),'; [0.7*cos(2*pi*0.3*(t- 0.5))] : -2*pi*0.3*0.5 =',num2str(-2*pi*0.3*0.5)])
%
%%%Corrupt the signal with zero-mean white noise with a variance of 0.1
%vNoiseVariance = 0.25; %can cope with 50% random signal variance
%vNoise = vNoiseVariance*randn(size(t));
%disp(['mean(vNoise) :: ',num2str(mean(vNoise))])
%disp(['rms(vNoise) :: ',num2str(rms(vNoise))])
%%mean(vNoise) :: 0.033133
%%rms(vNoise) :: 0.98838
%
%%x = vSignal;
%x = vSignal + vNoise;
%
%disp(['mean(x) :: ',num2str(mean(x))])
%x_sig = x-mean(x);
%disp(['rms(x-mean(x)) :: ',num2str(rms(x_sig))])
%
%iLogTo = 1;
%iLogFileIDs = [1];
%iFigureN = 1;
%vFrequencyTsRange = [-1/2-1/(SampleMul/Ts), 1/2];
%
%vMinRMSrefinement = 0.1; % each iteration has to reduce the RMS at lest by 10%
%vMaxHarmonicsN = 100; % maximum number of harmonics to evaluate
%vFrequencyTsPrecision = 0.1*(1/2); %super precise minimal frequency delta (mind the Nyquist freq)
%vMinAmplitude = 0.1*rms(x_sig); %rational minimal amplitude (mind the signal RMS)
%
%[vDominantAmplitudes, vDominantFrequencies, vDominantPhases, vDominantPhaseShifts, vDominantPhaseSteps, vDominantPeriodSteps, vDominantFrequencyTsRanges, y_means, y_harmonics, y_restored, y_restored_errors, y_restored_RMSs, iFigureN] = fDominantHarmonics(x, Ts, vMinRMSrefinement, vFrequencyTsPrecision, vMinAmplitude, vMaxHarmonicsN, iLogTo,iLogFileIDs,iFigureN);
%
%vDominantAmplitudes, vDominantFrequencies, vDominantPhases, vDominantPhaseShifts, vDominantPhaseSteps, vDominantPeriodSteps, vDominantFrequencyTsRanges, y_means, y_restored_RMSs
%figure(iFigureN), iFigureN=iFigureN+1;
%	plot(t,y_harmonics),title('harmonics')
%figure(iFigureN), iFigureN=iFigureN+1;
%	plot(t,y_restored),title('restored')
%figure(iFigureN), iFigureN=iFigureN+1;
%	plot(t,y_restored_errors),title('restore_errors')
%
%x_restore = sum(y_means) + sum(y_harmonics,2);
%figure(iFigureN), iFigureN=iFigureN+1;
%	plot(t,[x, x_restore]), title('x, restored')
%figure(iFigureN), iFigureN=iFigureN+1;
%	plot(t,[x - x_restore]), title('delta')
%disp(['rms(x - x_restore)) :: ',num2str(rms(x - x_restore))])
%
%%results in:
%%
%%t=(0:0.01:100-0.01)'
%%vSignal = 0.3*sin(2*pi*0.5*(t- 0.7)) + 0.7*cos(2*pi*0.3*(t- 0.5))
%%cos.vPhaseShift : vPhase :: [0.3*sin(2*pi*0.5*(t- 0.7))] : -2*pi*0.5*0.7+pi/2 =-0.62832; [0.7*cos(2*pi*0.3*(t- 0.5))] : -2*pi*0.3*0.5 =-0.94248
%%mean(vNoise) :: 0.0036629
%%rms(vNoise) :: 0.25229
%%mean(x) :: 0.0036629
%%rms(x-mean(x)) :: 0.59461
%%
%%. evaluating the [1]th dominant harmonic component
%%.. remaining signal has RMS :: [0.59465]
%%.. evaluating the [0]th frequency refinement
%%... the analysed frequency*Ts interval is :: [-0.4999		 0.5]
%%.... has updated parameters of the dominant harmonic
%%. found a [1]th significant dominant harmonic component
%%.. evaluated the [0]th frequency refinement
%%... reached frequency precision :: [0.0002]; precision limit[0.0005]
%%... dominant amplitude : [0.70024]; amplitude limit[0.059461]
%%... dominant frequency : [0.3]
%%... dominant phase : [-0.94079]
%%... dominant phase in Ts steps : [534.3185]
%%... dominant period in Ts steps : [333.3333]
%%. evaluating the [2]th dominant harmonic component
%%.. remaining signal has RMS :: [0.32925]
%%.. evaluating the [0]th frequency refinement
%%... the analysed frequency*Ts interval is :: [-0.4999		 0.5]
%%.... has updated parameters of the dominant harmonic
%%. found a [2]th significant dominant harmonic component
%%.. evaluated the [0]th frequency refinement
%%... reached frequency precision :: [0.0002]; precision limit[0.0005]
%%... dominant amplitude : [0.2992]; amplitude limit[0.059461]
%%... dominant frequency : [0.5]
%%... dominant phase : [2.5065]
%%... dominant phase in Ts steps : [251]
%%... dominant period in Ts steps : [200]
%%. evaluating the [3]th dominant harmonic component
%%.. remaining signal has RMS :: [0.25227]
%%.. evaluating the [0]th frequency refinement
%%... the analysed frequency*Ts interval is :: [-0.4999		 0.5]
%%.... has updated parameters of the dominant harmonic
%%. the [3]th dominant harmonic component is NOT significant
%%.. evaluated the [0]th frequency refinement
%%... reached frequency precision :: [0.0002]; precision limit[0.0005]
%%... dominant amplitude : [0.014754]; amplitude limit[0.059461]
%%... dominant frequency : [49.22]
%%... dominant phase : [-2.7513]
%%... dominant phase in Ts steps : [353.3185]
%%... dominant period in Ts steps : [2.0317]
%%. the final residual signal has RMS :: [0.25206]


%%% TEST
%	% Callisto's angular position is measured in minutes of arc.
%	% Missing data due to cloudy conditions are specified using NaNs. 
%	% The first observation is dated January 15. Generate a datetime array of observation times.
%clr
%	yg = [10.5 NaN 11.5 10.5 NaN NaN NaN -5.5 -10.0 -12.0 -11.5 -12.0 -7.5 ...
%		NaN NaN NaN NaN 8.5 12.5 12.5 10.5 NaN NaN NaN -6.0 -11.5 -12.5 ...
%		-12.5 -10.5 -6.5 NaN 2.0 8.5 10.5 NaN 13.5 NaN 10.5 NaN NaN NaN ...
%		-8.5 -10.5 -10.5 -10.0 -8.0]';
%
%	obsv = datetime(1610,1,14+(1:length(yg)));
%
%	plot(yg,'o')
%		ax = gca;
%		nights = [1 18 32 46];
%		ax.XTick = nights;
%		ax.XTickLabel = char(obsv(nights));
%		grid
%
%	% Estimate the power spectrum of the data using plomb. 
%	% 	Specify an oversampling factor of 10. 
%	% 	Express the resulting frequencies in inverse days.
%	[pxx,f] = plomb(yg,obsv,[],10,'power');
%	Ts = 86400; %=24*60*60 seconds in a day; for daily observations
%	f = f*Ts;
%
%	% Use findpeaks to determine the location of the only prominent peak of the spectrum. Plot the power spectrum and show the peak.
%	[pk,f0] = findpeaks(pxx,f,'MinPeakHeight',10);
%
%	figure(1)
%	plot(f,pxx,f0,pk,'o')
%	xlabel('Frequency (day^{-1})')
%	title('Power Spectrum and Prominent Peak')
%	grid
%
%	% Determine Callisto's orbital period (in days) as the inverse of the frequency of maximum energy. The result differs by less than 1% from the value published by NASA.
%	Period = 1/f0
%	% Period = 16.6454
%
%	NASA = 16.6890184;
%	PercentDiscrep = (Period-NASA)/NASA*100
%	% PercentDiscrep = -0.2613
%
%	%% DO no juggling with dates and days/seconds; simple take periods in number of samples
%
%	%top 3 peaks; plot to figure(2)
%	[iPeriod_top, vPeriodogramPower_top] = fPeriods(yg,[],3,[],2)
%
%iLogTo = 1;
%iLogFileIDs = [1];
%iFigureN = 1;
%vFrequencyTsRange = [-1/2-1/(1/Ts), 1/2];
%vMinRMSrefinement = 0.1; % each iteration has to reduce the RMS at lest by 10%
%vMaxHarmonicsN = 100; % maximum number of harmonics to evaluate
%vFrequencyTsPrecision = 0.1*(1/2)/Ts; %rational minimal frequency delta (mind the Nyquist freq)
%vMinAmplitude = 0.1*rms(yg); %rational minimal amplitude (mind the signal RMS)
%
%[vDominantAmplitudes, vDominantFrequencies, vDominantPhases, vDominantPhaseShifts, vDominantPhaseSteps, vDominantPeriodSteps, vDominantFrequencyTsRanges, y_means, y_harmonics, y_restored, y_restored_errors, y_restored_RMSs, iFigureN] = fDominantHarmonics(yg, Ts, vMinRMSrefinement, vFrequencyTsPrecision, vMinAmplitude, vMaxHarmonicsN, iLogTo,iLogFileIDs,iFigureN);
%1./vDominantFrequencies, vDominantPeriodSteps
%
%	% Determine Callisto's orbital period (in days) as the inverse of the frequency of maximum energy. The result differs by less than 1% from the value published by NASA.
%	Period_i = vDominantPeriodSteps(1)
%	% Period_i = 16.6484
%
%	NASA = 16.6890184;
%	PercentDiscrep = (Period_i-NASA)/NASA*100
%	% PercentDiscrep = -0.2432
%
%	%%same results by defaults
%[vDominantAmplitudes, vDominantFrequencies, vDominantPhases, vDominantPhaseShifts, vDominantPhaseSteps, vDominantPeriodSteps, vDominantFrequencyTsRanges, y_means, y_harmonics, y_restored, y_restored_errors, y_restored_RMSs, iFigureN] = fDominantHarmonics(yg, Ts);
%1./vDominantFrequencies, vDominantPeriodSteps
%
%	% Determine Callisto's orbital period (in days) as the inverse of the frequency of maximum energy. The result differs by less than 1% from the value published by NASA.
%	Period_i = vDominantPeriodSteps(1)
%	% = 16.6484
%