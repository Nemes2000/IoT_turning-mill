function [x, fTs] = fInverseFourierTransform(X, N)
%function [x, fTs] = fInverseFourierTransform(X, N)
% 	Inverse Fourier Transform is O(n*n)! opposed to ifft's O(n*log(n))
% NOTICE: must fTs/Ts to get the actual Frequency values!
%
% N is the number of samples to be restored
% 	the analysed  frequency spectrum is F = (1/Ts) * fTs; fTs=(0:L)'/(2*L), where L = length(X)-1;

%https://ccrma.stanford.edu/~jos/st/Matrix_Formulation_DFT.html
%https://blogs.uoregon.edu/seis/wiki/unpacking-the-matlab-fft/


	%N = length(x);
	%n = (0:(N-1));
	%k = (-N/2:(N/2-1))'/N; %no point in going above the Nyquist frequency
	%M = exp(-i*2*pi*k*n/N);
	%X = M*x;

	% better when the number of restoration points is explicitly defined as an input parameter
	if(~exist('N') || (~isempty(N) &&  ~isfinite(N))), N=length(X); end
	L = length(X); %number of frequency components
	n = (0:(N-1)); %sample time in Ts; n=(t-min(t))./Ts

%obsolete:fTs = (0:(N-1))'/N; %this is not frequency, but f*Ts
%obsolete: if(~exist('iFrequencyStepRefinement') || (~isempty(iFrequencyStepRefinement) &&  ~isfinite(iFrequencyStepRefinement))), iFrequencyStepRefinement=1; end
%obsolete:fTs = (0:1/iFrequencyStepRefinement:(N-1/iFrequencyStepRefinement))'/N;
%obsolete:fTs = (0:(L-1))'/L; %this is not frequency, but f*Ts; more simple format
%obsolete:fTs = (-L/2:(L/2-1))'/L; %this is not frequency, but f*Ts; no meaning in going above the Nyquist frequency fTs=1/2
%!!!FAILED INVERSE for L = L-1; %will start indexing from 0 :: fTs = (0:L)'/(2*L)
%!!!FAILED INVERSE for fTs = (0:L)'/(2*L); %this is not frequency, but f*Ts; in range [0:0.5] as no point in working with symmetric negative frequencies
%!!!FAILED INVERSE for LinvM = transp(exp(i*2*pi*fTs*n))/(L+1); %for M = exp(-i*2*pi*fTs*n); [N,L+1]

	fTs = [(0:(L/2)),(floor(-L/2+1):(-1))].'/L; %this is not frequency, but f*Ts; in "fftshift" order

	LinvM = transp(exp(i*2*pi*fTs*n))/L; %for M = exp(-i*2*pi*fTs*n); [N,L]
	x = real(LinvM*X);

	%M = exp(-i*2*pi*fTs*n);
	%pinvM = pinv(M);
	%x2 = real(pinvM*X); %X = M*x; pinvM*X = pinvM*M*x = 1*x
	%max(max(abs(x-x2))) % 2.1372e-15 vs 0.0043 for the failed inverse

return

%%% TST inverse
%clr
%N = 16
%n = (0:(N-1)); %sample time in Ts; n=(t-min(t))./Ts
%iFrequencyStepRefinement = 1 %3, 1/3
%L = floor(N*iFrequencyStepRefinement); %number of refined frequencies to check
%%!!!FAILED INVERSE (even with (*/(L+1))) for fTs = (0:L)'/(2*L); %this is not frequency, but f*Ts; in range [0:0.5] as no point in working with symmetric negative frequencies
%%fTs = [(0:ceil(L/2-1)),(ceil(-L/2):(-1))].'/L; %this is not frequency, but f*Ts; in "fftshift" order
%fTs = [(0:(L/2)),(floor(-L/2+1):(-1))].'/L; %this is not frequency, but f*Ts; in "fftshift" order
%M = exp(-i*2*pi*fTs*n); %size [N,L]
%
%RinvM = transp(exp(i*2*pi*fTs*n))/N;
%max(max(abs(diag(abs(M*RinvM))'-ones(L,1)))) %1.1102e-16
%
%LinvM = transp(exp(i*2*pi*fTs*n))/(L);
%max(max(abs(diag(abs(LinvM*M))'-ones(L,1)))) %1.1102e-16
%
%pinvM = pinv(M);
%max(max(abs(diag(abs(pinvM*M))'-ones(N,1)))) %2.2204e-16
%max(max(abs(diag(abs(M*pinvM*(L)/N))'-ones((L),1)))) %4.4409e-16
%
%
%N = 15
%n = (0:(N-1)); %sample time in Ts; n=(t-min(t))./Ts
%iFrequencyStepRefinement = 1 %3, 1/3
%L = floor(N*iFrequencyStepRefinement); %number of refined frequencies to check
%%!!!FAILED INVERSE (even with (*/(L+1))) for fTs = (0:L)'/(2*L); %this is not frequency, but f*Ts; in range [0:0.5] as no point in working with symmetric negative frequencies
%%fTs = [(0:ceil(L/2-1)),(ceil(-L/2):(-1))].'/L; %this is not frequency, but f*Ts; in "fftshift" order
%fTs = [(0:(L/2)),(floor(-L/2+1):(-1))].'/L; %this is not frequency, but f*Ts; in "fftshift" order
%M = exp(-i*2*pi*fTs*n); %size [N,L]
%
%RinvM = transp(exp(i*2*pi*fTs*n))/N;
%max(max(abs(diag(abs(M*RinvM))'-ones(L,1)))) %1.1102e-16
%
%LinvM = transp(exp(i*2*pi*fTs*n))/(L);
%max(max(abs(diag(abs(LinvM*M))'-ones(L,1)))) %2.2204e-16
%
%pinvM = pinv(M);
%max(max(abs(diag(abs(pinvM*M))'-ones(N,1)))) %2.2204e-16
%max(max(abs(diag(abs(M*pinvM*(L)/N))'-ones((L),1)))) %4.4409e-16
%

%% see fFourierTransform.m TST !!!

%% LEARN tstMyIFT.m !

%% TST
clr
iFigureN = 1;

%% create signal
disp(' ')
disp('----------')
disp('- signal -')

Tsampling = 1
PhaseShift = 11*Tsampling
PeriodDuration = 17
Frequency = 1/PeriodDuration %0.0588
Amplitude = 1
baseNsamples = 2^7 %128 %have N = power of 2
t = (0:(baseNsamples-1))'*Tsampling;
%t = (0:(baseNsamples))'*Tsampling; %to play around with even samples
Sig = Amplitude*cos(2*pi*Frequency*(t-PhaseShift));
RMS_signal = fRMS(Sig) %0.7093

%fDebugSTOP()

disp('Signal = Amplitude*cos(2*pi*Frequency*(t-PhaseShift))')
fMyFigure(iFigureN); plot(t,Sig), title('Signal')
iFigureN = iFigureN+1;

Sig_mean = mean(Sig); %-0.0360
y = Sig - Sig_mean; %adjusting for the DC component to become 0

%% FT with refined frequency spectrum
iFrequencyStepRefinement = 5; %keep it even to play around with even samples
[y_MyFFT, vMyFTs, vMyFTs_max, vMyFTs_max_isNyquist] = fFourierTransform(y,iFrequencyStepRefinement);
%	fTs = [(0:(L/2)),(floor(-L/2+1):(-1))].'/L; %this is not frequency, but f*Ts; in "fftshift" order
%	fTs_max = fTs(floor(L/2)+1);
%	fTs_max_isNyquist = (mod(L,2)==0);

MyNsamples = length(y)*iFrequencyStepRefinement; %analysed number a frequencies
												% == size(y_MyFFT,1)
Fs = 1/Tsampling; %==1							%sampling frequency
Fmax = Fs/2 %==0.5								%max observable 'Nyquist' frequency
vMyFTs_max %==0.5								% analysed highest frequency
vMyFTs_max_isNyquist %true						% is FTs_max same as Nyquist frequency?
MydFns = (1/MyNsamples) %==0.0020				%analysed frequency steps == 1/L
												% == (vMyFTs(2)-vMyFTs(1))/Tsampling
MyF_2sided = (2*Fmax)*[0:MydFns:1/2, (-1/2+MydFns*(mod(MyNsamples,2)==0)):MydFns:-MydFns]';	%analysed frequency spectrum for complete two sided interval, each side max Nyquist=Fmax/2
% max(max(abs(vMyFTs-MyF2s))) == 0

MyF = (2*Fmax)*(0:MydFns:1/2)';	%analysed frequency spectrum for useful one sided interval

%MyiI1 = (1:fix(MyNsamples/2+1*(mod(MyNsamples,2)==0)))'; 	%indexes of the one sided FFT frequency interval
MyiI1 = 1+(0:(MyNsamples/2))';

find(vMyFTs == vMyFTs_max) == max(MyiI1)

%%iFT with complete refined frequency spectrum
[y_MyiFFT, vMyInvFTs] = fInverseFourierTransform(y_MyFFT,length(y));
max(max(abs(vMyInvFTs-vMyFTs))) == 0
RMS_MyiFFTofMyFFT = fRMS(Sig-(y_MyiFFT+Sig_mean)) %2.7210e-15



%% analytically restored all harmonics
freq_y_MyFFT = MyF;
amp_y_MyFFT = abs(y_MyFFT(MyiI1))/MyNsamples; % (/MyNsamples) for the complete spectrum, and NOT (/Nsamples) !
amp_y_MyFFT(2:end-1) = 2*amp_y_MyFFT(2:end-1); %NOT*2 DC, NOT*2 Nyquist

phase_y_MyFFT = angle(y_MyFFT(MyiI1));
period_y_MyFFT = fFiniteDiv0(1,freq_y_MyFFT);
phaseshift_y_MyFFT = mod(fFiniteDiv0(-phase_y_MyFFT,(2*pi*freq_y_MyFFT)),period_y_MyFFT);

restored_y = fRestoreHarmonicSignal(t, amp_y_MyFFT, freq_y_MyFFT, phaseshift_y_MyFFT);

disp('RESTORED Signal = Amplitude*cos(2*pi*Frequency*(t-PhaseShift))')
RMS_RestoredYofMyFFT = fRMS(Sig-(restored_y+Sig_mean)) 
% 2.8469e-15
delta_Sig_RestoredYofMyFFT=max(max(abs(Sig - restored_y-Sig_mean))) 
%7.0360e-15

fMyFigure(iFigureN); plot(t,restored_y+Sig_mean), title('RESTORED.Signal')
iFigureN = iFigureN+1;
fMyFigure(iFigureN); plot(t,Sig-restored_y-Sig_mean), title('error = Sig - RESTORED.Signal')
iFigureN = iFigureN+1;



%% analytically restored single dominant harmonic
[vMyDominantAmplitude, iMyDominantAmplitude] = max(amp_y_MyFFT*iFrequencyStepRefinement); %(*iFrequencyStepRefinement) for individual samples!
%findpeaks(amp_y_MyFFT*iFrequencyStepRefinement,'MinPeakHeight',0.2,'SortStr','descend'); %(*iFrequencyStepRefinement) for individual samples! 
vMyDominantFrequency = MyF(MyiI1(iMyDominantAmplitude));
vMyDominantPhase = phase_y_MyFFT(MyiI1(iMyDominantAmplitude));
vMyDominantPeriod = 1/vMyDominantFrequency;
vMyDominantPhaseShift = mod(-vMyDominantPhase./(2*pi*vMyDominantFrequency),vMyDominantPeriod);

vMyDominantPeriod, vMyDominantAmplitude, vMyDominantFrequency, vMyDominantPhaseShift 
% 17.0667, 0.9951, 0.0586, 10.7848
% 16.8261, 0.9870, 0.0594, 11.4968

PeriodDuration, Amplitude, Frequency, PhaseShift 
% 17, 1, 0.0588, 11


MySig_dominantrestored = vMyDominantAmplitude*cos(2*pi*vMyDominantFrequency*(t-vMyDominantPhaseShift))+Sig_mean;
RMS_AnalyticalMyFFT_dominantrestored = fRMS(Sig-MySig_dominantrestored) %0.1791

RMS_signal 
%0.1050


%% InverseFT of the filtered FT
%the second side for negative frequencies
iFTDominant = zeros(length(y_MyFFT),1);
iFTDominant(MyiI1(iMyDominantAmplitude)) = 1;
% do NOT take the mean component, not the Nyquist frequency
if((iMyDominantAmplitude > 1) && (iMyDominantAmplitude < (length(y_MyFFT)+2)))
	MyiI1dominant_negative = length(y_MyFFT)+2 -iMyDominantAmplitude;
	iFTDominant(MyiI1dominant_negative) = 1;
end
find(iFTDominant), [~, iFTDominant2] = findpeaks(abs(y_MyFFT),'MinPeakHeight',mean(abs(y_MyFFT))+3*std(abs(y_MyFFT)),'SortStr','descend')


% ATTENTION! must use: MyIFTT of iFrequencyStepRefinement*iFTDominant.*y_FFT_refined
warning('when filtering the y_FT=fFT(y), and using selected individual frequency peaks then must use y_filtered=IFT(iFrequencyStepRefinement*indices*y_FT)')
y_FTdominant_iFFT = fInverseFourierTransform(iFrequencyStepRefinement*iFTDominant.*y_MyFFT, length(y));
RMS_DominantMyIFFT_restored = fRMS(Sig-(y_FTdominant_iFFT+Sig_mean)) %0.1791