function [X, fTs, fTs_max, fTs_max_isNyquist] = fFourierTransform(x, iFrequencyStepRefinement, vTsFrequencyies)
%function [X, fTs, fTs_max, fTs_max_isNyquist] = fFourierTransform(x, optional iFrequencyStepRefinement=1, optional , vTsFrequencyies=[])
%% fFourierTransform is O(n*n)! opposed to fft's O(n*log(n))
% operates on a single columns of x
% uses single matrix M = exp(i2pifTs*n) operation
%	is a faster than fNonuniformFourierTransform()
% copes with non-uniform samples, when some are NaN as in missing data like x(isnan(x)) = 0;
% does NOT cope with condensed data of uneven sampling periods (use the greatest common denominator of sampling periods for the full span, and fill the missing data with NaN)
% NaN x values are numerically treated as being 0
% the analysed  frequency spectrum is f = (1/Ts) * fTs; fTs = [(0:(L/2)),(floor(-L/2+1):(-1))].'/L; where L=iFrequencyStepRefinement*N, and N=length(x)
%	fTs is not frequency, but f*Ts; in "fftshift" order%
% NOTICE: must fTs/Ts to get the actual Frequency values!
%	fTs_max presents the max processed frequency <= Nyquist's frequency(*Ts)
%	fTs_max_isNyquist is true for fTs_max == Nyquist's frequency(*Ts)
% NOTICE: must have harmonic_amplitude = 2*abs(X)/(N*iFrequencyStepRefinement); %instead of fft(x) common abs(X)/N
%	do exclude from *2 the amplitude corresponding to the mean(DC) (fTs(1)) and the Nyquist frequency (fTs(fix(MyNsamples/2+1*(mod(MyNsamples,2)==0))), when (mod(MyNsamples,2)==0))
% fTs = [(0:(L/2)),(floor(-L/2+1):(-1))].'/L; %this is not frequency, but f*Ts; in "fftshift" order
%	fTs = n_df/L
%	L=11;n_df=[(0:(L/2)),(floor(-L/2+1):(-1))],length(n_df),n_df_max=n_df(floor(L/2)+1),wmax_is_Nyquist=(mod(L,2)==0)
%	%n_df = 0     1     2     3     4     5    -5    -4    -3    -2    -1
%	%11, 5, false
%	L=11;n_df=[(0:(L/2)),(floor(-L/2+1):(-1))],length(n_df),n_df_max=n_df(floor(L/2)+1),wmax_is_Nyquist=(mod(L,2)==0)
%	%n_df = 0     1     2     3     4     5     6    -5    -4    -3    -2    -1
%	%12, 6, true
% fTs_max = fTs(floor(L/2)+1);
% fTs_max_isNyquist = (mod(L,2)==0);

%https://ccrma.stanford.edu/~jos/st/Matrix_Formulation_DFT.html
%https://blogs.uoregon.edu/seis/wiki/unpacking-the-matlab-fft/
	%N = length(x);
	%n = (0:(N-1));
	%k = (-N/2:(N/2-1))'/N; %no point in going above the Nyquist frequency
	%M = exp(-i*2*pi*k*n/N);
	%X = M*x;

	if(~exist('iFrequencyStepRefinement') || (~isempty(iFrequencyStepRefinement) &&  ~isfinite(iFrequencyStepRefinement))), iFrequencyStepRefinement=1; end
	if(~exist('vTsFrequencyies') || any(~isfinite(vTsFrequencyies)) || isempty(vTsFrequencyies)), vTsFrequencyies=[]; end

	[N,K] = size(x); %samples
	n = (0:(N-1)); %sample time in Ts; n=(t-min(t))./Ts
	L = floor(N*iFrequencyStepRefinement); %frequencies to check

	if(isempty(vTsFrequencyies))
%fTs = (0:(N-1))'/N; %this is not frequency, but f*Ts
%fTs = (0:1/iFrequencyStepRefinement:(N-1/iFrequencyStepRefinement))'/N;
%fTs = (0:(L-1))'/L; %this is not frequency, but f*Ts; more simple format
%fTs = (-L/2:(L/2-1))'/L; %this is not frequency, but f*Ts; no meaning in going above the Nyquist frequency fTs=1/2
%!!!FAILED INVERSE for fTs = (0:L)'/(2*L); %this is not frequency, but f*Ts; in range [0:0.5]
		fTs = [(0:(L/2)),(floor(-L/2+1):(-1))].'/L; %this is not frequency, but f*Ts; in "fftshift" order
	else
		fTs = vTsFrequencyies;
	end
%	iPositiveSide_fTs = (1:(N/2)+1);
%	iNegativeSide_fTs = (floor(N/2)+2:N)


	x(isnan(x)) = 0; %ignoring NaN samples
	fTs_max = fTs(floor(L/2)+1); %the max analysed frequency component
	fTs_max_isNyquist = (mod(L,2)==0); % is the max analysed frequency component the same as the Nyquist's frequency == 1/2

	i2pifTs = -i*2*pi*fTs; %performance optimisation
	M = exp(i2pifTs*n);
	X = M*x;

return

%%%TST
%clr
%
%%% EVEN samples
%
%t=(0:0.01:10-0.01)'; disp(['t=(0:0.01:10-0.01)'''])
%x = 0.0 + 0.3*sin(2*pi*0.5*(t- 0.7)) + 0.7*cos(2*pi*0.3*(t- 0.5)); disp(['x = 0.3*sin(2*pi*0.5*(t- 0.7)) + 0.7*cos(2*pi*0.3*(t- 0.5))'])
%x_mean = mean(x);
%disp(['length(x) :: ',num2str(length(x))])
%disp(['mean(x) :: ',num2str(x_mean)])
%disp(['rms(x) :: ',num2str(rms(x))])
%disp([num2str(rms(fft(x)-fFourierTransform(x))),' :: rms(fft(x)-fFourierTransform(x))'])
%%t=(0:0.01:10-0.01)'
%%x = 0.3*sin(2*pi*0.5*(t- 0.7)) + 0.7*cos(2*pi*0.3*(t- 0.5))
%%length(x) :: 1000
%%mean(x) :: -4.7351e-17
%%rms(x) :: 0.53852
%%7.7128e-13 :: rms(fft(x)-fFourierTransform(x))
%
%disp([num2str(rms(x-ifft(fft(x)))),' :: rms(x-ifft(fft(x)))'])
%%1.6568e-16 :: rms(x-ifft(fft(x)))
%
%iFrequencyStepRefinement = 1;
%X_FT = fFourierTransform(x-x_mean,iFrequencyStepRefinement);
%disp(['length(X_FT) :: ',num2str(length(X_FT))])
%N = length(x);
%IX_FT = fInverseFourierTransform(X_FT,N)+x_mean;
%disp(['length(IX_FT) :: ',num2str(length(IX_FT))])
%disp([num2str(rms(x-fInverseFourierTransform(fFourierTransform(x,iFrequencyStepRefinement),length(x)))),' :: rms(x-fInverseFourierTransform(fFourierTransform(x,iFrequencyStepRefinement),length(x)))'])
%%length(X_FT) :: 1000
%%length(IX_FT) :: 1000
%%2.4371e-14 :: rms(x-fInverseFourierTransform(fFourierTransform(x,iFrequencyStepRefinement),length(x)))
%
%%D=[x,IX_FT, x-IX_FT, IX_FT./x, [0;diff(IX_FT./x)]];
%
%iFrequencyStepRefinement = 2;
%X_FT_refined = fFourierTransform(x-x_mean,iFrequencyStepRefinement);
%disp(['length(X_FT_refined) :: ',num2str(length(X_FT_refined))])
%IX_FT_refined = fInverseFourierTransform(X_FT_refined,length(x))+x_mean;
%disp(['length(IX_FT_refined) :: ',num2str(length(IX_FT_refined))])
%disp([num2str(rms(x-fInverseFourierTransform(fFourierTransform(x,iFrequencyStepRefinement),length(x)))),' :: rms(x-fInverseFourierTransform(fFourierTransform(x,iFrequencyStepRefinement),length(x)))'])
%%length(X_FT_refined) :: 2000
%%length(IX_FT_refined) :: 1000
%%1.8981e-14 :: rms(x-fInverseFourierTransform(fFourierTransform(x,iFrequencyStepRefinement),length(x)))
%
%iFrequencyStepRefinement = 5;
%X_FT_refined = fFourierTransform(x-x_mean,iFrequencyStepRefinement);
%disp(['length(X_FT_refined) :: ',num2str(length(X_FT_refined))])
%IX_FT_refined = fInverseFourierTransform(X_FT_refined,length(x))+x_mean;
%disp(['length(IX_FT_refined) :: ',num2str(length(IX_FT_refined))])
%disp([num2str(rms(x-fInverseFourierTransform(fFourierTransform(x,iFrequencyStepRefinement),length(x)))),' :: rms(x-fInverseFourierTransform(fFourierTransform(x,iFrequencyStepRefinement),length(x)))'])
%%length(X_FT_refined) :: 5000
%%length(IX_FT_refined) :: 1000
%%1.1458e-14 :: rms(x-fInverseFourierTransform(fFourierTransform(x,iFrequencyStepRefinement),length(x)))
%
%iFrequencyStepRefinement = 10;
%X_FT_refined = fFourierTransform(x-x_mean,iFrequencyStepRefinement);
%disp(['length(X_FT_refined) :: ',num2str(length(X_FT_refined))])
%IX_FT_refined = fInverseFourierTransform(X_FT_refined,length(x))+x_mean;
%disp(['length(IX_FT_refined) :: ',num2str(length(IX_FT_refined))])
%disp([num2str(rms(x-fInverseFourierTransform(fFourierTransform(x,iFrequencyStepRefinement),length(x)))),' :: rms(x-fInverseFourierTransform(fFourierTransform(x,iFrequencyStepRefinement),length(x)))'])
%%length(X_FT_refined) :: 10000
%%length(IX_FT_refined) :: 1000
%%8.5345e-15 :: rms(x-fInverseFourierTransform(fFourierTransform(x,iFrequencyStepRefinement),length(x)))
%
%iFrequencyStepRefinement = 100;
%X_FT_refined = fFourierTransform(x-x_mean,iFrequencyStepRefinement);
%disp(['length(X_FT_refined) :: ',num2str(length(X_FT_refined))])
%IX_FT_refined = fInverseFourierTransform(X_FT_refined,length(x))+x_mean;
%disp(['length(IX_FT_refined) :: ',num2str(length(IX_FT_refined))])
%disp([num2str(rms(x-fInverseFourierTransform(fFourierTransform(x,iFrequencyStepRefinement),length(x)))),' :: rms(x-fInverseFourierTransform(fFourierTransform(x,iFrequencyStepRefinement),length(x)))'])
%%length(X_FT_refined) :: 100000
%%length(IX_FT_refined) :: 1000
%%5.8242e-15 :: rms(x-fInverseFourierTransform(fFourierTransform(x,iFrequencyStepRefinement),length(x)))
%
%
%%% ODD samples
%
%t=(0:0.01:10)'; disp(['t=(0:0.01:10)'''])
%x = 0.0 + 0.3*sin(2*pi*0.5*(t- 0.7)) + 0.7*cos(2*pi*0.3*(t- 0.5)); disp(['x = 0.3*sin(2*pi*0.5*(t- 0.7)) + 0.7*cos(2*pi*0.3*(t- 0.5))'])
%x_mean = mean(x);
%disp(['length(x) :: ',num2str(length(x))])
%disp(['mean(x) :: ',num2str(x_mean)])
%disp(['rms(x) :: ',num2str(rms(x))])
%disp([num2str(rms(fft(x)-fFourierTransform(x))),' :: rms(fft(x)-fFourierTransform(x))'])
%%t=(0:0.01:10)'
%%x = 0.3*sin(2*pi*0.5*(t- 0.7)) + 0.7*cos(2*pi*0.3*(t- 0.5))
%%length(x) :: 1001
%%mean(x) :: 0.00016858
%%rms(x) :: 0.53827
%%6.8017e-13 :: rms(fft(x)-fFourierTransform(x))
%
%disp([num2str(rms(x-ifft(fft(x)))),' :: rms(x-ifft(fft(x)))'])
%%1.7457e-16 :: rms(x-ifft(fft(x)))
%
%iFrequencyStepRefinement = 1;
%X_FT = fFourierTransform(x-x_mean,iFrequencyStepRefinement);
%disp(['length(X_FT) :: ',num2str(length(X_FT))])
%N = length(x);
%IX_FT = fInverseFourierTransform(X_FT,N)+x_mean;
%disp(['length(IX_FT) :: ',num2str(length(IX_FT))])
%disp([num2str(rms(x-fInverseFourierTransform(fFourierTransform(x,iFrequencyStepRefinement),length(x)))),' :: rms(x-fInverseFourierTransform(fFourierTransform(x,iFrequencyStepRefinement),length(x)))'])
%%length(X_FT) :: 1001
%%length(IX_FT) :: 1001
%%2.1492e-14 :: rms(x-fInverseFourierTransform(fFourierTransform(x,iFrequencyStepRefinement),length(x)))
%
%%D=[x,IX_FT, x-IX_FT, IX_FT./x, [0;diff(IX_FT./x)]];
%
%iFrequencyStepRefinement = 2;
%X_FT_refined = fFourierTransform(x-x_mean,iFrequencyStepRefinement);
%disp(['length(X_FT_refined) :: ',num2str(length(X_FT_refined))])
%IX_FT_refined = fInverseFourierTransform(X_FT_refined,length(x))+x_mean;
%disp(['length(IX_FT_refined) :: ',num2str(length(IX_FT_refined))])
%disp([num2str(rms(x-fInverseFourierTransform(fFourierTransform(x,iFrequencyStepRefinement),length(x)))),' :: rms(x-fInverseFourierTransform(fFourierTransform(x,iFrequencyStepRefinement),length(x)))'])
%%length(X_FT_refined) :: 2002
%%length(IX_FT_refined) :: 1001
%%1.6457e-14 :: rms(x-fInverseFourierTransform(fFourierTransform(x,iFrequencyStepRefinement),length(x)))
%
%iFrequencyStepRefinement = 5;
%X_FT_refined = fFourierTransform(x-x_mean,iFrequencyStepRefinement);
%disp(['length(X_FT_refined) :: ',num2str(length(X_FT_refined))])
%IX_FT_refined = fInverseFourierTransform(X_FT_refined,length(x))+x_mean;
%disp(['length(IX_FT_refined) :: ',num2str(length(IX_FT_refined))])
%disp([num2str(rms(x-fInverseFourierTransform(fFourierTransform(x,iFrequencyStepRefinement),length(x)))),' :: rms(x-fInverseFourierTransform(fFourierTransform(x,iFrequencyStepRefinement),length(x)))'])
%%length(X_FT_refined) :: 5005
%%length(IX_FT_refined) :: 1001
%%1.142e-14 :: rms(x-fInverseFourierTransform(fFourierTransform(x,iFrequencyStepRefinement),length(x)))
%
%iFrequencyStepRefinement = 10;
%X_FT_refined = fFourierTransform(x-x_mean,iFrequencyStepRefinement);
%disp(['length(X_FT_refined) :: ',num2str(length(X_FT_refined))])
%IX_FT_refined = fInverseFourierTransform(X_FT_refined,length(x))+x_mean;
%disp(['length(IX_FT_refined) :: ',num2str(length(IX_FT_refined))])
%disp([num2str(rms(x-fInverseFourierTransform(fFourierTransform(x,iFrequencyStepRefinement),length(x)))),' :: rms(x-fInverseFourierTransform(fFourierTransform(x,iFrequencyStepRefinement),length(x)))'])
%%length(X_FT_refined) :: 10010
%%length(IX_FT_refined) :: 1001
%%7.9148e-15 :: rms(x-fInverseFourierTransform(fFourierTransform(x,iFrequencyStepRefinement),length(x)))
%
%iFrequencyStepRefinement = 100;
%X_FT_refined = fFourierTransform(x-x_mean,iFrequencyStepRefinement);
%disp(['length(X_FT_refined) :: ',num2str(length(X_FT_refined))])
%IX_FT_refined = fInverseFourierTransform(X_FT_refined,length(x))+x_mean;
%disp(['length(IX_FT_refined) :: ',num2str(length(IX_FT_refined))])
%disp([num2str(rms(x-fInverseFourierTransform(fFourierTransform(x,iFrequencyStepRefinement),length(x)))),' :: rms(x-fInverseFourierTransform(fFourierTransform(x,iFrequencyStepRefinement),length(x)))'])
%%length(X_FT_refined) :: 100100
%%length(IX_FT_refined) :: 1001
%%5.9372e-15 :: rms(x-fInverseFourierTransform(fFourierTransform(x,iFrequencyStepRefinement),length(x)))
%
%
%disp(' ')
%disp([num2str(rms(x-ifft(fft(x)))),' :: rms(x-ifft(fft(x)))'])
%disp([num2str(rms(ifft(fft(x))-fInverseFourierTransform(fFourierTransform(x)))), ' :: rms(ifft(fft(x))-fInverseFourierTransform(fFourierTransform(x))) :: '])
%disp([num2str(rms(x-fInverseFourierTransform(fFourierTransform(x)))),' :: rms(x-fInverseFourierTransform(fFourierTransform(x)))'])
%disp([num2str(rms(fft(x)-fFourierTransform(x))),' :: rms(fft(x)-fFourierTransform(x))'])
%disp([num2str(rms(x-ifft(fFourierTransform(x)))),' :: rms(x-ifft(fFourierTransform(x)))'])
%disp([num2str(rms(x-fInverseFourierTransform(fft(x)))),' :: rms(x-fInverseFourierTransform(fft(x)))'])
%%1.7457e-16 :: rms(x-ifft(fft(x)))
%%2.1487e-14 :: rms(ifft(fft(x))-fInverseFourierTransform(fFourierTransform(x))) :: 
%%2.1492e-14 :: rms(x-fInverseFourierTransform(fFourierTransform(x)))
%%6.8017e-13 :: rms(fft(x)-fFourierTransform(x))
%%2.1492e-14 :: rms(x-ifft(fFourierTransform(x)))
%%6.6909e-16 :: rms(x-fInverseFourierTransform(fft(x)))
%
%iFrequencyStepRefinement = 2;
%disp(' ')
%disp([num2str(rms(x-ifft(fft(x)))),' :: rms(x-ifft(fft(x)))'])
%disp([num2str(rms(ifft(fft(x))-fInverseFourierTransform(fFourierTransform(x,iFrequencyStepRefinement),length(x)))),' :: rms(ifft(fft(x))-fInverseFourierTransform(fFourierTransform(x,iFrequencyStepRefinement),length(x)))'])
%disp([num2str(rms(x-fInverseFourierTransform(fFourierTransform(x,iFrequencyStepRefinement),length(x)))),' :: rms(x-fInverseFourierTransform(fFourierTransform(x,iFrequencyStepRefinement),length(x)))'])
%%1.7457e-16 :: rms(x-ifft(fft(x)))
%%1.645e-14 :: rms(ifft(fft(x))-fInverseFourierTransform(fFourierTransform(x,iFrequencyStepRefinement),length(x)))
%%1.6457e-14 :: rms(x-fInverseFourierTransform(fFourierTransform(x,iFrequencyStepRefinement),length(x)))
%

%%DEV

%L=11;w=[(0:(L/2)),(floor(-L/2+1):(-1))],length(w),wmax=w(floor(L/2)+1),wmax_is_Nyquist=(mod(L,2)==0)
%%w = 0     1     2     3     4     5    -5    -4    -3    -2    -1
%%11, 5, false
%L=11;w=[(0:(L/2)),(floor(-L/2+1):(-1))],length(w),wmax=w(floor(L/2)+1),wmax_is_Nyquist=(mod(L,2)==0)
%%w = 0     1     2     3     4     5     6    -5    -4    -3    -2    -1
%%12, 6, true

%clr
%Ts = 1/7;			% Sampling period %%this is the key/obstacle for a precise dominant frequency identification!
%Fs = 1/Ts;			% Sampling frequency
%Ns = 1024+512;			% Number of signal samples
%t = (0:Ns-1)'*Ts;	% Time column vector
%
%amp  = [11, 7, 5]			%amplitude
%peri = [64, 129, 255]		%period
%freq = 1./peri				%frequency = 1/period
%fish = [32, 65, 127]		%phase shift
%
%x1 = amp(1)*cos(2*pi*freq(1)*(t-fish(1)));	% First column signal
%x2 = amp(2)*cos(2*pi*freq(2)*(t-fish(2)));	% Second column signal
%x3 = amp(3)*cos(2*pi*freq(3)*(t-fish(3)));	% Third column signal
%
%x = [x1, x2, x3]; %joined signals by columns
%x_mean = mean(x);
%
%fMyFigure(1);
%for iCol = 1:3
%    subplot(3,1,iCol)
%    plot(t,x(:,iCol))
%    title(['Column ',num2str(iCol),' in the Time Domain'])
%end
%
%%For algorithm performance purposes, fft allows you to pad the input with trailing zeros. 
%% pad each column of X with zeros so that the length of each row is the next higher power of 2 from the current length. 
%% Define the new length using the nextpow2 function.
%n = 2^nextpow2(Ns);
%
%disp(' ')
%disp('------------')
%disp('matlab.fft()')
%disp('------------')
%X = fft(x-x_mean,n);
%
%F = (0:(Fs/n):(Fs/2))'; %frequency domain
%iF = (1:n/2+1)'; %frequency index
%
%%Calculate the double-sided spectrum and single-sided spectrum of each signal.
%P2 = abs(X/Ns); %% no matter the padding, the amplitude still depends on the original sample size
%P1 = P2(1:n/2+1,:); %% analysed frequency spectrum depends on the padded sample size
%%% first is the DC (average) component; it is 0 when using as fft(x-mean(x))
%%% for even n (and preferably n is 2^k, so even) the n/2+1 is the Nyquist frequency
%%% the amplitude of these two components is NOT multiplied by *2
%P1(2:end-1,:) = 2*P1(2:end-1,:); 
%
%FI = angle(X(iF)); %phase shift
%
%F_max=[Fs/2, F(iF(n/2+1)), 1] %actual Nyquist, analysed last freq, is it the Nyquist
%[amp_MAX, iMAX] = max(P1)
%freq_MAX = [F(iMAX(1)), F(iMAX(2)), F(iMAX(3))]
%phase_MAX = [FI(iMAX(1)), FI(iMAX(2)), FI(iMAX(3))]
%Period_MAX = 1./freq_MAX
%PhaseShift_MAX = mod(-phase_MAX./(2*pi*freq_MAX),Period_MAX)
%
%fMyFigure(2);
%for iCol=1:3
%    subplot(3,1,iCol)
%    plot(F,P1(iF,iCol))
%    title(['Column ',num2str(iCol),' in the Frequency Domain'])
%end
%
%x_restored = ifft(X)+x_mean;
%
%fMyFigure(3);
%for iCol = 1:3
%    subplot(3,1,iCol)
%    plot(t,x_restored(1:Ns,iCol))
%    title(['Restored Column ',num2str(iCol),' in the Time Domain'])
%end
%
%fMyFigure(4);
%for iCol = 1:3
%    subplot(3,1,iCol)
%    plot(t,x(:,iCol)-x_restored(1:Ns,iCol))
%    title(['Original - Restored Column ',num2str(iCol),' ERROR in the Time Domain'])
%end
%
%disp(' ')
%disp('-----------------')
%disp('fFourierTransform')
%disp('-----------------')
%iRefined = 256
%[X_FT, fTs] = fFourierTransform(x-x_mean,iRefined);
%%the analysed frequency domain includes the mean component (0 frequency) and all only positive frequencies, the last is the Nyquist frequency
%F_FT = fTs./Ts; 
%%F_FT_max=[Fs/2, F_FT(end)] %actual Nyquist, analysed last freq
%
%%Calculate the spectrum
%P_FT = 2*abs(X_FT/Ns); %% no mater the refinement, the amplitude still depends on the original sample size
%
%%Calculate the phase shift
%FI_FT = angle(X_FT);
%
%%find amplitude of max power
%[amp_FT_MAX, i_FTMAX] = max(P_FT);amp_FT_MAX
%%present reference value
%amp
%
%%find frequency of max power
%freq_FT_MAX = [F_FT(i_FTMAX(1)), F_FT(i_FTMAX(2)), F_FT(i_FTMAX(3))]
%%present reference value
%freq
%
%%find phase of max power
%phase_FT_MAX = [FI_FT(i_FTMAX(1)), FI_FT(i_FTMAX(2)), FI_FT(i_FTMAX(3))];
%
%%period of max power
%Period_FT_MAX = 1./freq_FT_MAX
%%present reference value
%peri
%
%%phase shift of max power
%PhaseShift_FT_MAX = mod(-phase_FT_MAX./(2*pi*freq_FT_MAX),Period_FT_MAX)
%%present reference value
%fish
%
%fMyFigure(5);
%for iCol=1:3
%    subplot(3,1,iCol)
%    plot(1./F_FT,P_FT(:,iCol))
%    title(['Column_FT ',num2str(iCol),' in the Frequency Domain'])
%end
%
%x_FT_restored = fInverseFourierTransform(X_FT,length(t))+x_mean;
%
%fMyFigure(6);
%for iCol = 1:3
%    subplot(3,1,iCol)
%    plot(t,x_FT_restored(:,iCol))
%    title(['Restored Column_FT ',num2str(iCol),' in the Time Domain'])
%end
%
%fMyFigure(7);
%for iCol = 1:3
%    subplot(3,1,iCol)
%    plot(t,x(:,iCol)-x_FT_restored(:,iCol))
%    title(['Original - Restored Column_FT ',num2str(iCol),' ERROR in the Time Domain'])
%end