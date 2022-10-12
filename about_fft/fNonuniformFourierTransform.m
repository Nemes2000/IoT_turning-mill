function [X, vFreqTs, vFreqTs_max, vFreqTs_max_isNyquist, iTransformedSampleN] = fNonuniformFourierTransform(x, iFrequencyStepRefinement, vTsFrequencyies)
%function [X, vFreqTs, vFreqTs_max, vFreqTs_max_isNyquist, iTransformedSampleN] = fNonuniformFourierTransform(x, optional iFrequencyStepRefinement=1, optional vTsFrequencyies=[])
%% fNonuniformFourierTransform is O(n*n)! opposed to fft's O(n*log(n))
% operates on multiple columns of x in parallel
% brakes down single matrix M = exp(i2pivFreqTs*n) operation to series of vector column operations
%	is a bit slower than fFourierTransform(), but copes with longer signals
% copes with non-uniform samples, when some are NaN as in missing data like x(isnan(x)) = 0;
% does NOT cope with condensed data of uneven sampling periods (use the greatest common denominator of sampling periods for the full span, and fill the missing data with NaN)
% NaN x values are numerically treated as being 0
% the analysed  frequency spectrum is f = (1/Ts) * vFreqTs; vFreqTs = [(0:(L/2)),(floor(-L/2+1):(-1))].'/L; where L=iFrequencyStepRefinement*N, and N=length(x)
%	vFreqTs is not frequency, but f*Ts; in "fftshift" order%
% NOTICE: must vFreqTs/Ts to get the actual Frequency values!
%	vFreqTs_max presents the max processed frequency <= Nyquist's frequency(*Ts)
%	vFreqTs_max_isNyquist is true for vFreqTs_max == Nyquist's frequency(*Ts)
% NOTICE: must have harmonic_amplitude = 2*abs(X)/(N*iFrequencyStepRefinement); %instead of fft(x) common abs(X)/N
%	do exclude from *2 the amplitude corresponding to the mean(DC) (vFreqTs(1)) and the Nyquist frequency (vFreqTs(fix(MyNsamples/2+1*(mod(MyNsamples,2)==0))), when (mod(MyNsamples,2)==0))
% vFreqTs = [(0:(L/2)),(floor(-L/2+1):(-1))].'/L; %this is not frequency, but f*Ts; in "fftshift" order
%	vFreqTs = n_df/L
%	L=11;n_df=[(0:(L/2)),(floor(-L/2+1):(-1))],length(n_df),n_df_max=n_df(floor(L/2)+1),wmax_is_Nyquist=(mod(L,2)==0)
%	%n_df = 0     1     2     3     4     5    -5    -4    -3    -2    -1
%	%11, 5, false
%	L=11;n_df=[(0:(L/2)),(floor(-L/2+1):(-1))],length(n_df),n_df_max=n_df(floor(L/2)+1),wmax_is_Nyquist=(mod(L,2)==0)
%	%n_df = 0     1     2     3     4     5     6    -5    -4    -3    -2    -1
%	%12, 6, true
% vFreqTs_max = vFreqTs(floor(L/2)+1);
% vFreqTs_max_isNyquist = (mod(L,2)==0);

%https://ccrma.stanford.edu/~jos/st/Matrix_Formulation_DFT.html
%https://blogs.uoregon.edu/seis/wiki/unpacking-the-matlab-fft/
	%N = length(x);
	%n = (0:(N-1));
	%k = (0:(N-1))'/N;
	%M = exp(-i*2*pi*k*n/N);
	%X = M*x;

	if(~exist('iFrequencyStepRefinement') || (~isempty(iFrequencyStepRefinement) &&  ~isfinite(iFrequencyStepRefinement))), iFrequencyStepRefinement=1; end
	if(~exist('vTsFrequencyies') || any(~isfinite(vTsFrequencyies)) || isempty(vTsFrequencyies)), vTsFrequencyies=[]; end

	[N,K] = size(x); %samples
	n = (0:(N-1)); %sample time in Ts; n=(t-min(t))./Ts
	L = floor(N*iFrequencyStepRefinement); %frequencies to check

	if(isempty(vTsFrequencyies))
%vFreqTs = (0:(N-1))'/N; %this is not frequency, but f*Ts
%vFreqTs = (0:1/iFrequencyStepRefinement:(N-1/iFrequencyStepRefinement))'/N;
%vFreqTs = (0:(L-1))'/L; %this is not frequency, but f*Ts; more simple format
%vFreqTs = (-L/2:(L/2-1))'/L; %this is not frequency, but f*Ts; no meaning in going above the Nyquist frequency vFreqTs=1/2
%!!!FAILED INVERSE for vFreqTs = (0:L)'/(2*L); %this is not frequency, but f*Ts; in range [0:0.5]
		vFreqTs = [(0:(L/2)),(floor(-L/2+1):(-1))].'/L; %this is not frequency, but f*Ts; in "fftshift" order
		vFreqTs_max = vFreqTs(floor(L/2)+1); %the max analysed frequency component
		vFreqTs_max_isNyquist = (mod(L,2)==0); % is the max analysed frequency component the same as the Nyquist's frequency == 1/2
	else
		vFreqTs = vTsFrequencyies;
		vFreqTs_max = max(vFreqTs); %the max analysed frequency component
		vFreqTs_max_isNyquist = (vFreqTs_max==(1/2)); % is the max analysed frequency component the same as the Nyquist's frequency == 1/2
	end

	iTransformedSampleN = N-sum(isnan(x)); %count ignored NaNs for each column
	x(isnan(x)) = 0; %ignoring NaN samples

	%M = exp(i2pivFreqTs*n);
	%X = M*x;
	X = zeros(L,K, 'like',x);
	i2pivFreqTs = -i*2*pi*vFreqTs; %performance optimisation
	for iK = 1:K
		for iF = 1:L
			X(iF,iK) = exp(i2pivFreqTs(iF)*n)*x(:,iK);
		end %for iF
	end %for iK

return

%%TST
%clr
%t=(0:0.01:1-0.01)'; disp(['t=(0:0.01:1-0.01)''']) %#100
%%t=(0:0.01:1)'; disp(['t=(0:0.01:1)''']) %#101
%x = 0.3*sin(2*pi*0.5*(t- 0.7)) + 0.7*cos(2*pi*0.3*(t- 0.5)); disp(['x = 0.3*sin(2*pi*0.5*(t- 0.7)) + 0.7*cos(2*pi*0.3*(t- 0.5))'])
%disp(['rms(x) :: ',num2str(rms(x))])
%%rms(x) :: 0.51698
%
%disp(['max(max(abs(fFourierTransform(x)-fNonuniformFourierTransform(x)))) :: ',num2str(max(max(abs(fFourierTransform(x)-fNonuniformFourierTransform(x)))))]) 
%%max(max(abs(fFourierTransform(x)-fNonuniformFourierTransform(x)))) :: 5.0243e-15
%
%disp(['rms(fFourierTransform(x)-fNonuniformFourierTransform(x)) :: ',num2str(rms(fFourierTransform(x)-fNonuniformFourierTransform(x)))])
%%rms(fFourierTransform(x)-fNonuniformFourierTransform(x)) :: 1.0523e-15
%
%disp(['rms(x-fInverseFourierTransform(fFourierTransform(x))) :: ',num2str(rms(x-fInverseFourierTransform(fFourierTransform(x))))])
%%rms(x-fInverseFourierTransform(fFourierTransform(x))) :: 3.1606e-15
%
%disp(['rms(x-fInverseFourierTransform(fNonuniformFourierTransform(x))) :: ',num2str(rms(x-fInverseFourierTransform(fFourierTransform(x))))])
%%rms(x-fInverseFourierTransform(fNonuniformFourierTransform(x))) :: 3.1606e-15
%
%iFrequencyStepRefinement = 2;
%max(max(abs(fFourierTransform(x,iFrequencyStepRefinement)-fNonuniformFourierTransform(x,iFrequencyStepRefinement)))) 
%%   5.0243e-15
%
%disp(['rms(fFourierTransform(x,iFrequencyStepRefinement)-fNonuniformFourierTransform(x,iFrequencyStepRefinement)) :: ',num2str(rms(fFourierTransform(x,iFrequencyStepRefinement)-fNonuniformFourierTransform(x,iFrequencyStepRefinement)))])
%%rms(fFourierTransform(x,iFrequencyStepRefinement)-fNonuniformFourierTransform(x,iFrequencyStepRefinement)) :: 1.1203e-15
%
%disp(['fInverseFourierTransform(fNonuniformFourierTransform(x,iFrequencyStepRefinement),iFrequencyStepRefinement)-fInverseFourierTransform(fFourierTransform(x,iFrequencyStepRefinement),iFrequencyStepRefinement)) :: ',num2str(rms(fInverseFourierTransform(fNonuniformFourierTransform(x,iFrequencyStepRefinement),iFrequencyStepRefinement)-fInverseFourierTransform(fFourierTransform(x,iFrequencyStepRefinement),iFrequencyStepRefinement)))])
%%fInverseFourierTransform(fNonuniformFourierTransform(x,iFrequencyStepRefinement),iFrequencyStepRefinement)-fInverseFourierTransform(fFourierTransform(x,iFrequencyStepRefinement),iFrequencyStepRefinement)) :: 1.2413e-16
%