function Ix = fFTInt(x,t,Ix0,iOrderN)
%function Ix = fFTInt(x,t,Ix0,iOrderN)
% Spectral (Fourier Transform based) Integral
% works ok only for continuous signals; NOT good for discontinuous random and integer sequences
% even a single discontinuity in the first derivative of the signal will cause osculations aka. "Gibbs phenomenon"
% x = column vector of samples, 
% t = OPTIONAL column vector of sample time, DEFULT = (1:iSampleN)'
% Ix0 optional constant for integral
% iOrderN = OPTIONAL scalar, DEFAULT = 1
% Integral(exp(k*x)dx) = (1/k)*exp(k*x) + Ix0
%
% see "MIT_Notes on FFT-based differentiation.(fft-deriv).pdf"
%

	iSampleN = size(x,1);
	if(~exist('t') || isempty(t) || ~all(isfinite(t))), t=(1:iSampleN)'; end
	if(~exist('Ix0') || isempty(Ix0) || ~isfinite(Ix0)), Ix0=0; end
	if(~exist('iOrderN') || isempty(iOrderN) || ~all(isfinite(iOrderN))), iOrderN=1; end

	dt = (t(:,1) - [0;t(1:(end-1),1)]); %dt

	%x_ftt = fft(x(:,1));
	x_ftt = fNonuniformFourierTransform(x(:,1));

	%%w_fft = (2*pi)*(-1/2:1/iSampleN:(1/2-1/iSampleN))'./dt;
	%%w_fft = fftshift(w_fft); %[-3 -2 -1 0 1 2] -> [0 1 2 -3 -2 -1]; [-3 -2 -1 0 1 2 3] -> [1 2 3 -3 -2 -1 0]
	%w_fft = (2*pi)*(-iSampleN/2:1:(iSampleN/2-1))'./(dt*iSampleN);
	%w_fft = fftshift(w_fft); %[-3 -2 -1 0 1 2] -> [0 1 2 -3 -2 -1]; [-3 -2 -1 0 1 2 3] -> [1 2 3 -3 -2 -1 0]
	%
	%L = iSampleN
	%fTs = [(0:(L/2)),(floor(-L/2+1):(-1))].'/L; %this is not frequency, but f*Ts; in "fftshift" order
	%w_fft = (2*pi)*fTs/dt;
	%
	%L=11;n_df=[(0:(L/2)),(floor(-L/2+1):(-1))],length(n_df),n_df_max=n_df(floor(L/2)+1),wmax_is_Nyquist=(mod(L,2)==0)
	%	%n_df = 0     1     2     3     4     5    -5    -4    -3    -2    -1
	%	%11, 5, false
	%L=11;n_df=[(0:(L/2)),(floor(-L/2+1):(-1))],length(n_df),n_df_max=n_df(floor(L/2)+1),wmax_is_Nyquist=(mod(L,2)==0)
	%	%n_df = 0     1     2     3     4     5     6    -5    -4    -3    -2    -1
	%	%12, 6, true
	fTs = [(0:(iSampleN/2)),(floor(-iSampleN/2+1):(-1))].'/iSampleN; %this is not frequency, but f*Ts; in "fftshift" order

	rec_w_fft = dt./(i*(2*pi)*fTs);
	rec_w_fft(1) = 0; %fTs(1)==0;
	%see https://math.mit.edu/~stevenj/fft-deriv.pdf
	if((mod(iSampleN,2)==0) && (mod(iOrderN,2)==1)) %if(even sample number AND odd diff order), the X(N/2) has to vanish because sin(n*pi) = 0
		rec_w_fft(iSampleN/2+1) = 0; %so w_fft must be precisely 0 for k=N/2, k starts from 0, index is (N/2+1)
	end
	rec_w_fft = rec_w_fft.^iOrderN;
	%w_fft(find(abs(w_fft)==0)) = eps; %w_fft(1)==0; so 1/w_fft would be not defined (Inf in Matlab, thus Ix would be Inf also)

	Ix_fft = rec_w_fft.*x_ftt;

	%Ix = real(ifft(Ix_fft)) +Ix0;
	Ix = real(fInverseFourierTransform(Ix_fft)) +Ix0;

return

%%TEST
%pSpectralIntegral

%% see fFTDiff.m

%clr
%%iSampleN=127;
%iSampleN=128;
%x = (1:iSampleN)';
%w = (2*pi/iSampleN);
%t = w*x;
%dt2 = w;
%x2 = sin(t);
%t0 = 0;
%x0 = sin(t0);
%
%%% SpectralDiff is superb for reversing the real integral
%Ix2_rel = -cos(t);
%dIx2_dt_rel = fFTDiff(Ix2_rel,t);
%max(max(abs(x2-dIx2_dt_rel))) %7.3699e-13
%
%%% SpectralInt is superb for reversing the real differential
%dx2_rel = cos(t);
%Idx2_dt_rel = fFTInt(dx2_rel,t,x0);
%max(max(abs(x2-Idx2_dt_rel))) %1.9984e-15
%