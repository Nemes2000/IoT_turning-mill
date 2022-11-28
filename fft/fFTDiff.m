function dx_dt = fFTDiff(x,t,iOrderN)
%function dx_dt = fFTDiff(x,t,iOrderN)
% Spectral (Fourier Transform based) Diff of iOrderN
% works ok only for continuous signals; NOT good for discontinuous random and integer sequences
% even a single discontinuity in the first derivative of the signal will cause osculations aka. "Gibbs phenomenon"
% x = column vector of samples, 
% t = OPTIONAL column vector of sample time, DEFULT = (1:iSampleN)'
% iOrderN = OPTIONAL scalar, DEFAULT = 1
%	dx./dt =  ifft(i*w*fft(x(t))); where fft(t) ~ x(t)*exp(i*w)
%
% see "MIT_Notes on FFT-based differentiation.(fft-deriv).pdf"
%

	iSampleN = size(x,1);
	if(~exist('t') || isempty(t) || ~all(isfinite(t))), t=(1:iSampleN)'; end
	if(~exist('iOrderN') || isempty(iOrderN) || ~all(isfinite(iOrderN))), iOrderN=1; end

	dt = (t(:,1) - [0;t(1:(end-1),1)]); %dt

	%x_fft = fft(x(:,1));
	x_fft = fNonuniformFourierTransform(x(:,1));

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

	w_fft = i*(2*pi)*fTs./dt;
	%see https://math.mit.edu/~stevenj/fft-deriv.pdf
	if((mod(iSampleN,2)==0) & (mod(iOrderN,2)==1)) %if(even sample number AND odd diff order), the X(N/2) has to vanish because sin(n*pi) = 0
		w_fft(iSampleN/2+1) = 0; %so w_fft must be precisely 0 for k=N/2, k starts from 0, index is (N/2+1)
	end
	w_fft = w_fft.^iOrderN;
	
	dx_dt_fft = w_fft.*x_fft;

	%dx_dt = real(ifft(dx_dt_fft));
	dx_dt = real(fInverseFourierTransform(dx_dt_fft));

return

%			maxE		meanE		stdE
%fft/ifft: 1.7319e-14,  3.5760e-17, 5.7844e-15
%FT/IFT:   4.7073e-13, -1.0010e-17, 1.9626e-13
%

%%%TEST
%%pSpectralDerivative
%%
%%rehash path

%clr
%iSampleN=128; % 128 / 127
%
%t = (1:iSampleN)';
%w = (2*pi/iSampleN);
%x = w*t;
%dx = w;
%y = sin(x);
%dy = cos(x)*dx;
%Iy = -cos(x)/dx;
%x0 = 0;
%y0 = sin(x0);
%
%%% SpectralDiff is superb for continuous periodic functions
%dy_dt = fFTDiff(y,t);
%max(max(abs(dy_dt-dy))) %2.3169e-14 / 3.9010e-14
%max(mean(dy_dt-dy)) %-7.5620e-19 / -5.1905e-19
%max(std(dy_dt-dy)) %9.6332e-15 / 1.7341e-14
%
%%% whereas the simple num diff performs quite poorly...
%dy_dt_num = fNumDiff(y,t,y0,x0);
%max(max(abs(dy_dt_num-dy))) %0.0012 / 0.0012
%max(mean(dy_dt_num-dy)) %-1.0842e-19 / -1.2567e-18
%max(std(dy_dt_num-dy)) %8.5520e-04 / 8.6875e-04
%
%%% SpectralInt is superb for continuous periodic functions
%Iy_FFT = fFTInt(y,t,y0);
%max(max(abs(Iy-Iy_FFT))) %1.9540e-14
%max(mean(Iy-Iy_FFT)) %-2.8105e-16
%max(std(Iy-Iy_FFT)) %-2.8105e-16
%%
%% FFTdif and FFTint are inverses of each other
%Idy_dt_FFT = fFTInt(dy_dt,t,y0);
%max(max(abs(y-Idy_dt_FFT))) %1.5488e-14
%dIy_FFT_dt = fFTDiff(Iy_FFT,t);
%max(max(abs(y-dIy_FFT_dt))) %8.0347e-13

%clr
%
%%% SpectralDiff (FourierTransform based) is no good for simple integer and discontinuous random series
%iSampleN=128; %127
%iFeatureN=1;
%a=((1:iSampleN)').^2;
%da_n=fNumDiff(a,[],NaN,[]);
%t = (1:iSampleN)';
%%
%da=fFTDiff(a,t,1);
%any(a - (da_n+[0;a(1:end-1)])) %==zeros
%any(a - (da+[0;a(1:end-1)])) %NOT ==zeros !!!
%max(max(abs(a - (da+[0;a(1:end-1)])))) % 1.1483e+04 !!!
%%
%x = repmat(t,[1,iFeatureN]);
%dt = 1*ones(1,iFeatureN);
%x0 = zeros(1,iFeatureN);
%dx_dt_n = fNumDiff(x,[],x0,[]);
%max(max(abs(dx_dt_n-ones(iSampleN,iFeatureN)))) %0
%%
%warning('SpectralDiff is NOT good for functions with discontinuities, especially not for point-by-point discreet values!')
%dx_dt = fFTDiff(x,t);
%max(max(abs(dx_dt-ones(iSampleN,iFeatureN)))) %88.7164 !!!
%%
%%
%x0 = zeros(1,iFeatureN);
%x3 = rand(iSampleN,iFeatureN);
%dt3 = 0.001*ones(1,iFeatureN);
%%
%warning('SpectralDiff is NOT good for functions with discontinuities, especially not for point-by-point random values!')
%dx3_dt = fFTDiff(x3,0.001*x);
%%(Idx3_dt = fNumInt(dx3_dt,dt3,x0);)
%Idx3_dt = fFTInt(dx3_dt,0.001*x,x0);
%max(max(abs(x3-Idx3_dt))) %0.5599 !!! (0.9621)
%%
%%(Ix3 = fNumInt(x3,dt3,x0);)
%fftIx3 = fFTInt(x3,0.001*x,x0);
%dIx3_dt = fFTDiff(fftIx3,0.001*x);
%max(max(abs(x3-dIx3_dt))) %0.5599 !!! (45.2497)
%%
%disp('numInt is the inverse of numDif')
%dx3_dt_num = fNumDiff(x3,0.001*x,x0,[]);
%Idx3_dt_num = fNumInt(dx3_dt_num,dt3,x0);
%max(max(abs(x3-Idx3_dt_num))) %1.0103e-14
%%
%Ix3 = fNumInt(x3,dt3,x0);
%dIx3_dt_num = fNumDiff(Ix3,0.001*x,x0,[]);
%max(max(abs(x3-dIx3_dt_num))) %7.3275e-15
