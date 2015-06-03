function ACF = acf(Series,nLags) 

nFFT =  2^(nextpow2(length(Series)) + 1);
F    =  fft(Series-mean(Series) , nFFT);
F    =  F .* conj(F);
ACF  =  ifft(F);
ACF  =  ACF(1:nLags);         % Retain non-negative lags.
ACF  =  ACF ./ ACF(1);     % Normalize.
ACF  =  real(ACF);