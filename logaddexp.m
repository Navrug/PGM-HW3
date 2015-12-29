function [ logsum ] = logaddexp( loga, logb )
   if (abs(loga - logb) >= 36.043653389117155)		% 2^52-1 = 4503599627370495.	log of that is 36.043653389117155867651465390794
		logsum = max(loga, logb);		 % this branch is necessary, to avoid logb - logb having too big value	
	else
		logsum = loga + log(1 + exp(logb - loga));
   end
end