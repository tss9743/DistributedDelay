function [lowerBound,upperBound,k] = gammaKernel(a,p)
% Creating the shape of the delay dist to be a gamma kernel
global r N
tj = [0:N]*(r/N);

%k = a^b.*tj.^(b-1).*exp(-a.*tj)./factorial(b-1);
k = a.*(a.*tj).^(p-1)./gamma(p).*exp(-a.*tj);

% norm = 1/sum(k);
% k = k*norm;

% ensuring that the truncated Gamma dist is a true dist
mean1 = p/a;
sigma1 = sqrt(p/a^2);
% find 99.7% of Prob Dist with +-3sigma from mean
threeSigma = find(tj <= mean1+3*sigma1 & tj >= mean1-3*sigma1);
notThreeSigma = find(tj > mean1+3*sigma1 | tj < mean1-3*sigma1);

k(notThreeSigma) = 0;

% scale remaining values so that the area under is 1
numNotZero = length(threeSigma);
% find min and max index in the 3sigma range
lowerBound = min(threeSigma);
upperBound = max(threeSigma);
% do trapezoid rule for integral approx
if numNotZero > 2
    isum = k(lowerBound) + k(upperBound);
    for i = lowerBound+1:upperBound-1
        isum = isum + k(i)*2;
    end
    trapInt = ((tj(upperBound)-tj(lowerBound))/(2*(numNotZero-1)))*isum;
else
    trapInt = (tj(upperBound)-tj(lowerBound))*((k(lowerBound)+k(upperBound))/2);
end
% find the weight by which to scale the gamma function by so the area is 1
c = 1/trapInt;
k = k*c;

k = fliplr(k);
end