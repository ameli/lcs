function TestGaussKernel(KernelSize)

StandardDeviation = KernelSize / 3;
dim = 1;
Tolerance = 2;

n = (-KernelSize:KernelSize);
yn = GaussSecondDerivativeKernel(StandardDeviation,n,dim)
stem(n,yn)
grid on
hold on

%ylim([0,1.5*max(y)])
xlim([n(1)-Tolerance n(end)+Tolerance])

r = linspace (-KernelSize-Tolerance,KernelSize+Tolerance,1000);
y = GaussKernel(StandardDeviation,r,dim);
ddy = diff(y,2)./((r(2)-r(1)).^2);
plot(r(2:end-1),ddy);

% Gauss Kernel
function y = GaussKernel(sigma,r,dim)
y = (1/((sqrt(2*pi)*sigma)^dim))*exp(-(r.^2)./(2*sigma^2));

% Gauss First Derivative Kernel
function y = GaussFirstDerivativeKernel(sigma,r,dim)
y = (-r./(sigma^2)).*GaussKernel(sigma,r,dim);

% Gauss Second Derivative Kernel
function y = GaussSecondDerivativeKernel(sigma,r,dim)
y = ((r.^2 - sigma^2)./(sigma^4)).*GaussKernel(sigma,r,dim);