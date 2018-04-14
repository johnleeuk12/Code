function ys=gaussian_kern_lin_reg(xs,x,y,h)

% Gaussian kernel function
K=sqdist(diag(1./h)*x,diag(1./h)*xs);
K=exp(-K/2);

% linear kernel regression
ys=sum(K'.*y)/sum(K);
