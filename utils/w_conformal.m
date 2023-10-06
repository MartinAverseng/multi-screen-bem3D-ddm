function [r] = w_conformal(z)

root1 = sqrt(z.^2 - 1);
root2 = -root1;

root = root1;
root(imag(root1).*imag(z) < -1e-10) = root2(imag(root1).*imag(z) < -1e-10);
root(imag(z) == 0) = root1(imag(z)==0);

r = z + root;

dr = 1 + z./root; 


end

