function w = guass(dist,bandwidth)
% 
% w = exp(-(dist.^2)./(bandwidth.^2));
% w = exp(-0.5*((dist).^2)./(bandwidth.^2));
w = exp(-0.5*(bsxfun(@rdivide,3*dist,bandwidth)).^2);
end