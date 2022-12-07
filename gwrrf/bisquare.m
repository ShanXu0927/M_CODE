function w = bisquare(dist,bandwidth)

% w = (dist<bandwidth).*(1-(dist/bandwidth).^2).^2);
w = bsxfun(@le,dist,bandwidth).*...
    (1-(bsxfun(@rdivide,dist,bandwidth)).^2).^2;

end