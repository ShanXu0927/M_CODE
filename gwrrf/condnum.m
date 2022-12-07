function ci = condnum(x,w)
    w = sqrt(w);
    sw = sum(w);
    w = bsxfun(@rdivide,w,sw);
    xw = bsxfun(@times,x,w);
    sxw = sqrt(sum(xw.^2));
    sxw = bsxfun(@rdivide,xw,sxw);
    svdx = svd(sxw);
    ci = svdx(1)./svdx(end);
%     la = (svdx(1) - 30*svdx(end))/(30 - 1);
end