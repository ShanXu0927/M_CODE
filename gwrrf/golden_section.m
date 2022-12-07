function [opt_val,opt_score,vals,scores] = golden_section(a, c, delta, funs, tol, max_iter, int_score)
b = a + delta * abs(c-a);
d = c - delta * abs(c-a);
diff = 1.0e9;
iters  = 0;
vals = nan(max_iter*2,1);
scores = nan(max_iter*2,1);
icount = 0;
while abs(diff) > tol
    iters = iters + 1;
    if iters> max_iter
        break
    end
    
    if int_score
        b = round(b);
        d = round(d);
    end
    
    flag = vals==b;
    if any(flag==1)
        score_b = scores(flag);
    else
        icount = icount + 1;
        score_b = funs(b);
        vals(icount) = b;
        scores(icount) = score_b;
    end
    
    flag = vals==d;
    if any(flag==1)
        score_d = scores(flag);
    else
        icount = icount + 1;
        score_d = funs(d);
        vals(icount) = d;
        scores(icount) = score_d;
    end
    
    diff = (b-d)/b;
    
    if isinf(score_b) && isinf(score_d)
        if b>=d
            opt_val = b;
            opt_score = score_b;
            a = d;
            d = b;
            b = c - delta * abs(c-a);
        else
            opt_val = d;
            opt_score = score_d;
            a = b;
            b = d;
            d = c - delta * abs(c-a);
        end
    else
        if score_b <= score_d
            opt_val = b;
            opt_score = score_b;
            if b>=d
                a = d;
                d = b;
                b = c - delta * abs(c-a);
            else
                c = d;
                d = b;
                b = a + delta * abs(c-a);
            end
        else
            opt_val = d;
            opt_score = score_d;
            if b>=d
                c = b;
                b = d;
                d = a + delta * abs(c-a);
            else
                a = b;
                b = d;
                d = c - delta * abs(c-a);
            end
        end
    end
end
