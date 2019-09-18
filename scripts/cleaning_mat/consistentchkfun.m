    function [nchi1,nchi2]= consistentchkfun(chi1,chi2,lb,ub);
    ratio = chi1./chi2;
    bad = find(ratio <= lb | ratio >= ub);
    nchi1 = chi1; nchi2 = chi2;
    for i = 1:length(bad);
        if ~isnan(nchi1(bad(i)) ); 
           nchi1(bad(i)) = min([chi1(bad(i)) chi2(bad(i))]);
        end
        if ~isnan(nchi2(bad(i)) ); 
            nchi2(bad(i)) = min([chi1(bad(i)) chi2(bad(i))]);
        end
    end
