function [coefblm, coefclm] = gaussian_vshcoef(L, c1, c2, b0, c0, r, power)
%L = bw c1 for blm coeff c2 for clm coeff r correlation
%b0, c0 for variance at l=0

    B = L + 1;
    Len = (L+1)*(L+2)/2;
    
    mu1 = zeros(2*Len, 1);
    var0 = zeros(Len, 1);
    
    id = 2; %make id = 1 empty to put b0 or c0
    r1 = 2;

    while id < Len
        for l = 1:L
            tmp = 1/l^power;
            r2 = r1 + l;
            var0(id) = 2*tmp; % corresponds to m = 0, for the current l
            id = id + 1;
            for m = (r1+1):r2
                var0(id) = tmp;
                id = id + 1;
            end
            l = l+1;
            r1 = r2+1;
        end
    end

    var1 = c1*var0;
    var2 = c2*var0;
    var1(1) = b0;
    var2(1) = c0;
    cov11 = diag(var1); 
    cov12 = diag(var2);
    cov00 = r * sqrt(cov11) * sqrt(cov12);
    if r == 0
        cov00 = zeros(Len, Len);
    end
    
    covtmp1 = horzcat(cov11, cov00);
    covtmp2 = horzcat(cov00', cov12);
    cov1 = vertcat(covtmp1, covtmp2);
    almR = mvnrnd(mu1, cov1);
    almI = mvnrnd(mu1, cov1);
    almR1 = almR(1:Len); almI1 = almI(1:Len);
    almR2 = almR(Len+1: 2*Len); almI2 = almI(Len+1:2*Len);

    coefblm = zeros([B 2*B-1]);
    coefclm = zeros([B 2*B-1]);    
    %note the l, m= 0 only real value
    id = 1;
    while id < Len
        for l = 0:L
            for m = 0:l
                coefblm(l+1, B+m) = almR1(id) + 1i*almI1(id);
                coefclm(l+1, B+m) = almR2(id) + 1i*almI2(id);
                coefblm(l+1, B-m) = (-1)^m*(almR1(id) - 1i*almI1(id));
                coefclm(l+1, B-m) = (-1)^m*(almR2(id) - 1i*almI2(id)); 
                id = id +1;
            end
        end
    end
    coefblm(:,B) = real(coefblm(:,B));
    coefclm(:,B) = real(coefclm(:,B));
end