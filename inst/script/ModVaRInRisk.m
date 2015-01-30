
function val = ModVaRInRisk
  
    % ModVaRInRisk  Calculate variance formulas for modified VaR, modified
    % ETL, ordinary VaR and ordinary ETL using the variance matrix of moments
    %   risk = ModVaRInRisk returns variance in the following order
    %          1. modified VaR    2. VaR   3. modified ETL   4. ETL
    
    clear;

    syms mu ss mu3 mu4 mu5 mu6 mu7 mu8; 
    syms sigma g1 k; 
    syms z dz alpha_sig C;

    r1 = [ss  mu3                mu4                                mu5                                       ];
    r2 = [mu3 mu4-ss^2           mu5-4*ss*mu3                       mu6-ss*mu4-4*mu3^2                        ];
    r3 = [mu4 mu5-4*ss*mu3       mu6-mu3^2-6*ss*mu4+9*ss^3          mu7-mu3*mu4-3*ss*mu5-4*mu3*mu4+12*ss^2*mu3];
    r4 = [mu5 mu6-ss*mu4-4*mu3^2 mu7-3*ss*mu5-5*mu3*mu4+12*ss^2*mu3 mu8-mu4^2-8*mu3*mu5+16*ss*mu3^2           ];
    X  = [r1;r2;r3;r4];

    S = sym('S%d%d',5);

    S(1,1) = sigma^2;
    for i = 2:5 
        if i == 2
            A = X(2,2); f = ss^(1/2);
            g = gradient(f,ss);
        elseif i == 3
            A = X([2 3],[2 3]); f = mu3/ss;
            g = gradient(f,[ss mu3]);
        elseif i == 4
            A = X([2 4],[2 4]); f = mu4/(ss)^(3/2) - 3*(ss)^(1/2); 
            g = gradient(f,[ss mu4]);
        elseif i == 5
            A = X([2 3],[2 3]); f = (mu3^2)/(ss^(5/2)); 
            g = gradient(f,[ss mu3]);
        end

        v1 = g.'*A*g;
        x = simplify(subs(v1,ss,sigma^2));
        x = simplify(subs(x,mu3,g1*sigma^3));
        x = simplify(subs(x,mu4,(k+3)*sigma^4));
        S(i,i) = collect(x);
    end

    for i = 1:4 
        for  j = (i+1):5 

            if and(i == 1,j == 2)
                A = X([1 2],[1 2]); f = mu + ss^(1/2); 
                g = gradient(f,[mu ss]);
            elseif and(i == 1, j == 3)
                A = X([1 2 3],[1 2 3]); f = mu + mu3/ss; 
                g = gradient(f,[mu ss mu3]);
            elseif and(i == 1, j == 4)
                A = X([1 2 4],[1 2 4]); f = mu + mu4/(ss)^(3/2) - 3*(ss)^(1/2);
                g = gradient(f,[mu ss mu4]);
            elseif and(i == 1, j == 5)
                A = X([1 2 3],[1 2 3]); f = mu + mu3^2/(ss)^(5/2); 
                g = gradient(f,[mu ss mu3]);
            elseif and(i == 2, j == 3)
                A = X([2 3],[2 3]); f = ss^(1/2) + mu3/ss; 
                g = gradient(f,[ss mu3]);
            elseif and(i == 2, j == 4)
                A = X([2 4],[2 4]); f = ss^(1/2) + mu4/(ss)^(3/2) - 3*(ss)^(1/2);
                g = gradient(f,[ss mu4]);   
            elseif and(i == 2, j == 5) 
                A = X([2 3],[2 3]); f = ss^(1/2) + (mu3^2)/(ss^(5/2));
                g = gradient(f,[ss mu3]);    
            elseif and(i == 3, j == 4)
                A = X([2 3 4],[2 3 4]) ;f = mu3/ss + mu4/(ss)^(3/2) - 3*(ss)^(1/2);
                g = gradient(f,[ss mu3 mu4]);
            elseif and(i == 3, j == 5)
                A = X([2 3],[2 3]); f = mu3/ss + mu3^2/(ss)^(5/2); 
                g = gradient(f,[ss mu3]);
            elseif and(i == 4, j == 5)
                A = X([2 3 4],[2 3 4]); f = mu4/(ss)^(3/2) - 3*(ss)^(1/2) + mu3^2/(ss)^(5/2); 
                g = gradient(f,[ss mu3 mu4]);     
            end

            v1 = g.'*A*g; 
            x = simplify(subs(v1,ss,sigma^2));
            x = simplify(subs(x,mu3,g1*sigma^3));
            x = simplify(subs(x,mu4,(k+3)*sigma^4));
            S(i,j) = collect((1/2)*(x - S(i,i) - S(j,j)));
            S(j,i) = S(i,j);            
        end
    end

    %mVaR
    x = [1 z (z^2 - 1)/6 (z^3 - 3*z)/24 -(2*z^3 - 5*z)/36];
    V = simplify(x*S*x.'); 
    [N,D] = numden(V);
    var_mVaR = simplify(combine(collect(collect(collect(N,z),sigma),g1)))/D;

    %VaR
    M = S(1:2,1:2); y = [1 C*z];
    var_VaR = simplify(y*M*y.');

    %mETL
    ST  = 1/6* (ITerm(3) - 3*ITerm(1));
    KT  = 1/24*(ITerm(4) - 6*ITerm(2)  + 3*dz);
    S2T = 1/72*(ITerm(6) - 15*ITerm(4) + 45*ITerm(2) - 15*dz);
    x =   1/alpha_sig*[alpha_sig -dz -ST -KT -S2T];

    V = simplify(x*S*x.'); [N,D] = numden(V);
    var_mETL = simplify(combine(collect(collect(collect(collect(N,dz),z),sigma),alpha_sig))/D);

    %ETL
    M = S(1:2,1:2); y = 1/alpha_sig*[alpha_sig -C*dz];
    N = simplify(y*M*y.');
    var_ETL = simplify(collect(collect(N,dz),sigma));
    
    val = [var_mVaR, var_VaR, var_mETL, var_ETL];
end

function val = ITerm(k)

    syms z dz
    isodd = mod(k,2);
    sum = 0;

    if (isodd)    
        for i = 0:(k-1)/2            
            sum = sum + z^(2*i)*prod(2*(0:(k-1)/2) + 1)/prod(2*(0:i) + 1);
        end    

        term = prod(2*(0:(k-1)/2) + 1);
        val = dz*(sum - term);
    else
        for i = 1:(k/2)                    
            sum = sum + z^(2*i)*prod(2*(1:k/2))/prod(2*(1:i));
        end    

        term = prod(2*(1:k/2));
        val = dz*(sum + term);
    end
end