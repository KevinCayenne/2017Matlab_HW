function fv=ex2(a,n)

    record = ones(n+1,7)*-9999;
    record(1,1) = 0;    
    record(1,2) = a;
    record(1,5) = a;

    for i=1:1:n
        a1 = record(i,2);
        a2 = record(i,5);
        
        fn = 5*a1.^5-3*a1+6;
        fnd = 25*a1.^4-3;

        fndopt = 25*a2.^4-3;
        fnd2opt = 100*a2.^3;

        nx=record(i,2)-fn/fnd;
        nx2=record(i,5)-fndopt/fnd2opt;

        record(i+1,1)=i;
        record(i+1,2)=nx;
        record(i+1,3)=fn;
        record(i+1,4)=fnd;
        record(i+1,5)=nx2;
        record(i+1,6)=fndopt;
        record(i+1,7)=fnd2opt;
    end

save record;
fv=record;

end