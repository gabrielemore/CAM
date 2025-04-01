function [A_call,A_call_last] = getA_calligrafico(A,N)
    A_call = [A^0];
    for i = 1:N
        A_call = [A_call;A^i];
        if(i==N)
            A_call_last = A^i;
        end
    end