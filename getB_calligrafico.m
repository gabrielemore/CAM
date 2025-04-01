function [B_call,B_call_last] = getB_calligrafico(A,B,N)
   zeri = zeros(size(B));
   B_call = [zeros(size(B)); B];
    [n,~] = size(A);

   for i = 1:N
        ultima_riga = [];
        for k = i:-1:1
            ultima_riga = [ultima_riga, A^k*B];
        end
        B_call = [B_call;ultima_riga];

        zeri = [zeros(size(B)); zeri];
        B_call = [B_call, [zeri;B]];

   end

   B_call_last = B_call(end-3:end,:);



