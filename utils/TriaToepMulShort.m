function M = TriaToepMulShort(c,P,rev)
% computes the matrix multiplication of a lower triangular toeplitz matrix
% with first column c with its transposed by the computationally efficient
% recursive formula

    if rev == false
        M = zeros(P,P);
        for i = [1:length(c)]
            for j = [1:length(c)]
                if i == 1
                    M(i,i+j-1) = c(i) * c(i+j-1);
                    M(i+j-1,i) = M(i,i+j-1);
                else
                    M(i,min(i+j-1,length(c))) = M(i-1,min(i+j-1-1,length(c)-1)) + c(i)*c(min(i+j-1,length(c)));
                    M(min(i+j-1,length(c)),i) = M(i,min(i+j-1,length(c)));
                end
            end
        end
        for i = [length(c)+1:P]
            for j = [1:length(c)]
                M(min(i-j+1,P),i) = M(min(i-j,P-1),i-1);
                M(i,min(i-j+1,P)) = M(min(i-j+1,P),i);
            end
        end
    else
        cflip = flip(c(2:end));
        M = zeros(P,P);
        Mfull = zeros(length(c)-1,length(c)-1);
        for i = [1:length(c)-1]
            for j = [i:length(c)-1]
                if i == 1
                    Mfull(i,j) = cflip(i) * cflip(j);
                    Mfull(j,i) = Mfull(i,j);
                else
                    Mfull(i,j) = Mfull(i-1,j-1) + cflip(i)*cflip(j);
                    Mfull(j,i) = Mfull(i,j);
                end
            end
        end
        M(end-length(c)+2:end,end-length(c)+2:end) = Mfull;
    end
end