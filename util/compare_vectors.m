function delta = compare_vectors(A,B)
%compares across two array of vectors (matrices) and return the smallest
%angle of rotation between them i.e. the angle of rotation of A(:,i) onto
%B(:,i)about a perpendicular axis.
A=A./repmat(vecnorm(A),size(A,1),1);
B=B./repmat(vecnorm(B),size(B,1),1);
delta = acos(dot(A,B));
for i=1:size(delta,2)
    if norm(delta(i)) < 10^-3 || isnan(norm(A(:,i))) || isnan(norm(B(:,i)))
        delta(i) = 0;
    end
    if delta(i)>pi %rotate in other direction instead
        delta(i) = 2*pi-delta(i);
    end
end