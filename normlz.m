function [M,Q] = normlz(ref,i,M,Q,R)
x=size(M{ref(1)}) %changed from M[i] to m2 to make it work for dimensions across pt and mo
for k = 1:length(M)
    if length(M{k})< x(1)
        x(1) = length(M{k})
    end
    if width(M{k})< x(2)
        x(2) = width(M{k})
    end    
end %this now makes in both dimensions the smallest length the one that is divided for continuity


normel = R{ref(1)}(1:x(1),ref(2)+1); %Constant normalisation reference that doesn't change throughout the normalisation (as though M is changed, R is not)

for j = 1:x(2)-1     %have changed from numel(Q{ref(1)}) 
    M{i}(1:x(1),j+1) = (M{i}(1:x(1),j+1))./normel; %have changed : to 1:x(1) in a number of places to get dimensions to agree but something slicker needs putting in that finds the smallest dimension and bases it off that.
    %M{i}(~isfinite(M{i})) = 1
    M{i}(isinf(M{i})|isnan(M{i})) = 1.0
end
end


