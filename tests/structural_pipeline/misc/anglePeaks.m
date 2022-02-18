function angDeg = anglePeaks(peaksTest, peaksRef, indx)

if nargin < 3
    indx = true(size(peaksTest,1),1);
end

for i = 1:size(peaksTest, 3)
    D1 = peaksTest(indx, :, 1);
    D2 = peaksRef(indx, :, :);
    ang = acos(sum(D1 .* D2, 2));
    
    ang = 0.5*pi - abs(0.5*pi - ang);
    angDeg = rad2deg(ang);
    
    angDeg = min(angDeg, [], 3);
    
end

end