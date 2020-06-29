function coeff = treat(EmWL, ExcWL, images)
close all
A = size(images)
average = mean(images,3);
imagesc(EmWL, ExcWL, 10*log10(average))
set(gca,'YDir','normal')
figure
T = adaptthresh(average, 0.2);

bnw = imbinarize(average,T)
imagesc(EmWL, ExcWL, bnw)

size(average)
end

