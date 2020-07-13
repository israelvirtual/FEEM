
%Final tasks
% - encapsulate regressions
% - solve problem with svm ok
% - write error function ok
% - tests on the correction matrix to see corrections ok
    % result, good for all except coupee
% - make ihm more user friendly ok


clear all
close all


chem0=[pwd '\DonnéesFluorescence'];
chem1=[pwd '\DonnéesFluorescence\MEEF Pb Seul'];
chem3=[pwd '\NosFonctions'];


addpath(chem0,chem1,chem3);

%%% Si on commence veut charger le tampon
%[EmWaveLength1, ExcWaveLength1, MeefTampon] = getImage3DFromFichier('M_3D_Layglon_S_Tampon_C_N_T_Rbis.txt');
% c'est ok même si on lit une matrice rectangulaire !!!
[EmWaveLength1, ExcWaveLength1, MeefTampon] = getImage3DFromFichier('M_3D_Layglon_S_Tampon_C_N_T_R.txt');
figure(1)
imagesc(EmWaveLength1, ExcWaveLength1, log10(MeefTampon))
set(gca,'YDir','normal')
title('MEEF du tampon, en échelle logarithmique (dB)');
xlabel('Longueur d''onde \lambda d''émission (nm)');
ylabel('Longueur d''onde \lambda d''excitation (nm)');
colorbar



%%% Si on veut charger une MEEF quelconque
[EmWaveLength1, ExcWaveLength1, MeefPlomb10] = getImage3DFromFichier('M_3D_Layglon_S_10eq_C_N_T_R.txt');
figure(2)
imagesc(EmWaveLength1, ExcWaveLength1, log10(MeefPlomb10))
set(gca,'YDir','normal')
%hold off
title('MEEF en présence de métal lourd (Plomb), en échelle logarithmique (dB)');
xlabel('Longueur d''onde \lambda d''émission (nm)');
ylabel('Longueur d''onde \lambda d''excitation (nm)');
colorbar


%%% Si on veut charger plusieurs MEEF de même taille pour construire le cube à décomposer...
[FileName,PathName] = uigetfile('*.txt','Select the .txt file','MultiSelect','on');

if(~iscell(FileName) && (isequal(FileName,0) || isequal(PathName,0)))
    msgbox('Erreur : aucun fichier ouvert');
    return
end

chaine = '';
if iscell(FileName) == 0
    chaine = [PathName, FileName];
else
    L = min(length(EmWaveLength1),length(ExcWaveLength1));
    h = waitbar(0,'Please wait...');
    for ii = 1 : length(FileName)
        images = getImage3DFromFichier([PathName FileName{ii}]);
       [EmWaveLength1, ExcWaveLength1, MeefPlomb(:,:,ii)]=getImage3DFromFichier([PathName FileName{ii}]);
       waitbar(ii/length(FileName),h)
    end
    delete(h)     
end
% puis afficher la dernière de celles que l'on a chargé
figure(3)
imagesc(EmWaveLength1, ExcWaveLength1, 10*log10(MeefPlomb(:,:,end)))
set(gca,'YDir','normal')
title('MEEF en présence de métal lourd (Plomb), en échelle logarithmique (dB)');
xlabel('Longueur d''onde \lambda d''émission (nm)');
ylabel('Longueur d''onde \lambda d''excitation (nm)');
colorbar


%%%  Affichage en courbes de niveau avec contourf
figure(4)
contourf(EmWaveLength1, ExcWaveLength1, 10*log10(MeefPlomb(:,:,end)),20)
set(gca,'YDir','normal')
title('MEEF en présence de métal lourd (Plomb), en échelle logarithmique (dB)');
xlabel('Longueur d''onde \lambda d''émission (nm)');
ylabel('Longueur d''onde \lambda d''excitation (nm)');
colorbar

%%%  Si on souhaite couper un bout de la MEEF car les faibles excitations
%%%  sont e fait du bruit -> on observe la même chose dans la solution
%%%  tampon alors qu'il n'y a rien dedans...

MeefCoupee=MeefPlomb(27:end,:,end);
ExcCoupee=ExcWaveLength1(27:end);
figure(5)
contourf(EmWaveLength1, ExcCoupee, 10*log10(MeefCoupee),20)
set(gca,'YDir','normal')
title('MEEF en présence de métal lourd (Plomb), en échelle logarithmique (dB)');
xlabel('Longueur d''onde \lambda d''émission (nm)');
ylabel('Longueur d''onde \lambda d''excitation (nm)');
colorbar

%%%%%
%% Pour éliminer les raies de diffusion Raman et Rayleigh
close all
%Average of images

MeefPlombC=MeefPlomb(27:end,:,:);                       %Cut image
average = mean(MeefPlombC,3)/max(MeefPlombC,[],'all');
I = imagesc(EmWaveLength1, ExcWaveLength1, 10*log10(average));
set(gca,'YDir','normal')

%Thresholding : The adaptive feature finds the correct thresehold

bnw = imbinarize(average,'adaptive');

% Finding the angle of the structuring element (here we use the fact that
% the rays are roughly linear for a later opening

angle = angle1(bnw, ExcCoupee)
 
% Preprocessing the image to isolate the rays

%opening


% Average filtering
kernel = 1/4*[1 1 ;1 1 ]
bnw = imfilter(bnw, kernel);

se = strel('line', 7, angle);               %Here, just angle has to be precise, the size 7 can be different

bnw = imopen(bnw,se);
figure
imagesc(EmWaveLength1, ExcCoupee, bnw)
set(gca,'YDir','normal')
title('Final processing')
pause
% We count the number of rays and initiallize coeff and correct as empty

CC = bwconncomp(bnw);
coeff = [];
correct = [];

numObjects = CC.NumObjects
%Idea: we see that the hand correction takes into consideration the
%intensity of the the fluorescence. Maybe ask later about that.

% Here we isolate each ray in its own image, we map them the wavelengths
% and interpolate
for i = 1:CC.NumObjects
    bn(:,:,i) = bwpropfilt(logical(bnw),'perimeter',1);
    bnw = bnw - bn(:,:,i);
     figure
     imagesc(EmWaveLength1, ExcWaveLength1, bn(:,:,i))
         set(gca,'YDir','normal')
    [data(:,:,i) corr] = findcenterpixel(bn(:,:,i));
    correct = [correct; corr corr];
    figure
    imagesc(EmWaveLength1, ExcCoupee, data(:,:,i))
    
    set(gca,'YDir','normal')
%     figure
%     plot(donnees(:,:,i));
x = []; y = [];
for k = 1:size(data(:,:,i),1)
    for L = 1:size(data(:,:,i),2)
        if data(k,L,i) ==1
            y = [y ExcCoupee(k)];
            x = [x EmWaveLength1(L)];
        end
    end
end
% y = y';
% p = polyfit(x,y,2);
% coeff = [coeff;(polyfit(x,y,2))];
% y1 = polyval(p,EmWaveLength1);
% plot(x,y);
% error = rmse(y',y1(1:size(y)))
% pause
tic
pModel = fitrsvm(x',y','KernelFunction', 'polynomial','PolynomialOrder',2,'KernelScale','auto','Standardize',true,'OptimizeHyperparameters','auto',...
    'HyperparameterOptimizationOptions',struct('AcquisitionFunctionName',...
    'expected-improvement-plus'))
y1 = predict(pModel,EmWaveLength1');
y2 = predict(pModel, x');
error = rmse(y', y2(1:size(y)))
coeff = [coeff ; polyfit(EmWaveLength1,y1,2)]
plot (EmWaveLength1,y1)
toc
end

close all

%% Second method: GUI
%Average of images
close all
MeefPlombC=MeefPlomb(27:end,:,:);                       %Cut image
average = mean(MeefPlombC,3)/max(MeefPlombC,[],'all');
I = imagesc(EmWaveLength1, ExcWaveLength1, 10*log10(average));
set(gca,'YDir','normal')

numObjects = inputdlg({'Type number of rays'}, 'Number of rays', [1 10])
numObjects = str2double(numObjects);
%regressionLearner

%% Regression 1 SVM
correct = []; coeff=[];

for i = 1:numObjects
I = imagesc(EmWaveLength1, ExcCoupee, 10*log10(average));
set(gca,'YDir','normal')
title(['Choose points on the ', num2str(i), 'th biggest rays for regression. Press enter to validate'])
[X Y] = getpts;
pModel = fitrsvm(X,Y,'KernelFunction', 'polynomial','PolynomialOrder',2,'KernelScale','auto','Standardize',true,'OptimizeHyperparameters','auto',...
    'HyperparameterOptimizationOptions',struct('AcquisitionFunctionName',...
    'expected-improvement-plus'))
pModel.ModelParameters
y1 = predict(pModel,EmWaveLength1');
coeff = [coeff ; polyfit(EmWaveLength1,y1,2)];
y2 = predict(pModel, X);
error = rmse(Y, y2)
plot(EmWaveLength1,y1)
pause

I = imagesc(EmWaveLength1, ExcWaveLength1, 10*log10(average));
set(gca,'YDir','normal')
title(['Choose 2 points on the ', num2str(i), 'th ray for correction. Press enter to validate'])
[X Y] = getpts;
correct = [correct; [abs(Y(1)-Y(2)) abs(Y(1)-Y(2))]];
end

close all
  
%% Regression 2 Polyfit

coeff = []; correct = [];
for i = 1:numObjects
I = imagesc(EmWaveLength1, ExcCoupee, 10*log10(average));
set(gca,'YDir','normal')
title(['Choose points on the ', num2str(i), 'th biggest rays for regression. Press enter to validate'])
[X Y] = getpts;

p = polyfit(X,Y,2);
coeff = [coeff;p];
y1 = polyval(p,EmWaveLength1');
y2 = polyval(p, X);
error = rmse(Y, y2)
figure
plot(EmWaveLength1,y1)
pause

I = imagesc(EmWaveLength1, ExcCoupee, 10*log10(average));
set(gca,'YDir','normal')
title(['Choose 2 points on the ', num2str(i), 'th biggest rays for correction. Press enter to validate'])
[X Y] = getpts;
correct = [correct; [abs(Y(1)-Y(2)) abs(Y(1)-Y(2))]];
end

close all

%% Regression 3 fitnlm

coeff = []; correct = [];
for i = 1:numObjects
I = imagesc(EmWaveLength1, ExcCoupee, 10*log10(average));
set(gca,'YDir','normal')
title(['Choose points on the ', num2str(i), 'th biggest rays for regression. Press enter to validate'])
[X Y] = getpts;
TB = table(X,Y);
modelfun = @(b,x) b(1)+b(2)*x+b(3)*x.^2;
beta0 = [0 1 1];
pModel = fitnlm(TB,modelfun, beta0);
coeff = [coeff; (flip(pModel.Coefficients.Estimate))']
y1 = predict(pModel,EmWaveLength1');
y2 = predict(pModel, X);
error = rmse(Y, y2)
plot(EmWaveLength1,y1)
pause

I = imagesc(EmWaveLength1, ExcCoupee, 10*log10(average));
set(gca,'YDir','normal')
title(['Choose 2 points on the ', num2str(i), 'th biggest rays for correction. Press enter to validate'])
[X Y] = getpts;
correct = [correct; [abs(Y(2)-Y(1)) abs(Y(2)-Y(1))]];
end

close all
  
%% Extra Regression 1 Polyfit g

coeff = []; correct = [];
for i = 1:numObjects 
I = imagesc(EmWaveLength1, ExcCoupee, 10*log10(average));
set(gca,'YDir','normal')
title(['Choose points on the ', num2str(i), 'th biggest rays for regression. Press enter to validate'])
[X Y] = getpts;

p = polyfit(X,Y,1);
coeff = [coeff;p];
y1 = polyval(p,EmWaveLength1');
y2 = polyval(p, X);
error = rmse(Y, y2)
figure
plot(EmWaveLength1,y1)
pause

I = imagesc(EmWaveLength1, ExcCoupee, 10*log10(average));
set(gca,'YDir','normal')
title(['Choose 2 points on the ', num2str(i), 'th biggest rays for correction. Press enter to validate'])
[X Y] = getpts;
correct = [correct; [abs(Y(1)-Y(2)) abs(Y(1)-Y(2))]];
end

close all

%%
% coeff = [ 0      1    4; 
%           0.000    0.81 35;
%           0      2.025      0
%         
%          ];
% 
%    
%%% polynome en ax^2+bx+c
%%% première ligne : raie la plus forte au milieu de la MEEF
%%% 2ème ligne raie située en dessous
%%% 3ème ligne raie bord inférieur droit de la MEEF
%%% 4ème ligne  : raie bord supérieur gauche de la MEEF
      
% Largeur de chaque raie.
% 4 raies ici
% Première colonne, la largeur 'en-dessous' du pic, deuxième colonne,
% la largeur en dessus

% correct = [ 35 35;
%             15 15;
%             25 25];


 figure(22)
 imagesc( EmWaveLength1,ExcWaveLength1, log10(MeefPlomb(:,:,end)))
 set(gca, 'YDir', 'normal')
 h = gca;
 hold on

 
%%% pour superposer les raies sur l'image des MEEFs



% Le polynome modélisant les pics
pic = zeros(4, length(EmWaveLength1));

% Pour chaque pic
for ii = 1 : numObjects
    %     Valeurs du polynome aux points correspondant à la longueur d'onde
    %     d'excitation
    pic(ii, :) = polyval(coeff(ii, :), EmWaveLength1);
    %
    %     Affichage par dessus l'image
    plot(EmWaveLength1, pic(ii, :), 'o', 'parent', h);
end

title('MEEF en présence de métal lourd (Plomb), en échelle logarithmique (dB) + position Raman et Rayleigh');
xlabel('Longueur d''onde \lambda d''émission (nm)');
ylabel('Longueur d''onde \lambda d''excitation (nm)');
colorbar

%%% Meme chose mais en considérant des MEEFS que l'on a tronqué....
%%%

figure(25)
 imagesc(EmWaveLength1,ExcCoupee,10*log10(MeefCoupee))
 set(gca, 'YDir', 'normal')
 h = gca;
 hold on

% Le polynome modélisant les pics
pic = zeros(4, length(EmWaveLength1));

% Pour chaque pic
for ii = 1 : numObjects
    %     Valeurs du polynome aux points correspondant à la longueur d'onde d'excitation
    pic(ii, :) = polyval(coeff(ii, :), EmWaveLength1);
    %     Affichage par dessus l'image 
    plot(EmWaveLength1, pic(ii, :), 'o', 'parent', h);
end

title('MEEF en présence de métal lourd (Plomb), en échelle logarithmique (dB) + position Raman et Rayleigh');
xlabel('Longueur d''onde \lambda d''émission (nm)');
ylabel('Longueur d''onde \lambda d''excitation (nm)');
colorbar

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%   Corrections de la MEEF tronquée par Zepp
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

[eem_cor,correct,eem_filter] = traitementZepp(MeefCoupee(:,:,end), correct, coeff,  EmWaveLength1,ExcCoupee);

figure(50)
contourf(EmWaveLength1,ExcCoupee,10*log10(eem_cor),20)
set(gca,'YDir','normal')
title('MEEFCoupee en présence de métal lourd (Plomb), en échelle logarithmique (dB) Corrigee (IMPOTANTE)');
xlabel('Longueur d''onde \lambda d''émission (nm)');
ylabel('Longueur d''onde \lambda d''excitation (nm)');
colorbar

figure(51)
imagesc(EmWaveLength1,ExcCoupee,eem_filter(2:end,2:end))
set(gca,'YDir','normal')
title('Masque utilisé');
xlabel('Longueur d''onde \lambda d''émission (nm)');
ylabel('Longueur d''onde \lambda d''excitation (nm)');
colorbar


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%   Corrections de la MEEF on tronquée par Zepp
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
[eem_cor,correct,eem_filter] = traitementZepp(MeefPlomb(:,:,end), correct, coeff, EmWaveLength1, ExcWaveLength1);
figure(60)
contourf(EmWaveLength1,ExcWaveLength1,10*log10(eem_cor),20)
set(gca,'YDir','normal')
title('MEEF Normal en présence de métal lourd (Plomb), en échelle logarithmique (dB) Corrigee');
xlabel('Longueur d''onde \lambda d''émission (nm)');
ylabel('Longueur d''onde \lambda d''excitation (nm)');
colorbar

figure(61)
imagesc(EmWaveLength1,ExcWaveLength1,eem_filter(2:end,2:end))
set(gca,'YDir','normal')
title('Masque utilisé');
xlabel('Longueur d''onde \lambda d''émission (nm)');
ylabel('Longueur d''onde \lambda d''excitation (nm)');
colorbar




%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%   Corrections de la MEEF du tampon par Zepp
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

[eem_cor,correct,eem_filter] = traitementZepp(MeefTampon(:,:,end), correct, coeff,  EmWaveLength1,ExcWaveLength1);

figure(70)
contourf(EmWaveLength1,ExcWaveLength1,10*log10(eem_cor),20)
set(gca,'YDir','normal')
title('MEEF du Tampon, en échelle logarithmique (dB)');
xlabel('Longueur d''onde \lambda d''émission (nm)');
ylabel('Longueur d''onde \lambda d''excitation (nm)');
colorbar

figure(71)
imagesc(EmWaveLength1,ExcWaveLength1,eem_filter(2:end,2:end))
set(gca,'YDir','normal')
title('Masque utilisé');
xlabel('Longueur d''onde \lambda d''émission (nm)');
ylabel('Longueur d''onde \lambda d''excitation (nm)');
colorbar

function [data, correction] = findcenterpixel(BinImage)
data = zeros(size(BinImage));
for i = 1:size(BinImage,1)
    taille = 0; min = 0;
    for j = 1:size(BinImage,2)
        if BinImage(i,j) ~=0
            if BinImage(i,j-1) == 0
                min = j;
            end
            taille = taille+1;
        end
    end
    correct(i) = taille*5;
    %zero any non central pixel
    for j = 1:size(BinImage,2)
        if BinImage(i,j) ~=0 && j ==min+round(taille/2)
            data(i,j-1) = 1;
        end
    end 
end
correction = max(correct);

end

function angle = angle1 (bnw,ExcCoupee)
    Y1 = 0; Y2=0;
    sep = bwpropfilt(logical(bnw),'perimeter',1);
    [sep corr] = findcenterpixel(sep);
    for i = 1:size(ExcCoupee,2)
        if sep(5,i) == 1
            Y1 = i;
        end
        if sep (10,i) == 1
            Y2= i;
        end
    end
    angle = atand((Y2-Y1)/5) + 90;
end

function rootmean = rmse(set1,set2)
    rootmean = sqrt(sum((set1 - set2).^2)/length(set1));
end
