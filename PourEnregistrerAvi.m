%% Faire un film .avi à partir du tenseur de données enregistré dans la
%% workspace
%%%

 figure(70);
 fig=gcf;
 %aviobj = avifile('MEEFWR085CinétiqueCorr-t5.avi');
 %aviobj = avifile('MEEFWR085ReplicatCorrCoupeeEch2Mes3.avi');
 aviobj = avifile('MEEFWR085InterferantCorr2.avi');
 %aviobj = avifile('MEEFWR085Cinétique-t2.avi');
 %aviobj = avifile('MEEFWR085CinétiqueCoupee-t2.avi');
 aviobj.fps=3;
 aviobj.quality=100;
 for i=1:length(FileName)
     %contourf(EmWaveLength1,ExcCoupee,10*log10(eem_corCoupe(:,:,i)),20) 
     %contourf(EmWaveLength1,ExcCoupee,10*log10(MeefCoupee(:,:,i)),20)
     %contourf(EmWaveLength1,ExcWaveLength1,10*log10(MeefPlomb(:,:,i)),20) 
     contourf(EmWaveLength1,ExcWaveLength1,10*log10(eem_cor(:,:,i)),20) 
     
%       h = gca;
% %      
%       hold on
% % 
% % % Le polynome modélisant les pics
%      pic = zeros(4, length(EmWaveLength1));
% % 
% % % Pour chaque pic
%      for ii = 1 : 4
%          %     Valeurs du polynome aux points correspondant à la longueur d'onde d'excitation
%          pic(ii, :) = polyval(coeff(ii, :), EmWaveLength1);
%          %     Affichage par dessus l'image
%          plot(EmWaveLength1, pic(ii, :), 'o', 'parent', h);
%      end
   
     set(gca,'YDir','normal')
     title('MEEF en présence de métal lourd (Plomb), en échelle logarithmique (dB)');
     xlabel('Longueur d''onde \lambda d''émission (nm)');
     ylabel('Longueur d''onde \lambda d''excitation (nm)');  
     caxis([-5 20])
     colorbar 
     %name=strcat({FileName{i}(1:end-12)},{'Coupee+posRaies'})
     name=strcat({FileName{i}(1:end-12)},{'corrCoupe'})
     %name=strcat({FileName{i}(1:end-12)},{'corr'})
     %name=strcat({FileName{i}(1:end-12)},{'+posRaies'})
     %h=gca;
     %figs = findobj(0, 'type', 'figure') 
     %fig = get(groot,'CurrentFigure');
     
     nameFin=strcat({'SauveVero/MEEF WR085 + adjuvants/'},name)
     print(fig,'-djpeg100',[nameFin{1} '.jpg'],'-r96') ;
       
     F(i) = getframe(gcf);
     aviobj = addframe(aviobj,F(i));
 end
 %h=gcf;

 aviobj = close(aviobj);
