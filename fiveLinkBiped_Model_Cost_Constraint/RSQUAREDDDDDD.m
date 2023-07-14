Data = load('Fukuchi_Features_DTW_RMS_2023-01-19.mat');

PCA_X = Data.Cost_Components;

% {'\int u^2' '\int du^2' '\int u.dq' '\int (d^2q)^2' '\int (d^3q)^3'}


Variables_Names = {'1. \int (d^2q)^2' '2. \int (d^3q)^2' '3. \int (d^2P)^2' '4. \int (d^3P)^2'...
    '5. \int u^2' '6. \int |u|' '7. \int du^2' '8. \int |du|' '9. \int (d^2u)^2' ...
    '10. \int |d^2u|' '11. \int (d(yGRF))^2' '12. \int (d(xGRF))^2' ...
    '13. \int (u.dq)^+' '14. \int |u.dq|' '15. max(1/2mV^2)' '16. range(I\omega)'};
Variables_Names = {'1. Angular Acceleration' '2. Angular Jerk' '3. Linear Acceleration' '4. Linear Jerk'...
    '5. Torque Squared' '6. Absolute Torque' '7. Torque Rate Squared' '8. Absolute Torque Rate' '9. Torque 2nd Rate Squared' ...
    '10. Absolute Torque 2nd Rate' '11. yGRF Rate' '12. xGRF Rate' ...
    '13. Positive Work Rate' '14. Absolute Work Rate' '15. Kinetic Energy' '16. Pk-to-pk Angular Momentum'};

len_cost = size(PCA_X,2);
idd_subj = 0;
Rsq = ones(16,16)*nan;
% for id_subj = 1:No_Subj
for i = 1:len_cost
    for j=i:len_cost
        Rsq(i,j) = corr2(PCA_X(:,i),PCA_X(:,j))^2;
    end
end
CorrelationPlot(PCA_X, Variables_Names);
set(gca, 'XTIck', [1:len_cost])
set(gca, 'YTIck', [1:len_cost])
set(gca, 'XAxisLocation', 'top')
set(gca, 'YAxisLocation', 'right')

textStrings = num2str(abs(Rsq), '%0.2f');       % Create strings from the matrix values
% textStrings = [num2str(abs(Rsq_Mean), '%0.2f') '+/-' num2str(abs(Rsq_SD), '%0.2f')];       % Create strings from the matrix values
% textStrings = num2str(Rsq_Mean, '%0.2f');       % Create strings from the matrix values
% textStrings_SD = num2str(Rsq_SD, '%0.2f');       % Create strings from the matrix values
% textStrings = strtrim(cellstr(textStrings));  % Remove any space padding
text_modif = strings(len_cost,len_cost);
text_modif_Mean = strings(len_cost,len_cost);
text_modif_SD = strings(len_cost,len_cost);
for k = 1:len_cost
    for i = 1:k
        if i~=k & Rsq(i,k)>.3
            if Rsq(i,k)<0
                txt_sgn = '-';
            else
                txt_sgn = '';
            end
            text_modif (k,i)= [txt_sgn, textStrings(i,4*(k-1)+1:4*k)];
            %     text_modif (i,k)= [txt_sgn, textStrings(4*(k-1)+1:4*k),i];
            % text_modif (i,k)= [txt_sgn, textStrings(i,4*(k-1)+1:4*k), char(10), char(177), textStrings_SD(i,4*(k-1)+1:4*k)];
            % text_modif_Mean (i,k)= [txt_sgn, textStrings(i,4*(k-1)+1:4*k)];
            % text_modif_SD (i,k)= [char(177), textStrings_SD(i,4*(k-1)+1:4*k)];
        end
    end
end
[xs, ys] = meshgrid(1:len_cost);  % Create x and y coordinates for the strings
text(xs(:), ys(:), text_modif(:), ...  % Plot the strings
'HorizontalAlignment', 'center');
set( gca, 'XTickLabel', Variables_Names','FontSize',12 )
xtickangle(70)

X_scaled = (PCA_X - min(PCA_X)) ./ std(PCA_X) + 1;
figure(322)
% subplot(1,2,2)
% corrplot(X_scaled, 'VarNames', Variables_Names)



len_cost = size(PCA_X,2);
idd_subj = 0;
Rsq = ones(16,16)*nan;
% for id_subj = 1:No_Subj
for i = 1:len_cost
    for j=i:len_cost
        Rsq(i,j) = corr2(PCA_X(:,i),PCA_X(:,j))^2;
    end
end


% Rsq_Mean = squeeze(mean(Rsquared,1));
% Rsq_SD = squeeze(std(Rsquared,[],1));
% Rsq = Rsq_Mean;


figure(13)
% % imagesc(abs(Rsq))
imagesc(Rsq)
% impixelregion(imagesc(Rsq))
set(gca, 'XTIck', [1:len_cost])
set(gca, 'YTIck', [1:len_cost])
set(gca, 'XAxisLocation', 'top')
set(gca, 'YAxisLocation', 'right')
colorbar('Location','eastoutside')
colorbar('Location','southoutside')
colormap Bone
textStrings = num2str(abs(Rsq), '%0.2f');       % Create strings from the matrix values
% textStrings = [num2str(abs(Rsq_Mean), '%0.2f') '+/-' num2str(abs(Rsq_SD), '%0.2f')];       % Create strings from the matrix values
% textStrings = num2str(Rsq_Mean, '%0.2f');       % Create strings from the matrix values
% textStrings_SD = num2str(Rsq_SD, '%0.2f');       % Create strings from the matrix values
% textStrings = strtrim(cellstr(textStrings));  % Remove any space padding
text_modif = strings(len_cost,len_cost);
text_modif_Mean = strings(len_cost,len_cost);
text_modif_SD = strings(len_cost,len_cost);
for k = 1:len_cost
    for i = 1:len_cost
        if i~=k
            if Rsq(i,k)<0
                txt_sgn = '-';
            else
                txt_sgn = '';
            end
            text_modif (i,k)= [txt_sgn, textStrings(i,4*(k-1)+1:4*k)];
            %     text_modif (i,k)= [txt_sgn, textStrings(4*(k-1)+1:4*k),i];
            % text_modif (i,k)= [txt_sgn, textStrings(i,4*(k-1)+1:4*k), char(10), char(177), textStrings_SD(i,4*(k-1)+1:4*k)];
            % text_modif_Mean (i,k)= [txt_sgn, textStrings(i,4*(k-1)+1:4*k)];
            % text_modif_SD (i,k)= [char(177), textStrings_SD(i,4*(k-1)+1:4*k)];
        end
    end
end


[xs, ys] = meshgrid(1:len_cost);  % Create x and y coordinates for the strings
text(xs(:), ys(:), text_modif(:), ...  % Plot the strings
'HorizontalAlignment', 'center');
set( gca, 'XTickLabel', Variables_Names','FontSize',14 )
xtickangle(70)
set( gca, 'YTickLabel', Variables_Names','FontSize',14 )
% set( gca, 'YTickLabel', {'angl accel sqrd' 'angl jerk sqrd' 'cartes accel sqrd' 'cartes jerk sqrd' 'torqs sqrd' 'torqs abs' 'torq rate sqrd' 'torq rate abs' 'torq rate change sqrd' 'torq rate change abs' 'yGRF rate sqrd' 'xGRF rate sqrd' 'positive work' 'absolute work' 'kinetic energy' 'angular momentum'}','FontSize',12 )
















len_cost = size(PCA_X,2);
idd_subj = 0;
% for id_subj = 1:No_Subj
for id_subj = 1:31
    if id_subj~=5 && id_subj~=17 && id_subj~=41
        idd_subj = idd_subj+1;
        Rsquared(idd_subj,:,:) = reshape(R_Squared_Bootstrap(Data.Cost_Components_Subjects{idd_subj,1}),len_cost,len_cost)';
    end
end


Rsq_Mean = squeeze(mean(Rsquared,1));
Rsq_SD = squeeze(std(Rsquared,[],1));
% Rsq = Rsq_Mean;




% Variables_Names = {'1.angl accel sqrd' '2.angl jerk sqrd' '3.cartes accel sqrd' '4.cartes jerk sqrd' '5.torqs sqrd' '6.torqs abs' '7.torq rate sqrd' '8.torq rate abs' '9.torq rate change sqrd' '10.torq rate change abs' '11.yGRF rate sqrd' '12.xGRF rate sqrd' '13.positive work' '14.absolute work' '15.pk kinetic energy' '16.pk2pk angular momentum'};


figure(12)
set(gcf,'Color','white','renderer','Painters','Units','Centimeters','Position',[1 1 36 50]);
% % imagesc(abs(Rsq))
imagesc(Rsq)
% % imagesc(Rsq_Mean-Rsq_SD)
% impixelregion(imagesc(Rsq))
set(gca, 'XTIck', [1:len_cost])
set(gca, 'YTIck', [1:len_cost])
set(gca, 'XAxisLocation', 'top')
set(gca, 'YAxisLocation', 'right')
colorbar('Location','eastoutside')
colorbar('Location','southoutside')
colormap Bone
textStrings = num2str(abs(Rsq), '%0.2f');       % Create strings from the matrix values
% textStrings = [num2str(abs(Rsq_Mean), '%0.2f') '+/-' num2str(abs(Rsq_SD), '%0.2f')];       % Create strings from the matrix values
% textStrings = num2str(Rsq_Mean, '%0.2f');       % Create strings from the matrix values
textStrings_SD = num2str(Rsq_SD, '%0.2f');       % Create strings from the matrix values
% textStrings = strtrim(cellstr(textStrings));  % Remove any space padding
text_modif = strings(len_cost,len_cost);
text_modif_Mean = strings(len_cost,len_cost);
text_modif_SD = strings(len_cost,len_cost);
for i = 1:len_cost
    for k = i:len_cost
        if i~=k
            
            if Rsq(i,k)<0
                txt_sgn = '-';
            else
                txt_sgn = '';
            end
            %     text_modif (i,k)= [txt_sgn, textStrings(i,4*(k-1)+1:4*k)];
            text_modif (i,k)= [txt_sgn, textStrings(i,4*(k-1)+1:4*k), char(10), char(177), textStrings_SD(i,4*(k-1)+1:4*k)];
            text_modif_Mean (i,k)= [txt_sgn, textStrings(i,4*(k-1)+1:4*k)];
            text_modif_SD (i,k)= [char(177), textStrings_SD(i,4*(k-1)+1:4*k)];
        end
    end
end
[xs, ys] = meshgrid(1:len_cost);  % Create x and y coordinates for the strings
text(xs(:), ys(:), text_modif(:), ...  % Plot the strings
'HorizontalAlignment', 'center');
set( gca, 'XTickLabel', Variables_Names','FontSize',12 )
xtickangle(70)
set( gca, 'YTickLabel', Variables_Names','FontSize',12 )
title('Correlation Matrix for Cost Components')
% set( gca, 'YTickLabel', {'angl accel sqrd' 'angl jerk sqrd' 'cartes accel sqrd' 'cartes jerk sqrd' 'torqs sqrd' 'torqs abs' 'torq rate sqrd' 'torq rate abs' 'torq rate change sqrd' 'torq rate change abs' 'yGRF rate sqrd' 'xGRF rate sqrd' 'positive work' 'absolute work' 'kinetic energy' 'angular momentum'}','FontSize',12 )

print('Name_of_file','-dpdf')

print(['filename.pdf'])
















function Rsq = R_Squared_Bootstrap (Cost_Components_Eval_Fakuchi)

len_cost = size(Cost_Components_Eval_Fakuchi,2);
for i = 1:len_cost
    for j = i:len_cost
        

        Rsq(i,j) = corr2(Cost_Components_Eval_Fakuchi(:,i),Cost_Components_Eval_Fakuchi(:,j))^2;
        
        
    end
end
Rsq = reshape(Rsq',1,len_cost*len_cost);
end
