function TestDraw(  )

close all;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Eigenvalues = importdata('Eigenvalues.txt','\t');
% [m n] = size(Eigenvalues);
% startend = sprintf('[1,%d]', m);
% disp(strcat(startend,'rows of data are imported.'));
% valStartEndRow = input(strcat('Input training [start end]rows ', startend,':'));
% startRow = 1; endRow = m;
% if ~isempty(valStartEndRow)
%     startRow = valStartEndRow(1);
%     endRow = valStartEndRow(2);
% end;
% 
% %S = cov(Eigenvalues(startRow:endRow, :));
% 
% [coeff score latent] = princomp(Eigenvalues(startRow:endRow, :));
% 
% valStartEndRow = input(strcat('Input testing [start end]rows ', startend,':'));
% startRow = 1; endRow = m;
% if ~isempty(valStartEndRow)
%     startRow = valStartEndRow(1);
%     endRow = valStartEndRow(2);
% end;
% standardizedData = zscore(Eigenvalues(startRow:endRow,:));
% transferedData = standardizedData*coeff;
% s1 = transferedData(:,1);
% s2 = transferedData(:,2);
% s3 = transferedData(:,3);
% avgS1 = s1; %getAverage(s1);
% avgS2 = s2; %getAverage(s2);
% avgS3 = s3; %getAverage(s3);
% figure;
% subplot(221); plot(avgS1);
% subplot(222); plot(avgS2);
% subplot(223); plot(avgS3);
% 
% avgSsum=sqrt(avgS1.*avgS1+avgS2.*avgS2+avgS3.*avgS3);
% subplot(224); plot(avgSsum);
% 
%     function avg = getAverage(si)
%         avg=zeros(floor(m/30),1);
%         j=1;
%         for i=1:30:m-29
%             avg(j) = mean(si(i:(i+29)));
%             j = j + 1;
%         end
%     end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% load 1-nbin
% load 1-xbin
% figure, bar(xbin_p,nbin_p);
% set(gca,'YLim', [0 max(nbin_p)]);
% xlabel('s','fontsize', 20); ylabel('f','fontsize', 20);
% 
% load 2-nbin
% load 2-xbin
% figure, bar(xbin_p,nbin_p);
% set(gca,'YLim', [0 max(nbin_p)]);
% xlabel('s','fontsize', 20); ylabel('f','fontsize', 20);
% 
% load 3-nbin
% load 3-xbin
% figure, bar(xbin_p,nbin_p);
% set(gca,'YLim', [0 max(nbin_p)]);
% xlabel('s','fontsize', 20); ylabel('f','fontsize', 20);
% 
% load 4-nbin
% load 4-xbin
% figure, bar(xbin_p,nbin_p);
% set(gca,'YLim', [0 max(nbin_p)]);
% xlabel('s','fontsize', 20); ylabel('f','fontsize', 20);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Eigenvalues = importdata('Eigenvalues.txt','\t');
d1 = Eigenvalues(:,1)'; d2 = Eigenvalues(:,2)'; d3 = Eigenvalues(:,3)';

% figure, hist(d1);
% figure, hist(d2);
% figure, hist(d3);

% d1 = zscore(d1);
% d2 = zscore(d2);
% d3 = zscore(d3);

figure, plot3(d1,d2,d3,'.'); axis equal
xlabel('d1');
ylabel('d2');
zlabel('d3');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

end