function [ distance xbin nbin confidence] = ...
    st_covmat_analyse(videofile, duration, filtersize, handles, conf, d, x, n)

    % Video information
    video = mmreader(videofile);
    numberOfFrames = video.NumberOfFrames;
    Height = video.Height;
    Width = video.Width;
    length = video.duration;

    distance = []; xbin =[]; nbin = []; confidence = []; occurProb = [0];
    xx = x; nn = n; dd = d; cc = conf; detect_distance = [];
    posImagShow = get(handles.imag, 'Position');   
    
    % Add 9/29/2011 Mahalanobis distance
    mmu = []; mcov = []; malpha = 0.1; minX=0; maxX=Inf;
  
    fid_results = fopen('Eigenvalues.txt','wt+');
    if fid_results == 0
        disp('Create log file failed!\n');
        return;
    end
    fid_agnles = fopen('Anglevalues.txt','wt+');
    if fid_agnles == 0
        disp('Create log file failed!\n');
        return;
    end
    
    % Parameter setting
    szDirevative = filtersize;   szWeightFunc = 5;          % Dimension of video
    szSampledDataSize = 23;                                 % Sample Image Size (Temporal)
    szHalfWindow = ceil(szWeightFunc/2);                    % Sample Image Size (Spatial)
    posImg=get(handles.posi,'Value'); 
    posX=posImg(1); posY=posImg(2);                         % The postion to observe

    % Image sample region
    topX = max(1, posX-szHalfWindow);
    bottomX = min(Height, posX+szHalfWindow);
    spanX = bottomX-topX+1;
    if topX==1
        newX = posX;
    else
        newX = posX-topX+1;
    end
    leftY = max(1, posY-szHalfWindow);
    rightY = min(Width, posY+szHalfWindow);
    spanY = rightY-leftY+1;
    if leftY==1
        newY = posY;
    else
        newY = posY-leftY+1;
    end

    % Get reference structure tensor
    ST0 = getReferenceST();
    if isempty(ST0)
        error('Cannot get reference structure tensor. Check if the video is long enough.');
    end

    data = zeros(spanX, spanY, szSampledDataSize);

    playStopped = get(handles.stop,'Value');
    bTrai = get(handles.trai,'Value');
    bDete = get(handles.dete,'Value');
    
    curPosition = szSampledDataSize + 1;
    while ((curPosition < numberOfFrames) && (~playStopped))
        if curPosition >= szSampledDataSize
            j = 1; start = curPosition-szSampledDataSize+1;
            
            % Get the images
            for i= start:curPosition
                matImage = rgb2gray(read(video,i));
                data(:,:,j) = matImage(topX:bottomX, leftY:rightY);
                j = j + 1;
            end
            matImage = imresize(matImage, [ posImagShow(4)  posImagShow(3)]);
            axes(handles.imag); imshow(matImage);    
            
            if ( (bTrai==1) || (bDete==1))
                % Compute the gradients
                [IX, IY, IT] = partial_derivative_3D(data, szDirevative);

                % Compute the images of Ixx, Ixy, Ixt, Iyy, Iyt, Itt            
                IXX=IX.*IX; IXY=IX.*IY; IXT=IX.*IT;
                            IYY=IY.*IY; IYT=IY.*IT;
                                        ITT=IT.*IT;

                % Convolve spatially each of these images with a larger Gaussian
                IXX2 = convole_3D(IXX, szWeightFunc);
                IXY2 = convole_3D(IXY, szWeightFunc);
                IXT2 = convole_3D(IXT, szWeightFunc);
                IYY2 = convole_3D(IYY, szWeightFunc);
                IYT2 = convole_3D(IYT, szWeightFunc);
                ITT2 = convole_3D(ITT, szWeightFunc);

                % Get covariance matrix, eigenvalues and distance between them
                Deigenvalues = zeros(szSampledDataSize,5); jEigvalues = 1;
                for j=1:szSampledDataSize 
                %for j=ceil(szWeightFunc/2):szSampledDataSize-floor(szWeightFunc/2)
                        % Structure Tensor
                        ST = zeros(3,3);
                        ST(1,1) = IXX2(newX,newY,j);    ST(1,2) = IXY2(newX,newY,j);  ST(1,3) = IXT2(newX,newY,j);
                        ST(2,1) = ST(1,2);              ST(2,2) = IYY2(newX,newY,j);  ST(2,3) = IYT2(newX,newY,j);
                        ST(3,1) = ST(1,3);              ST(3,2) = ST(2,3);            ST(3,3) = ITT2(newX,newY,j);

                        % Get the distance between st using generalized eigenvalue
                        [e1,e2,e3,d1,d2,d3] = eigen_decomposition(ST);%/ST0);
                        
                        [e_theta e_rho e_z] = cart2pol(e3(1), e3(2), e3(3));
                        e_phi = atan(e_z/e_rho);
                        Deigenvalues(jEigvalues,:) = [d1,d2,d3, e_theta, e_phi];
                        jEigvalues = jEigvalues + 1;
                        % dist(iFrame) = sqrt(log(d1)*log(d1)+log(d2)*log(d2)+log(d3)*log(d3));
                        fprintf(fid_results, '%.8f\t%.8f\t%.8f\n', d1, d2, d3);
                        fprintf(fid_agnles, '%.2f\t%.2f\n', e_theta, e_phi);
                end
            end
            
            % Update hisogram if training is checked
            if (bTrai && (~bDete))            
                dist = getSTDistance(Deigenvalues, bTrai);
                szDist = size(dist,2);
                
                if size(dd,2)<20000
                    distance = [dd dist]; dd = distance;
                    [nbin,xbin] = hist(distance, 50); 
                    %nbin = nbin/sum(nbin(:)).*100;
                    if isempty(xx) || (max(xx)<=max(xbin))
                        xx = xbin;
                        nn = nbin;
                    end
                else
                    for j=1:szDist
                        ind = find(xx>dist(j), 1);
                        if ~isempty(ind)
                            nn(ind(1)) = nn(ind(1)) + 1;
                        end
                    end
                    xbin = xx;
                    nbin = nn;
                end
                sum_sample = sum(nbin(:));
                fprintf('Total valid samples: %d ', sum_sample);
                nbin_p = nbin/sum_sample.*100;
                nbin_p = nbin_p;%(nbin_p>1);
                xbin_p = xbin;%(nbin_p>1);

                % Draw histogram
                if ~isempty(xbin_p)
                    bar(handles.hist, xbin_p, nbin_p, 'b');
                    xlabel(handles.hist, 'Distance'); ylabel(handles.hist,'Percentage');
                    title(handles.hist, 'Histogram of Distance Between Structure Tensors');
                end
            end
            
            % Computer and draw confidence if detecing is checked
            if (bDete && (~bTrai))
                dist = getSTDistance(Deigenvalues, bTrai);
                bNorm = get(handles.norm,'Value');
                if bNorm         % Normalize distance to [0,1]
                    if isempty(detect_distance)
                       normalizeHist();
                       nbin_p = nbin; xbin_p = xbin;
                       bar(handles.hist, xbin_p, nbin_p, 'b');
                       xlim(handles.hist, [0 1]);
                    end
                    dist = bsxfun(@rdivide, bsxfun(@minus, dist,minX), maxX-minX);
                end
                detect_distance = [detect_distance dist];
                tempConf = getConfidence(dist);
                cc = [cc tempConf]; confidence = cc;
                dConf = 0.1/(abs(occurProb(end)-tempConf)+1);
                %anomaly accumulation penalty (Commented out 10/6/2011)
                %{if size(confidence,2)>=5
                %    dConf = sum(abs(confidence(end-3:end-1) - confidence(end-4:end-2)),2)/3;
                %end
                %occurProb = [occurProb tempConf./(dConf+0.01)];
                occurProb = [occurProb tempConf+dConf*occurProb(end)];
                plot(handles.conf, 1:size(occurProb,2), occurProb);
                xlabel(handles.conf, 'Frame #'); ylabel(handles.conf, 'Occurrence Rate');  
                
                figure(11);
                plot(1:size(occurProb,2), occurProb);
                
                occurProb2 = occurProb < -3;
                figure(12);
                plot(1:size(occurProb2,2), occurProb2);
                
                figure(10);
                subplot(3,1,1), plot(detect_distance);        title('dist');
                subplot(3,1,2), plot(confidence);  title('confidence');
                subplot(3,1,3), plot(occurProb);   title('occruProb');
            end

            drawnow;
        end
        
        fprintf('Approximate position 00:%d\n', floor(curPosition/numberOfFrames*length));
        
        playStopped = get(handles.stop,'Value');
        bTrai = get(handles.trai,'Value');
        bDete = get(handles.dete,'Value');
        curPosition = curPosition + szSampledDataSize;
    end
    
    [f,xi] = ksdensity(distance);
    figure, plot(xi, f);
    
        
    fclose(fid_results);
    
    % Get reference structure tesor
    function [ST0] = getReferenceST()
      % Generate reference ST according to first ** frames        
%         if numberOfFrames < szSampledDataSize
%             ST0=[];
%             return;
%         else
%             referenceData = zeros(spanX, spanY, szSampledDataSize);
%             for ii = 1:szSampledDataSize
%                 matImage = rgb2gray(read(video,ii));
%                 referenceData(:,:,ii) = matImage(topX:bottomX, leftY:rightY);
%             end;
%             % Compute the gradients
%             [IX0, IY0, IT0] = partial_derivative_3D(referenceData);
% 
%             % Compute the images of Ixx, Ixy, Ixt, Iyy, Iyt, Itt
%             IXX0 = zeros(spanX, spanY, szSampledDataSize);
%             IXY0 = zeros(spanX, spanY, szSampledDataSize);
%             IXT0 = zeros(spanX, spanY, szSampledDataSize);
%             IYY0 = zeros(spanX, spanY, szSampledDataSize);
%             IYT0 = zeros(spanX, spanY, szSampledDataSize);
%             ITT0 = zeros(spanX, spanY, szSampledDataSize);
%             for jjFrame=1:szSampledDataSize
%                 for ii=1:spanX
%                     for jj=1:spanY
%                        IXX0(ii,jj,jjFrame) = IX0(ii,jj,jjFrame)*IX0(ii,jj,jjFrame); 
%                        IXY0(ii,jj,jjFrame) = IX0(ii,jj,jjFrame)*IY0(ii,jj,jjFrame);
%                        IXT0(ii,jj,jjFrame) = IX0(ii,jj,jjFrame)*IT0(ii,jj,jjFrame);
%                        IYY0(ii,jj,jjFrame) = IY0(ii,jj,jjFrame)*IY0(ii,jj,jjFrame);
%                        IYT0(ii,jj,jjFrame) = IY0(ii,jj,jjFrame)*IT0(ii,jj,jjFrame);
%                        ITT0(ii,jj,jjFrame) = IT0(ii,jj,jjFrame)*IT0(ii,jj,jjFrame);
%                     end
%                 end
%             end     
% 
%             % Convolve each of these images with a larger Gaussian
%             sigma = 1.0;
%             IXX0_C = convole_3D(IXX0, 'Gauss', 11, sigma);
%             IXY0_C = convole_3D(IXY0, 'Gauss', 11, sigma);
%             IXT0_C = convole_3D(IXT0, 'Gauss', 11, sigma);
%             IYY0_C = convole_3D(IYY0, 'Gauss', 11, sigma);
%             IYT0_C = convole_3D(IYT0, 'Gauss', 11, sigma);
%             ITT0_C = convole_3D(ITT0, 'Gauss', 11, sigma);  
% 
%             % Return the center point structure tensor as reference
%             ST0(1,1) = IXX0_C(newX,newY,posT);  ST0(1,2) = IXY0_C(newX,newY,posT);  ST0(1,3) = IXT0_C(newX,newY,posT);
%             ST0(2,1) = ST0(1,2);                ST0(2,2) = IYY0_C(newX,newY,posT);  ST0(2,3) = IYT0_C(newX,newY,posT);
%             ST0(3,1) = ST0(1,3);                ST0(3,2) = ST0(2,3);                ST0(3,3) = ITT0_C(newX,newY,posT);
%         end;
        ST0(1,1) = 1;                       ST0(1,2) = 0;                       ST0(1,3) = 0;
        ST0(2,1) = ST0(1,2);                ST0(2,2) = 1;                       ST0(2,3) = 0;
        ST0(3,1) = ST0(1,3);                ST0(3,2) = ST0(2,3);                ST0(3,3) = 1;
        ST0 = ST0.*1000;
    end

    % Get average confidence during #szSampledDataSize# frames
    function [cf] = getConfidence(dist)
        cf = 0.0;
        %ccNN = nbin; ccSum_sample = sum(ccNN(:));
        %ccNN_p = ccNN/ccSum_sample.*100; 
        ccNN_p = nbin;
        ccXX_p = xbin;
        deltaX = ccXX_p(2)-ccXX_p(1);
        
        [mu_dist sigma_dist] = getStat(xbin, nbin);
        
        for jj=1:szDist
            tmpInd = find(ccXX_p>dist(jj),1);
            if ((~isempty(tmpInd)) && (ccNN_p(tmpInd(1))>0))
                  cf = cf + log(ccNN_p(tmpInd(1)));
            else
                x1 = 0.0; x2 = 0.0;
                if dist(jj)<mu_dist%xbin(1)
                    x1 = xbin(1) - deltaX.*(floor((xbin(1)-dist(jj))/deltaX)+1);
                    x2 = xbin(1) - deltaX.* floor((xbin(1)-dist(jj))/deltaX);
                elseif dist(jj)>=mu_dist%xbin(end)
                    x1 = xbin(end) + deltaX.* floor((dist(jj)-xbin(end))/deltaX);
                    x2 = xbin(end) + deltaX.*(floor((dist(jj)-xbin(end))/deltaX)+1);
                else % Within [minX, maxX] but probability is zero.
                    if dist(jj)<mu_dist
                        x1 = xbin(1) - deltaX.*(floor((xbin(1)-dist(jj))/deltaX)+1);
                        x2 = xbin(1) - deltaX.* floor((xbin(1)-dist(jj))/deltaX);
                    elseif dist(jj)>=mu_dist
                        x1 = xbin(end) + deltaX.* floor((dist(jj)-xbin(end))/deltaX);
                        x2 = xbin(end) + deltaX.*(floor((dist(jj)-xbin(end))/deltaX)+1);
                    end
                end
                Eps1 = abs(x1-mu_dist);
                Eps2 = abs(x2-mu_dist);
                tmpProb = 0.5*(sigma_dist^2)*abs((1/(Eps1^2)-1/(Eps2^2)))*100; %percetage%
                cf = cf + log(tmpProb);
                disp('====Chebysheve====');
                disp(tmpProb);
            end            
        end
        cf = cf./szDist;
        
        %meanProb = mean(ccNN_p); %sum(ccNN_p/100.*ccNN_p, 2);
        %cf = (cf - meanProb)./std(ccNN_p); 
        tmp_ccNN_p = ccNN_p; tmp_ccNN_p(tmp_ccNN_p==0) = 1;
        multipProb = sum(log(tmp_ccNN_p'))/50;
        cf = cf-multipProb;
    end
    
    % Get the statistical mean and standard deviation of the distances
    function [mu,sigma] = getStat(xbin, nbin)
        xMid = xbin+ones(1,size(xbin,2))*(xbin(2)-xbin(1))/2;
        dataStat = [];
        for iStat=1:size(xbin,2)
            if nbin(iStat)>0
                dataStat = [dataStat ones(1,uint16(nbin(iStat)))*xMid(iStat)];
            end
        end
        mu = mean(dataStat,2);
        sigma = std(dataStat);
    end
    
    % Get distances of eigenvalues to their mean
    function [dist_st] = getSTDistance(dvalues, bTrain)
        nDvalues = size(dvalues,1);
        dist_st = zeros(1, nDvalues);
        if bTrain
            mmean = mean(dvalues);
            msigma = cov(dvalues);
            
            if isempty(mmu) || isempty(mcov)
                mmu = mmean;
                mcov = msigma;                
            else
                mmu = (1-malpha)*mmu + malpha*mmean;
                mcov = (1-malpha)*mcov + malpha*msigma;
            end
        end
        if rcond(mcov) > 1e-15
            for di=1:nDvalues
                dmval = dvalues(di,:) - mmu;
                dist_st(di) = sqrt(dmval/mcov*dmval');
            end
        end
    end
    
    function normalizeHist()
        minX = min(xbin);
        maxX = max(xbin);
        xbin = bsxfun(@minus, xbin, minX);
        xbin = bsxfun(@rdivide, xbin, maxX-minX);
        nbin = filter([1 2 1]/4, 1, [nbin(1) nbin nbin(end)]);
        nbin = nbin(2:end-1);
        sum_nbin = sum(nbin(:));
        nbin = nbin/sum_nbin.*100;
    end
end
        

