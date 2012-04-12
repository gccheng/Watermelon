function [ distance xbin nbin confidence] = st_multiple_validate(videodir, ...
        filtersize, handles, HistTrained)
    
    cla(handles.precision_recall);
    P = []; % precision
    R = []; % recall
       
    picturelist = get_files(videodir);
    gtlist = get_files(get(handles.eb_groundtruth, 'String'));
    
    % Parameter setting
    szDerivative = filtersize;   szConvFunc = 5;             % Derivative and Convolution kernel size
    szSampledTemporalSize = 23;                              % Sample Size (Temporal)
    szSampledSpatialSize = szConvFunc+1;                   % Sample Size (Spatial)   
    %abbreviation
    ts = szSampledTemporalSize;
    ss = szSampledSpatialSize;
    
    % Video information
    curPosition = 1;
    numberOfFrames = size(dir(videodir), 1)-2;
    firstFrame = imread(picturelist{1});
    [Height, Width, ~] = size(firstFrame); 
    
    distance = cell(floor(Height/szSampledSpatialSize), floor(Width/szSampledSpatialSize));
    confidence = zeros(floor(Height/szSampledSpatialSize),floor(Width/szSampledSpatialSize));      
    qFrames = zeros(Height, Width, szSampledTemporalSize);
    qDist = zeros(floor(Height/szSampledSpatialSize), floor(Width/szSampledSpatialSize), szSampledTemporalSize);
    
    % weighting scheme of probability within the window
    Mw = floor(szSampledTemporalSize/2+1);
    Sw = 11;
    Xw = 1:23;
    Wt = mvnpdf(Xw', Mw, Sw);
    
    if (1==HistTrained)
        Histograms = importdata('Histograms.mat');
        [tB tX tY] = size(Histograms);
        mmu = importdata('mmu.mat');
        mcov = importdata('mcov.mat');
        xbin = permute(repmat(repmat(0:0.1:0.9, [tX,1]), [1,1,tY]), [2 1 3]);
        nbin = Histograms(3:end,:,:);
        minX = squeeze(Histograms(1,:,:));
        maxX = squeeze(Histograms(2,:,:));
        malpha = ones(floor(Height/szSampledSpatialSize),floor(Width/szSampledSpatialSize))*0.1;
    end
            
    % Get reference structure tensor
    % ST0 = getReferenceST();
    
    playStopped = get(handles.stop,'Value');
    bTrai = 0;
    
    curPosition = 1;
    flag = 1;
    while (curPosition+1 < numberOfFrames) && (playStopped ~= 1.0)       
        if (1==flag) % first szSampledTemporalSize frames
            % load frames and show current frame
            readNFrames(curPosition, szSampledTemporalSize);
            curPosition = curPosition + szSampledTemporalSize;
            
            flag = 0;
            [qIX, qIY, qIT] = partial_derivative_3D(qFrames, szDerivative);
            
            % Compute the images of Ixx, Ixy, Ixt, Iyy, Iyt, Itt
            qIXX=qIX.*qIX;  qIXY=qIX.*qIY;  qIXT=qIX.*qIT;
                            qIYY=qIY.*qIY;  qIYT=qIY.*qIT;
                                            qITT=qIT.*qIT;
            
            % Convolve spatially each of these images with a larger Gaussian
            % GH = importdata('DGKernel.mat');
            qIXX2 = convole_3D(qIXX, szConvFunc);
            qIXY2 = convole_3D(qIXY, szConvFunc);
            qIXT2 = convole_3D(qIXT, szConvFunc);
            qIYY2 = convole_3D(qIYY, szConvFunc);
            qIYT2 = convole_3D(qIYT, szConvFunc);
            qITT2 = convole_3D(qITT, szConvFunc);
            
            for posX=szSampledSpatialSize:szSampledSpatialSize:Height
                for posY=szSampledSpatialSize:szSampledSpatialSize:Width
                    indX = floor(posX/szSampledSpatialSize); indY = floor(posY/szSampledSpatialSize);
                    % Get covariance matrix, eigenvalues and distance between them
                    Deigenvalues = zeros(szSampledTemporalSize,3); jEigvalues = 1;
                    for j=1:szSampledTemporalSize
                        ST = zeros(3,3);
                        ST(1,1) = qIXX2(posX,posY,j);    ST(1,2) = qIXY2(posX,posY,j);  ST(1,3) = qIXT2(posX,posY,j);
                        ST(2,1) = ST(1,2);               ST(2,2) = qIYY2(posX,posY,j);  ST(2,3) = qIYT2(posX,posY,j);
                        ST(3,1) = ST(1,3);               ST(3,2) = ST(2,3);            ST(3,3) = qITT2(posX,posY,j);
                        
                        % Get the distance between st using generalized eigenvalue
                        [~,~,~,d1,d2,d3] = eigen_decomposition(ST);
                        Deigenvalues(jEigvalues,:) = [d1,d2,d3];
                        jEigvalues = jEigvalues + 1;
                    end
                    qDist(indX, indY, :) = getSTDistance(Deigenvalues, bTrai, indX, indY);
                end
            end    
        else
            readNFrames(curPosition, 1);
            curPosition = curPosition + 1;
            
            % Remove the 1st, update the last 4, add 1
            qIX(:,:,1) = []; qIY(:,:,1) = []; qIT(:,:,1) = [];
            [qIX_3, qIY_3, qIT_3] = partial_derivative_3D(qFrames(:,:,ts-2:ts), szDerivative);
            qIX(:,:,ts-1:ts) = qIX_3(:,:,2:3);
            qIY(:,:,ts-1:ts) = qIY_3(:,:,2:3);
            qIT(:,:,ts-1:ts) = qIT_3(:,:,2:3);
            
            % Compute the images of Ixx, Ixy, Ixt, Iyy, Iyt, Itt
            qIXX=qIX.*qIX; qIXY=qIX.*qIY; qIXT=qIX.*qIT;
            qIYY=qIY.*qIY; qIYT=qIY.*qIT;
            qITT=qIT.*qIT;
            
            % Convolve spatially each of these images with a larger Gaussian
            qIXX2(:,:,1) = []; qIXY2(:,:,1) = []; qIXT2(:,:,1) = [];
            qIYY2(:,:,1) = []; qIYT2(:,:,1) = []; qITT2(:,:,1) = [];
            qIXX2_5 = convole_3D(qIXX(:,:,ts-4:ts), szConvFunc);
            qIXY2_5 = convole_3D(qIXY(:,:,ts-4:ts), szConvFunc);
            qIXT2_5 = convole_3D(qIXT(:,:,ts-4:ts), szConvFunc);
            qIYY2_5 = convole_3D(qIYY(:,:,ts-4:ts), szConvFunc);
            qIYT2_5 = convole_3D(qIYT(:,:,ts-4:ts), szConvFunc);
            qITT2_5 = convole_3D(qITT(:,:,ts-4:ts), szConvFunc);
            qIXX2(:,:,ts-2:ts) = qIXX2_5(:,:,3:5);
            qIXY2(:,:,ts-2:ts) = qIXY2_5(:,:,3:5);
            qIXT2(:,:,ts-2:ts) = qIXT2_5(:,:,3:5);
            qIYY2(:,:,ts-2:ts) = qIYY2_5(:,:,3:5);
            qIYT2(:,:,ts-2:ts) = qIYT2_5(:,:,3:5);
            qITT2(:,:,ts-2:ts) = qITT2_5(:,:,3:5);
            
            % Only the last 3 frames are updated, + 1st frame removed
            qDist(:, :, 1) = [];
            for posX=szSampledSpatialSize:szSampledSpatialSize:Height
                for posY=szSampledSpatialSize:szSampledSpatialSize:Width
                    indX = floor(posX/szSampledSpatialSize); indY = floor(posY/szSampledSpatialSize);
                    % Get covariance matrix, eigenvalues and distance between them
                    Deigenvalues = zeros(3,3); jEigvalues = 1;
                    for j=ts-2:ts
                        ST = zeros(3,3);
                        ST(1,1) = qIXX2(posX,posY,j);    ST(1,2) = qIXY2(posX,posY,j);  ST(1,3) = qIXT2(posX,posY,j);
                        ST(2,1) = ST(1,2);               ST(2,2) = qIYY2(posX,posY,j);  ST(2,3) = qIYT2(posX,posY,j);
                        ST(3,1) = ST(1,3);               ST(3,2) = ST(2,3);            ST(3,3) = qITT2(posX,posY,j);
                        
                        % Get the distance between st using generalized eigenvalue
                        [~,~,~,d1,d2,d3] = eigen_decomposition(ST);
                        Deigenvalues(jEigvalues,:) = [d1,d2,d3];
                        jEigvalues = jEigvalues + 1;
                    end
                    qDist(indX, indY, ts-2:ts) = getSTDistance(Deigenvalues, bTrai, indX, indY);
                end
            end
            
        end
        
        newImage = squeeze(qFrames(:,:,1));
        axes(handles.imag); imshow(newImage, [min(newImage(:))  max(newImage(:))]);
        fprintf('Current Position:%d/%d\n',curPosition,numberOfFrames);
        
        for indX=1:floor(Height/szSampledSpatialSize)
            for indY=1:floor(Width/szSampledSpatialSize)
                if sum(Histograms(:,indX,indY))==0
                    normalizeHist(indX,indY);
                end
                dist = bsxfun(@rdivide, bsxfun(@minus, qDist(indX,indY,:),minX(indX,indY)), maxX(indX,indY)-minX(indX,indY));
                tempConf = getConfidence(dist,indX,indY);
                incCoff = 0.2/(abs(confidence(indX,indY)-tempConf)+1);
                confidence(indX,indY)=confidence(indX,indY)*incCoff + tempConf;
            end
        end
        
        theta = 3.5;
        confidenceShow = -confidence;
        %H = fspecial('average', 3); confidenceShow = imfilter(confidenceShow, H, 'replicate');
        confidenceShow(confidenceShow<theta) = 0;
        confidenceShow(confidenceShow>=theta) = 1;
        if(sum(confidenceShow(:))<2)
            confidenceShow = zeros(size(confidenceShow));
        end;
        
        confidenceShow = imresize(confidenceShow,[Height Width],'bilinear');
        confidenceShow(confidenceShow>0.7) = 1;
        confidenceShow(confidenceShow<=0.7) = 0;
        
        figure(1);
        imshow(confidenceShow, [min(confidenceShow(:)) max(confidenceShow(:))]);
        %axes(handles.precision_recall); imshow(confidenceShow, [min(confidenceShow(:)) max(confidenceShow(:))]);
        
        [p, r] = computePrecisionRecall(uint8(confidenceShow), curPosition);
        P = [P, p];
        R = [R, r];
        
        axes(handles.roc); plot(P, R, '+');
       
        playStopped = get(handles.stop,'Value');
        
    end % while
    
    disp(mean(P));
    disp(mean(R));
    
    function [precision, recall] = computePrecisionRecall(anomaly_detect, pos)   
        % get the corresponding groundtruth
        prec = zeros(1, szSampledTemporalSize);
        rec = zeros(1, szSampledTemporalSize);
        for i=1:szSampledTemporalSize
            gt = imread(gtlist{pos-i});
            gt(gt<theta) = 0;
            gt(gt>=theta) = 1;

            % compare and compute
            umat = anomaly_detect.*gt;
            uval = sum(umat(:));
            truthval = sum(gt(:));
            detectval = sum(anomaly_detect(:));
            prec(szSampledTemporalSize-i+1) = (uval+0.01)/(detectval+0.01);
            rec(szSampledTemporalSize-i+1) = (uval+0.01)/(truthval+0.01);
        end
        
        gt = imread(gtlist{pos});
        gt(gt<theta) = 0; gt(gt>=theta) = 1;
        figure(11);  imshow(gt, [min(gt(:)) max(gt(:))]);
        
        axes(handles.precision_recall);
        plot(prec, 'o');
        
        precision = max(prec(:));
        recall = max(rec(:));
    end
    
    % Get reference structure tesor
    function [ST0] = getReferenceST()
        ST0(1,1) = 1;                       ST0(1,2) = 0;                       ST0(1,3) = 0;
        ST0(2,1) = ST0(1,2);                ST0(2,2) = 1;                       ST0(2,3) = 0;
        ST0(3,1) = ST0(1,3);                ST0(3,2) = ST0(2,3);                ST0(3,3) = 1;
        ST0 = ST0.*1000;
    end

    % Get average confidence during #szSampledTemporalSize# frames
    function [cf] = getConfidence(dist,indX,indY)
        cf = 0.0; szdist = size(dist,2);
        ccNN_p = Histograms(3:end, indX, indY);
        ccXX_p = xbin(:,indX,indY);
        deltaX = ccXX_p(2)-ccXX_p(1);
        
        [mu_dist sigma_dist] = getStat(ccXX_p, ccNN_p);
        if (isnan(sigma_dist))
            return;
        end;
        
        for jj=1:szdist
            tmpInd = find(ccXX_p>dist(jj),1);
            if ((~isempty(tmpInd)) && (tmpInd(1)>1) && (ccNN_p(tmpInd(1))>0))
                  cf = cf + log(ccNN_p(tmpInd(1)));
            else
                x1 = 0.0; x2 = 0.0;
                if dist(jj)<xbin(1,indX,indY)
                    x1 = xbin(1,indX,indY) - deltaX.*(floor((xbin(1,indX,indY)-dist(jj))/deltaX)+1);
                    x2 = xbin(1,indX,indY) - deltaX.* floor((xbin(1,indX,indY)-dist(jj))/deltaX);
                elseif dist(jj)>xbin(end,indX,indY)
                    x1 = xbin(end,indX,indY) + deltaX.* floor((dist(jj)-xbin(end,indX,indY))/deltaX);
                    x2 = xbin(end,indX,indY) + deltaX.*(floor((dist(jj)-xbin(end,indX,indY))/deltaX)+1);
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
                % disp('====Chebysheve===='); disp(tmpProb);
            end            
        end
        cf = cf./szdist;
        
        %meanProb = mean(ccNN_p); %sum(ccNN_p/100.*ccNN_p, 2);
        %cf = (cf - meanProb)./std(ccNN_p); 
        tmp_ccNN_p = ccNN_p; tmp_ccNN_p(tmp_ccNN_p==0) = 1;
        multipProb = sum(log(tmp_ccNN_p'))/50;
        cf = cf - multipProb;
    end
    
    % Get the statistical mean and standard deviation of the distances
    function [mu,sigma] = getStat(xbin_xy, nbin_xy)
        xMid = xbin_xy+ones(size(xbin_xy,1),1)*(xbin_xy(2)-xbin_xy(1))/2;
        dataStat = [];
        for iStat=1:size(xbin_xy,1)
            if nbin_xy(iStat)>0
                dataStat = [dataStat ones(1,uint16(nbin_xy(iStat)))*xMid(iStat)];
            end
        end
        mu = mean(dataStat,2);
        sigma = std(dataStat);
    end
    
    % Get distances of eigenvalues to their mean
    function [dist_st] = getSTDistance(dvalues, bTrain, indX, indY)
        nDvalues = size(dvalues,1);
        dist_st = zeros(1, nDvalues);
        if bTrain
            mmean = mean(dvalues);
            msigma = cov(dvalues);
            
            if (sum(mmu(:,indX,indY))==0) || (sum(sum(mcov(:,:,indX,indY)))==0)
                mmu(:,indX,indY) = mmean';
                mcov(:,:,indX,indY) = msigma;                
            else
                mmu(:,indX,indY) = (1-malpha(indX,indY))*mmu(:,indX,indY)...
                                + malpha(indX,indY)*mmean';
                mcov(:,:,indX,indY) = (1-malpha(indX,indY))*mcov(:,:,indX,indY)...
                                + malpha(indX,indY)*msigma;
            end
        end
        if rcond(mcov(:,:,indX,indY)) > 1e-15
            for di=1:nDvalues
                dmval = dvalues(di,:) - mmu(:,indX,indY)';
                dist_st(di) = sqrt(dmval/mcov(:,:,indX,indY)*dmval');
            end
        end
        %dist_st = dist_st(dist_st~=0);
    end
    
    function normalizeHist(indX,indY)
        minX(indX,indY) = min(xbin(:,indX,indY));
        maxX(indX,indY) = max(xbin(:,indX,indY));
        xbin(:,indX,indY) = bsxfun(@minus, xbin(:,indX,indY), minX(indX,indY));
        xbin(:,indX,indY) = bsxfun(@rdivide, xbin(:,indX,indY), maxX(indX,indY)-minX(indX,indY));
        Histograms(1:2,indX,indY) = [minX(indX,indY);maxX(indX,indY)];
        tmpHistogram = filter([1 2 1]/4, 1, [nbin(1,indX,indY);nbin(:,indX,indY);nbin(end,indX,indY)]);
        sum_tmpHistogram = sum(tmpHistogram(:));
        tmpHistogram = tmpHistogram/sum_tmpHistogram.*100;
        Histograms(3:end,indX,indY) = tmpHistogram(2:end-1);
    end    

    function readNFrames(curPosition,N)
        if curPosition>1
            for i=1:N
                qFrames(:,:,1) = [];
            end;
        end;
        for pos=0:N-1
            sFile = picturelist{curPosition+pos};
            qFrames(:,:,szSampledTemporalSize-N+1+pos) = rgb2gray(imread(sFile));
        end
    end

    function readNFramesVideo(curPosition,N)
        videoreader = mmreader('C:\\Users\\Guangchun\\Dropbox\\Videos\\LeftBag_xvid.avi');
        for pos=0:N-1
            tmpFrame = rgb2gray(read(videoreader, curPosition+pos));
            qFrames(:,:,23-N+1+pos) = tmpFrame;            
        end
    end

    function readNFramesDir(curPosition,N)
        for pos=0:N-1
            sFile = sprintf('%04d.tif', curPosition+pos);
            qFrames(:,:,szSampledTemporalSize-N+1+pos) = imread([videodir '/' sFile]);
        end
    end

    % Get all the pictures under a directory
    function [picturelist] = get_files(videodir)
        listing = dir(videodir);
        count = size(listing, 1);
        
        picturelist = cell(0);
        
        % get all valid file names
        for i=1:count
            name = listing(i).name;
            isdir = listing(i).isdir;
            if ((~isdir) && (name(1)~='.'))
                picturelist = [picturelist; [videodir, '/', name]];
            elseif (isdir && (~strcmp('.',name)) && (~strcmp('..',name)))
                sub_files = files_in_folder([videodir '/', name]);
                picturelist = [picturelist; sub_files];
            end;
        end;
    end

    function [files] = files_in_folder(folder)
        files = cell(0);
        listing = dir(folder);
        count = size(listing, 1);
        for i=1:count
            name = listing(i).name;
            isdir = listing(i).isdir;
            if ((~isdir) && (name(1)~='.'))
                files = [files; [folder '/', name]];
            end;
        end;
    end

end

