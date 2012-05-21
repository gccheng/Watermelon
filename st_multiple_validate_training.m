function [ distance xbin nbin confidence] = st_multiple_validate_training(videodir, ...
        filtersize, handles, HistTrained)
    
    picturelist = get_files(videodir);
      
    % Parameter setting
    szDerivative = filtersize;   szConvFunc = 5;             % Derivative and Convolution kernel size
    szSampledTemporalSize = 23;                              % Sample Size (Temporal)
    szSampledSpatialSize = 2*szConvFunc+1;                   % Sample Size (Spatial)   
    %abbreviation
    ts = szSampledTemporalSize;
    ss = szSampledSpatialSize;
    
    % Video information
    numberOfFrames = size(picturelist, 1);
    readNFrames(1,1);
    [Height, Width, ~] = size(qFrames); 
    
    distance = cell(floor(Height/szSampledSpatialSize), floor(Width/szSampledSpatialSize));
    confidence = zeros(floor(Height/szSampledSpatialSize),floor(Width/szSampledSpatialSize));
    if (0==HistTrained)
        Histograms = zeros(52,floor(Height/szSampledSpatialSize),floor(Width/szSampledSpatialSize));
        xbin =zeros(50, floor(Height/szSampledSpatialSize), floor(Width/szSampledSpatialSize));
        nbin = zeros(50, floor(Height/szSampledSpatialSize), floor(Width/szSampledSpatialSize));
        
        % Add 9/29/2011 Mahalanobis distance
        mmu = zeros(3,floor(Height/szSampledSpatialSize),floor(Width/szSampledSpatialSize));
        mcov = zeros(3,3,floor(Height/szSampledSpatialSize),floor(Width/szSampledSpatialSize));
        minX = zeros(floor(Height/szSampledSpatialSize),floor(Width/szSampledSpatialSize));
        maxX = ones(floor(Height/szSampledSpatialSize),floor(Width/szSampledSpatialSize))*Inf;
        malpha = ones(floor(Height/szSampledSpatialSize),floor(Width/szSampledSpatialSize))*0.1;
    end
    qFrames = zeros(Height, Width, szSampledTemporalSize);   
    
    % Get reference structure tensor
    % ST0 = getReferenceST();

    playStopped = get(handles.stop,'Value');
    curPosition = 1; 
    while (curPosition+szSampledTemporalSize < numberOfFrames) && (~playStopped)
        % read szSampledTemporalSize frames and show current frame
        readNFrames(curPosition, szSampledTemporalSize);
        curPosition = curPosition + szSampledTemporalSize;
        
        newImag = squeeze(qFrames(:,:,1));
        axes(handles.imag); imshow(newImag,[min(newImag(:)) max(newImag(:))]);
        xlabel('video frame/23');
        fprintf('Current Position:%d/%d\n',curPosition,numberOfFrames);
                
        % Compute the gradients
        [IX, IY, IT] = partial_derivative_3D(qFrames, szDerivative);
        
        % Compute the images of Ixx, Ixy, Ixt, Iyy, Iyt, Itt
        IXX=IX.*IX; IXY=IX.*IY; IXT=IX.*IT;     clear IX;
        IYY=IY.*IY; IYT=IY.*IT;     clear IY;
        ITT=IT.*IT;     clear IT;
        
        % Convolve spatially each of these images with a larger Gaussian
        % GH = importdata('DGKernel.mat');
        IXX2 = convole_3D(IXX, szConvFunc); clear IXX;
        IXY2 = convole_3D(IXY, szConvFunc); clear IXY;
        IXT2 = convole_3D(IXT, szConvFunc); clear IXT;
        IYY2 = convole_3D(IYY, szConvFunc); clear IYY;
        IYT2 = convole_3D(IYT, szConvFunc); clear IYT;
        ITT2 = convole_3D(ITT, szConvFunc); clear ITT;
        
        for posX=szSampledSpatialSize:szSampledSpatialSize:Height
            for posY=szSampledSpatialSize:szSampledSpatialSize:Width
                indX = floor(posX/szSampledSpatialSize); indY = floor(posY/szSampledSpatialSize);
                % Get covariance matrix, eigenvalues and distance between them
                Deigenvalues = zeros(szSampledTemporalSize,3); jEigvalues = 1;
                for j=1:szSampledTemporalSize
                    ST = zeros(3,3);
                    ST(1,1) = IXX2(posX,posY,j);    ST(1,2) = IXY2(posX,posY,j);  ST(1,3) = IXT2(posX,posY,j);
                    ST(2,1) = ST(1,2);              ST(2,2) = IYY2(posX,posY,j);  ST(2,3) = IYT2(posX,posY,j);
                    ST(3,1) = ST(1,3);              ST(3,2) = ST(2,3);            ST(3,3) = ITT2(posX,posY,j);
                    
                    % Get the distance between st using generalized eigenvalue
                    [~,~,~,d1,d2,d3] = eigen_decomposition(ST);
                    Deigenvalues(jEigvalues,:) = [d1,d2,d3];
                    jEigvalues = jEigvalues + 1;
                end
                
                % compute and save distance
                dist = getSTDistance(Deigenvalues, 1, indX, indY);
                distance{indX,indY} = [distance{indX,indY} dist];
                %[nTmp,xTmp] = hist(distance{indX,indY}, 50);
                %nbin(:,indX,indY) = nTmp; xbin(:,indX,indY)=xTmp;
            end
        end

        % proceed to process other frames
        playStopped = get(handles.stop,'Value');
    end % while

    % normalize histograms
    for posX=szSampledSpatialSize:szSampledSpatialSize:Height
        for posY=szSampledSpatialSize:szSampledSpatialSize:Width
            indX = floor(posX/szSampledSpatialSize); indY = floor(posY/szSampledSpatialSize);
            [nTmp,xTmp] = hist(distance{indX,indY}, 50);
            nbin(:,indX,indY) = nTmp; xbin(:,indX,indY)=xTmp;
            normalizeHist(indX,indY);
        end;
    end;

    delete('Histograms.mat'); delete('mmu.mat'); delete('mcov.mat'); delete('distance.mat');
    save('distance.mat', 'distance');
    save('Histograms.mat', 'Histograms');
    save('mmu.mat', 'mmu');
    save('mcov.mat', 'mcov');
    
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
        dist_st = dist_st(dist_st~=0);
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

    function readNFrames(curPosition, N)
        for pos=0:N-1
            %sFile = sprintf('%03d.tif', curPosition+pos);
            sFile = picturelist{curPosition+pos};
            pTmp = imread(sFile);
            if ndims(pTmp)==3
                pTmp = rgb2gray(pTmp);
            end;
            qFrames(:,:,szSampledTemporalSize-N+1+pos) = pTmp;
        end
    end

    function readNFramesVideo(curPosition, N)
        videoreader = mmreader('C:\\Users\\Guangchun\\Dropbox\\Videos\\LeftBag_xvid.avi');
        for pos=0:N-1
            tmpFrame = rgb2gray(read(videoreader, curPosition+pos));
            qFrames(:,:,szSampledTemporalSize-N+1+pos) = tmpFrame;            
        end
    end

    function frame = readNFramesDir(curPosition, N)
        for pos=0:N-1
            sFile = sprintf('%04d.tif', curPosition+pos);
            qFrames(:,:,szSampledTemporalSize-N+1+pos) = imread([videodir '/' sFile]);
        end
        frame = qFrames;
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