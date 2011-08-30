function [ distance xbin nbin confidence] = ...
    st_covmat_analyse(videofile, duration, filtersize, handles, conf, d, x, n)

    distance = []; xbin =[]; nbin = []; confidence = [];
    xx = x; nn = n; dd = d; cc = conf;
    posImagShow = get(handles.imag, 'Position');

    % Video information
    video = mmreader(videofile);
    read(video,inf);
    numberOfFrames = video.NumberOfFrames;
    Height = video.Height;
    Width = video.Width;
    length = video.duration;

    % Parameter setting
    szDirevative = 5;   szWeightFunc = 11;          % Dimension of video
    szSampledDataSize = szWeightFunc*4+1;           % Sample Image Size (Temporal)
    szHalfWindow = floor(szSampledDataSize/2);
    nLength = szWeightFunc; nFrames=nLength;        % Total length of video 
    posImg=get(handles.posi,'Value'); 
    posX=posImg(1); posY=posImg(2); posT = szWeightFunc;           % The postion to observe

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
    PlayerControl = get(handles.figure1, 'UserData');

    playStopped = get(handles.stop,'Value');
    bTrai = get(handles.trai,'Value');
    bDete = get(handles.dete,'Value');
    
    curPosition = szSampledDataSize + 1;
    while (curPosition < numberOfFrames) && (~playStopped) && (bTrai==1 || bDete==1) 
        
        if curPosition >= szSampledDataSize
            j = 1; start = curPosition-szSampledDataSize+1;
            
            for i= start:curPosition
                matImage = rgb2gray(read(video,i));
                data(:,:,j) = matImage(topX:bottomX, leftY:rightY);
                j = j + 1;
            end
            newImag = imresize(read(video,curPosition), [ posImagShow(4)  posImagShow(3)]);
            axes(handles.imag); imshow(newImag);    
            
            % Compute the gradients
            [IX, IY, IT] = partial_derivative_3D(data);

            % Compute the images of Ixx, Ixy, Ixt, Iyy, Iyt, Itt
            for jFrame=1:szSampledDataSize
                for i=1:spanX
                    for j=1:spanY
                       IXX(i,j,jFrame) = IX(i,j,jFrame)*IX(i,j,jFrame); 
                       IXY(i,j,jFrame) = IX(i,j,jFrame)*IY(i,j,jFrame);
                       IXT(i,j,jFrame) = IX(i,j,jFrame)*IT(i,j,jFrame);
                       IYY(i,j,jFrame) = IY(i,j,jFrame)*IY(i,j,jFrame);
                       IYT(i,j,jFrame) = IY(i,j,jFrame)*IT(i,j,jFrame);
                       ITT(i,j,jFrame) = IT(i,j,jFrame)*IT(i,j,jFrame);
                    end
                end
            end     

            % Convolve each of these images with a larger Gaussian
            sigma = 1.0;
            IXX2 = convole_3D(IXX, 'Gauss', 11, sigma);
            IXY2 = convole_3D(IXY, 'Gauss', 11, sigma);
            IXT2 = convole_3D(IXT, 'Gauss', 11, sigma);
            IYY2 = convole_3D(IYY, 'Gauss', 11, sigma);
            IYT2 = convole_3D(IYT, 'Gauss', 11, sigma);
            ITT2 = convole_3D(ITT, 'Gauss', 11, sigma);

            % Get covariance matrix, eigenvalues and distance between them
            dist = zeros(1,nFrames); iFrame = 1;
            for j=1:szSampledDataSize %for j=ceil(szWeightFunc/2):szSampledDataSize-floor(szWeightFunc/2)
                    % Structure Tensor
                    ST = zeros(3,3);
                    ST(1,1) = IXX2(newX,newY,j);    ST(1,2) = IXY2(newX,newY,j);  ST(1,3) = IXT2(newX,newY,j);
                    ST(2,1) = ST(1,2);              ST(2,2) = IYY2(newX,newY,j);  ST(2,3) = IYT2(newX,newY,j);
                    ST(3,1) = ST(1,3);              ST(3,2) = ST(2,3);            ST(3,3) = ITT2(newX,newY,j);

                    % Get the distance between st using generalized eigenvalue
                    [e1,e2,e3,d1,d2,d3] = eigen_decomposition(ST/ST0);
                    dist(iFrame) = sqrt(log(d1)*log(d1)+log(d2)*log(d2)+log(d3)*log(d3));
                    iFrame = iFrame + 1;
            end
            dist = dist(~isinf(dist));
            szDist = size(dist,2);
            
            % Update hisogram if training is checked
            if (bTrai)
                if size(dd,2)<2000
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
                nbin_p = nbin_p(nbin_p>1);
                xbin_p = xbin(nbin_p>1);

                % Draw histogram
                if ~isempty(xbin_p)
                    bar(handles.hist, xbin_p, nbin_p, 'b');
                    xlabel(handles.hist, 'Distance'); ylabel(handles.hist,'Percentage');
                    title(handles.hist, 'Histogram of Distance Between Structure Tensors');
                end
            end
            
            % Computer and draw confidence if detecing is checked
            if (bDete==1) && (size(dd,2)>500)
                tempConf = getConfidence();
                cc = [cc tempConf]; confidence = cc;
                plot(handles.conf, 1:size(cc,2), cc);
                xlabel(handles.conf, 'Window Frame #'); ylabel(handles.conf, 'Normality Confidence');
            end

            drawnow;
        end
        
        fprintf('Approximate position 00:%d\n', floor(curPosition/numberOfFrames*length));
        
        playStopped = get(handles.stop,'Value');
        bTrai = get(handles.trai,'Value');
        bDete = get(handles.dete,'Value');
        curPosition = curPosition + szSampledDataSize;
    end
    
    % Get reference structure tesor
    function [ST0] = getReferenceST()
%       % Generate reference ST according to first ** frames        
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
        ST0 = ST0.*100000;
    end

    % Get average confidence during #szSampledDataSize# frames
    function [cf] = getConfidence()
        cf = 0.0;
        ccNN = nbin; ccXX = xbin;
        ccSum_sample = sum(ccNN(:));
        ccNN_p = ccNN/ccSum_sample.*100; 
        %ccNN_p = ccNN(ccNN>1);
        %ccXX_p = ccXX(ccNN>1);
        ccXX_p = ccXX;
        
        avgProb = sum((ccNN_p./100).*ccNN_p, 2);
        covProb = std(ccNN_p);
        for jj=1:szDist
            tmpInd = find(ccXX_p>dist(jj),1);
            if ~isempty(tmpInd)
                if tmpInd(1)>1
                    cf = cf + ccNN_p(tmpInd(1))-avgProb;
                else
                    cf = cf + 0 - avgProb;
                end
                %cf = cf.*exp(-cf/10000);
            end
        end
        cf = cf./szDist; cf = cf./(2*covProb);
    end
    
    function [mu,sigma] = getStat(xbin, nbin)
        xMid = xbin-ones(1,size(xbin,2))*(xbin(2)-xbin(1));
        dataStat = [];
        for i=1:size(xbin,2)
            dataStat = [dataStat ones(1,nbin(i))*xMid(i)];
        end
        mu = mean(dataStat,2);
        sigma = std(dataStat);
    end

end
        

