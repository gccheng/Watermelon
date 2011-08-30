function GMM(videofile, duration, filtersize)

%close all;

%X = st_covmat_analyse(videofile, duration);
load distance_bagleftbehind;
X = distance;
Y = X;
clear distance;

% [nbin,xbin] = hist(X, 100); 
% sizeX = size(X);
% for ii = 1:sizeX(2)
%     ind = find(xbin>Y(ii), 1);
%     if ~isempty(ind)
%         if nbin(ind)/sizeX(2)*100 < 5
%             X = X(X~=Y(ii));
%         end
%     end
% end

X = Y(Y<15);

init.mu = [0, max(X)];
init.Sigma(:,:,1) = 1;
init.Sigma(:,:,2) = 2;
init.weight = [0.5, 0.5];

fprintf('EM for Gaussian mixture: running ... \n');
R = initialization(X,init);
[~,label(1,:)] = max(R,[],2);
R = R(:,unique(label));

tol = 1e-6;
maxiter = 100;
llh = -inf(1,maxiter);
converged = false;
t = 1;
while ~converged && t < maxiter
    t = t+1;
    model = maximization(X,R);
    [R, llh(t)] = expectation(X,model);
    
    [~,label(1,:)] = max(R,[],2);
    idx = unique(label);   % non-empty components
    if size(R,2) ~= size(idx,2)
        R = R(:,idx);   % remove empty components
    else
        converged = llh(t)-llh(t-1) < tol*abs(llh(t));
    end

end
llh = llh(2:t);
if converged
    fprintf('Converged in %d steps.\n',t-1);
    u = model.mu;
    Sigma = model.Sigma;
    weight = model.weight;
    
    dist = X;
    [n,xout] = hist(dist, 100); 
    n = n/sum(n(:));
    figure, bar(xout, n, 0.95, 'y');
    ylim([0 1.0]); 
    xlabel('Distance'); ylabel('Probability');
    title('Histogram of Distance Between Structure Tensors');
    
    hold on;
    xg = 0:0.1:max(dist);
    %y1 = exp(loggausspdf(xg, u(1), Sigma(1)))./(2*pi*Sigma(1));
    y1 = mvgauss(xg, u(1), Sigma(1));
    plot(xg, y1, 'm-.');
    %y2 = exp(loggausspdf(xg, u(2), Sigma(2)))./(2*pi*Sigma(2));
    y2 = mvgauss(xg, u(2), Sigma(2));
    plot(xg, y2, 'b-.');
    
    ys = weight(1)*y1 + weight(2)*y2;
    plot(xg, ys, 'r');
    
    hold off;
    
else
    fprintf('Not converged in %d steps.\n',maxiter);
end

function R = initialization(X, init)
[d,n] = size(X);
if isstruct(init)  % initialize with a model
    R  = expectation(X,init);
elseif length(init) == 1  % random initialization
    k = init;
    idx = randsample(n,k);
    m = X(:,idx);
    [~,label] = max(bsxfun(@minus,m'*X,sum(m.^2,1)'/2),[],1);
    while k ~= unique(label)
        idx = randsample(n,k);
        m = X(:,idx);
        [~,label] = max(bsxfun(@minus,m'*X,sum(m.^2,1)'/2),[],1);
    end
    R = full(sparse(1:n,label,1,n,k,n));
elseif size(init,1) == 1 && size(init,2) == n  % initialize with labels
    label = init;
    k = max(label);
    R = full(sparse(1:n,label,1,n,k,n));
elseif size(init,1) == d  %initialize with only centers
    k = size(init,2);
    m = init;
    [~,label] = max(bsxfun(@minus,m'*X,sum(m.^2,1)'/2),[],1);
    R = full(sparse(1:n,label,1,n,k,n));
else
    error('ERROR: init is not valid.');
end

function [R, llh] = expectation(X, model)
mu = model.mu;
Sigma = model.Sigma;
w = model.weight;

n = size(X,2);
k = size(mu,2);
logR = zeros(n,k);

for i = 1:k
    logR(:,i) = loggausspdf(X,mu(:,i),Sigma(:,:,i));
end
logR = bsxfun(@plus,logR,log(w));
T = logsumexp(logR,2);
llh = sum(T)/n; % loglikelihood
logR = bsxfun(@minus,logR,T);
R = exp(logR);


function model = maximization(X, R)
[d,n] = size(X);
k = size(R,2);

s = sum(R,1);
w = s/n;
mu = bsxfun(@times, X*R, 1./s);
Sigma = zeros(d,d,k);
for i = 1:k
    Xo = bsxfun(@minus,X,mu(:,i));
    Xo = bsxfun(@times,Xo,sqrt(R(:,i)'));
    Sigma(:,:,i) = Xo*Xo'/s(i);
    Sigma(:,:,i) = Sigma(:,:,i)+eye(d)*(1e-6); % add a prior for numerical stability
end

model.mu = mu;
model.Sigma = Sigma;
model.weight = w;

function y = loggausspdf(X, mu, Sigma)
d = size(X,1);
X = bsxfun(@minus,X,mu);
[R,p]= chol(Sigma);
if p ~= 0
    error('ERROR: Sigma is not PD.');
end
q = sum((R'\X).^2,1);  % quadratic term (M distance)
c = d*log(2*pi)+2*sum(log(diag(R)));   % normalization constant
y = -(c+q)/2;


function y = mvgauss(X, mu, Sigma)
[d n] = size(X);
X = bsxfun(@minus, X, mu);
y = zeros(1, n);
for i=1:n
    q = exp(-1/2*X(i)'/Sigma*X(i));
    c = 1/(abs(2*pi)^(d/2)*det(Sigma)^(1/2));
    y(i) = c*q;
end

