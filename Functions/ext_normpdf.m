function p = ext_normpdf(data, mu, sigma_)

% X = randn(10,1); mu = mean(X); SIGMA = std(X);
% p1 = normpdf(X, mu, SIGMA);
% p2 = ext_normpdf(X, mu, SIGMA); [p1, p2]
% 
% mu = [1 -1]; SIGMA = [.9 .4; .4 .3]; X = mvnrnd(mu,SIGMA,100); 
% p1 = mvnpdf(X,mu,SIGMA); 
% p2 = ext_normpdf(X, mu, SIGMA); [p1, p2]

    if size(data, 2) == 1
        p = custom_normal_pdf(data, mu, sigma_);
    else
        p = custom_multivariate_normal_pdf(data, mu, sigma_);
    end

end

function probability = custom_normal_pdf(X, mu, sigma_)
    % X: 数据矩阵，每列表示一个变量，每行表示一个观测值
    % mu: 均值
    % sigma: 方差

    % 计算概率密度
    constant_term = 1 / (sqrt(2 * pi) * sigma_);
    exponent_term = -0.5 * ((X - mu) / sigma_) .^ 2;
    probability = constant_term * exp(exponent_term);
end

function probability = custom_multivariate_normal_pdf(X, mu, sigma_)
    % X: 数据矩阵，每列表示一个变量，每行表示一个观测值
    % mu: 均值向量，长度与变量数量相同
    % sigma: 协方差矩阵，大小为 [M, M]
    
    % 检查输入的数据尺寸是否匹配
    assert(size(X, 2) == length(mu), '数据点的维度与均值向量的长度不匹配');
    assert(isequal(size(sigma_), [length(mu), length(mu)]), '协方差矩阵的尺寸不正确');
    
    % 计算概率密度
    constant_term = 1 / sqrt((2*pi)^length(mu) * det(sigma_));
    exponent_term = -0.5 * sum((X - mu) / sigma_ .* (X - mu), 2);
    probability = constant_term * exp(exponent_term);
end

function p = custom_normal_pdf_v2(X, mu, sigma_) %#ok
    % X: 数据矩阵，每列表示一个变量，每行表示一个观测值
    % mu: 均值
    % sigma: 方差
    
    % 计算标准化的数据矩阵
    X_minus_mu = bsxfun(@minus, X, mu);

    % 计算指数部分的系数
    coefficient = 1 / (sqrt(2 * pi) * sigma_);
    
    % 计算指数部分的幂次
    exponent = -0.5 * (X_minus_mu / sigma_) .^ 2 ;
    
    % 计算概率密度函数值
    p = coefficient * exp(exponent);
end

function p = custom_multivariate_normal_pdf_v2(X, mu, sigma_) %#ok
    % X: 数据矩阵，每列表示一个变量，每行表示一个观测值
    % mu: 均值向量，长度与变量数量相同
    % sigma: 协方差矩阵，大小为 [M, M]
    
    % 检查输入的数据尺寸是否匹配
    assert(size(X, 2) == length(mu), '数据点的维度与均值向量的长度不匹配');
    assert(isequal(size(sigma_), [length(mu), length(mu)]), '协方差矩阵的尺寸不正确');
    
    % 计算数据的维度
    M = length(mu);
    
    % 计算标准化的数据矩阵
    X_minus_mu = bsxfun(@minus, X, mu);
    
    % 计算协方差矩阵的逆矩阵和行列式
    inv_sigma = inv(sigma_);
    det_sigma = det(sigma_);
    
    % 计算指数部分的系数
    coefficient = (2 * pi) ^ (M/2) * sqrt(det_sigma);
    
    % 计算指数部分的幂次
    exponent = -0.5 * sum((X_minus_mu * inv_sigma) .* X_minus_mu, 2);
    
    % 计算概率密度函数值
    p = 1 / coefficient * exp(exponent);
 
end
