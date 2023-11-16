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
    % X: ���ݾ���ÿ�б�ʾһ��������ÿ�б�ʾһ���۲�ֵ
    % mu: ��ֵ
    % sigma: ����

    % ��������ܶ�
    constant_term = 1 / (sqrt(2 * pi) * sigma_);
    exponent_term = -0.5 * ((X - mu) / sigma_) .^ 2;
    probability = constant_term * exp(exponent_term);
end

function probability = custom_multivariate_normal_pdf(X, mu, sigma_)
    % X: ���ݾ���ÿ�б�ʾһ��������ÿ�б�ʾһ���۲�ֵ
    % mu: ��ֵ���������������������ͬ
    % sigma: Э������󣬴�СΪ [M, M]
    
    % �����������ݳߴ��Ƿ�ƥ��
    assert(size(X, 2) == length(mu), '���ݵ��ά�����ֵ�����ĳ��Ȳ�ƥ��');
    assert(isequal(size(sigma_), [length(mu), length(mu)]), 'Э�������ĳߴ粻��ȷ');
    
    % ��������ܶ�
    constant_term = 1 / sqrt((2*pi)^length(mu) * det(sigma_));
    exponent_term = -0.5 * sum((X - mu) / sigma_ .* (X - mu), 2);
    probability = constant_term * exp(exponent_term);
end

function p = custom_normal_pdf_v2(X, mu, sigma_) %#ok
    % X: ���ݾ���ÿ�б�ʾһ��������ÿ�б�ʾһ���۲�ֵ
    % mu: ��ֵ
    % sigma: ����
    
    % �����׼�������ݾ���
    X_minus_mu = bsxfun(@minus, X, mu);

    % ����ָ�����ֵ�ϵ��
    coefficient = 1 / (sqrt(2 * pi) * sigma_);
    
    % ����ָ�����ֵ��ݴ�
    exponent = -0.5 * (X_minus_mu / sigma_) .^ 2 ;
    
    % ��������ܶȺ���ֵ
    p = coefficient * exp(exponent);
end

function p = custom_multivariate_normal_pdf_v2(X, mu, sigma_) %#ok
    % X: ���ݾ���ÿ�б�ʾһ��������ÿ�б�ʾһ���۲�ֵ
    % mu: ��ֵ���������������������ͬ
    % sigma: Э������󣬴�СΪ [M, M]
    
    % �����������ݳߴ��Ƿ�ƥ��
    assert(size(X, 2) == length(mu), '���ݵ��ά�����ֵ�����ĳ��Ȳ�ƥ��');
    assert(isequal(size(sigma_), [length(mu), length(mu)]), 'Э�������ĳߴ粻��ȷ');
    
    % �������ݵ�ά��
    M = length(mu);
    
    % �����׼�������ݾ���
    X_minus_mu = bsxfun(@minus, X, mu);
    
    % ����Э������������������ʽ
    inv_sigma = inv(sigma_);
    det_sigma = det(sigma_);
    
    % ����ָ�����ֵ�ϵ��
    coefficient = (2 * pi) ^ (M/2) * sqrt(det_sigma);
    
    % ����ָ�����ֵ��ݴ�
    exponent = -0.5 * sum((X_minus_mu * inv_sigma) .* X_minus_mu, 2);
    
    % ��������ܶȺ���ֵ
    p = 1 / coefficient * exp(exponent);
 
end
