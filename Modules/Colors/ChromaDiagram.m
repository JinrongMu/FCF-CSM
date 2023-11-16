classdef ChromaDiagram  < handle 
    properties
        cieType     % CIE����
        waveLength  % ����
        x           % ɫ������ x
        y           % ɫ������ y
    end
    
    methods
        function obj = ChromaDiagram(cieType)
            % ���캯��
            if nargin < 1
                cieType = 'CIE1931';  % Ĭ��ֵΪ 'CIE1931'
            end
            obj.cieType = cieType;
        end
        
        function obj = initCoordinates(obj)
            % ��ʼ��ɫ������
            contourData = load('CIE1931.txt');
            obj.waveLength = contourData(:, 1);
            x_ = contourData(:, 2); 
            y_ = contourData(:, 3);
            if strcmpi(obj.cieType, 'CIE1960')
                obj.x = 4 * x_ ./ (-2 * x_ + 12 * y_ + 3);
                obj.y = 6 * y_ ./ (-2 * x_ + 12 * y_ + 3);
            elseif strcmpi(obj.cieType, 'CIE1976')
                obj.x = 4 * x_ ./ (-2 * x_ + 12 * y_ + 3);
                obj.y = 9 * y_ ./ (-2 * x_ + 12 * y_ + 3);
            else
                obj.x = x_; 
                obj.y = y_;
            end
        end
        
        function chromDiagImage = getChromaDiagram(obj)
            filename = sprintf("%s_1000_divmax_.png", upper(obj.cieType));
            chromDiagImage = imread(filename);
        end
                 
        function refWhitePoint = getRefWhitePoint(obj)
            % ��ȡ�ο�������
            wp_rgb = uint8([255 255 255]);
            % wp_xyz = whitepoint('d65');
            refWhiteYxy = rgb2Yxy(wp_rgb, 'cieType', obj.cieType);
            refWhitePoint = refWhiteYxy(2:3);
        end
        
        function srgbVertices = getSRGBGamutVertices(obj)
            rgbVertices = uint8([255 0 0; 0 255 0; 0 0 255]);
            YxyVertices = rgb2Yxy(rgbVertices, 'cieType', obj.cieType);
            srgbVertices = YxyVertices(:,2:3);  
        end
    end
end


% % ���� ChromaMap ����
% chromDiag = ChromaDiagram('CIE1976');
% 
% % ��ȡ�ο�������
% whitePoint = chromDiag.getRefWhitePoint();
% 
% % ��ȡsRGB������ɫ��Ķ�������
% srgbVertices = chromDiag.getSRGBGamutVertices();
% 
% % ����ɫƷͼ
% chromDiag.initCoordinates();
% chromEdges = [chromDiag.x, chromDiag.y; chromDiag.x(1), chromDiag.y(1)];
% srgbGamut = [srgbVertices; srgbVertices(1,:)];
% axisRange = [0 max(chromEdges(:))+0.05];
% 
% figure
% plot(chromEdges(:,1), chromEdges(:,2), 'r','DisplayName','������')
% hold on
% plot(srgbGamut(:,1), srgbGamut(:,2), 'g','DisplayName','sRGBɫ��')
% scatter(whitePoint(1), whitePoint(2), 'b*','DisplayName','�ο��׵�')
% title(chromDiag.cieType); legend()
% axis equal; xlim(axisRange); ylim(axisRange);
% hold off