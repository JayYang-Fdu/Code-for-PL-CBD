function matr = givensB(h1,h2,arg)
%UNTITLED2 此处显示有关此函数的摘要
%   此处显示详细说明
% '-u'表示 第二个元素为0
% '-d'表示 第一个元素为0
angle1 = angle(h1);
angle2 = angle(h2);
deltaAngle = angle2-angle1;
if arg == "RowGivens-u"

    c = sqrt(abs(h1)^2/(abs(h2)^2+abs(h1)^2));
    s = sqrt(abs(h2)^2/(abs(h2)^2+abs(h1)^2))*exp(-1j*deltaAngle);
    matr = [c, -conj(s);s, c];
    if (abs(angle1)<0.001 ||abs(angle1-pi)<0.001) && (abs(angle2)<0.001 ||abs(angle2-pi)<0.001)
        matr = real([c, -conj(s);s, c]);
    end

elseif arg == "ColGivens-u"

    c = sqrt(abs(h1)^2/(abs(h2)^2+abs(h1)^2));
    s = sqrt(abs(h2)^2/(abs(h2)^2+abs(h1)^2))*exp(-1j*deltaAngle);
    matr = [c, s;-conj(s), c];
    if (abs(angle1)<0.001 ||abs(angle1-pi)<0.001) && (abs(angle2)<0.001 ||abs(angle2-pi)<0.001)
        matr = real([c, s;-conj(s), c]);
    end

elseif arg == "RowGivens-d"

    c = sqrt(abs(h2)^2/(abs(h2)^2+abs(h1)^2));
    s = sqrt(abs(h1)^2/(abs(h2)^2+abs(h1)^2))*exp(-1j*deltaAngle);
    matr = [c, conj(s);-s, c];
    if (abs(angle1)<0.001 ||abs(angle1-pi)<0.001) && (abs(angle2)<0.001 ||abs(angle2-pi)<0.001)
        matr = real([c, conj(s);-s, c]);
    end

elseif arg == "ColGivens-d"

    c = sqrt(abs(h2)^2/(abs(h2)^2+abs(h1)^2));
    s = sqrt(abs(h1)^2/(abs(h2)^2+abs(h1)^2))*exp(-1j*deltaAngle);
    matr = [c, -s;conj(s), c];
    if (abs(angle1)<0.001 ||abs(angle1-pi)<0.001) && (abs(angle2)<0.001 ||abs(angle2-pi)<0.001)
        matr = real([c, -s;conj(s), c]);
    end
end


end