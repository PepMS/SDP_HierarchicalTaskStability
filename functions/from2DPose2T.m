function T = from2DPose2T(pos, ori)

T = [...
    cos(ori) -sin(ori) 0 pos(1);...
    sin(ori)  cos(ori) 0 pos(2);...
    0         0        1 0     ;...
    0         0        0 1]; 