clear;
global fil datasetStart datasetEnd;
fil='20140222_01_01_03_250lm';
datasetStart = 15000; %Elimino il tratto a vx = 0.
datasetEnd   = 70000;
    
load(sprintf('%s.mat',fil));
[~,~]=mkdir(fil);

saveOne(insData.ayCG.time, 't');
saveOne(insData.ayCG.value, 'ayCG');
saveOne(insData.ayCGFilt.value, 'ayCGFilt');
saveOne(insData.axCG.value, 'axCG');
saveOne(insData.axCGFilt.value, 'axCGFilt');
saveOne(insData.yawRate.value* pi / 180, 'yawRate');
saveOne(insData.yawRateFilt.value, 'yawRateFilt');

saveOne(insData.yawAngAcc.value* pi / 180, 'yawAngAcc');
saveOne(insData.yawAngAccFilt.value* pi / 180, 'yawAngAccFilt');

saveOneDecim(tireData.roadWheelAngle.value* pi / 180, 'delta');
saveOneDecim(tireData.roadWheelAngleFL.value* pi / 180, 'deltaFL');
saveOneDecim(tireData.roadWheelAngleFR.value* pi / 180, 'deltaFR');

saveOne(insData.vxCG.value, 'vx');
saveOne(insData.vyCG.value, 'vy');
saveOne(insData.sideSlip.value*pi/180, 'beta_true');

function []=saveOne(d,filename)
    global fil datasetStart datasetEnd;
    d=d(datasetStart:datasetEnd); 
    save(sprintf('%s/%s.txt',fil,filename),'d','-ascii');
end

function []=saveOneDecim(d,filename)
    saveOne(d(1:10:end),filename);
end
