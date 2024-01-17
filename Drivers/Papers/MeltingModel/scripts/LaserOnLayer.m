xMax = 0.00055;
meanRadius = 3.58591e-05;
laserPower = 2.3; 
spotRadius = 250e-6;
porosity = 0.55;
extinctionCoefficient = 3.0*(1-porosity)/(4.0*porosity*meanRadius);
peakIntensity = 2.0*laserPower/pi/spotRadius^2;
surfaceHeight = 8.25e-05;
hatchWidth = 0.7*spotRadius;
hatchLength = xMax-2.0*hatchWidth;
laserSpeed = 0.1;
halfPeriod = (hatchWidth+hatchLength)/laserSpeed;
maxTime = floor(hatchLength/hatchWidth)*halfPeriod+corner;

time=linspace(0,maxTime,1000);
line = floor(time/halfPeriod) + 1;
progress = mod(time,halfPeriod);
corner = hatchLength/(hatchLength+hatchWidth)*halfPeriod; 
centerX = hatchWidth + min(progress,corner)*laserSpeed;
centerXEven = hatchWidth + (corner-min(progress,corner))*laserSpeed;
id = mod(line,2)==0;
centerX(id) = centerXEven(id);
centerY = line*hatchWidth + max(progress-corner,0.0) * laserSpeed;
plot(centerX, centerY, '.')
xlim([0 xMax])
ylim([0 xMax])
