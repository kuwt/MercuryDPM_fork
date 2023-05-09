./CalibrationHeap -psd 0.01 0 0.01 1 -species LinearViscoelasticFrictionReversibleAdhesiveSpecies -density 1000 -collisionTime 0.01 -restitutionCoefficient 0.5 -slidingFriction 0.5 -param X -output -constantRestitution -bondNumber 1 -rollingFriction 0.1
./CalibrationDrum -psd 0.01 0 0.01 1 -species LinearViscoelasticFrictionReversibleAdhesiveSpecies -density 1000 -collisionTime 0.01 -restitutionCoefficient 0.5 -slidingFriction 0.5 -param X -output -constantRestitution -bondNumber 1 -rollingFriction 0.1
./CalibrationShearCell -psd 0.01 0 0.01 1 -species LinearViscoelasticFrictionReversibleAdhesiveSpecies -density 1000 -collisionTime 0.01 -restitutionCoefficient 0.5 -slidingFriction 0.5 -param X -output -constantRestitution -bondNumber 1 -rollingFriction 0.1 -normalStress 1000 800 600 400

./CalibrationHeap -psdVolumetricCumulativeDiameter 0.01 0 0.01 1 -species LinearViscoelasticFrictionReversibleAdhesiveSpecies -density 1000 -collisionTime 0.01 -restitutionCoefficient 0.5 -slidingFriction 0.5 -param X -output -constantRestitution -bondNumber 1 -rollingFriction 0.1

species LinearViscoelasticFrictionReversibleAdhesiveSpecies
density 1500
psdVolumetricCumulativeDiameter
0.0005 0
0.0005 1
cutOffMin 0
cutOffMax 1
scaleUp 1
squeeze 1
collisionTime 0.000067997
torsionFriction 0
# Parameters to be calibrated (names and ranges)
restitutionCoefficient    0.5 1
slidingFriction           0   1
rollingFriction           0   1
cohesion                  0   0.7
# Experimental data to be tested (name, values and weights)
CalibrationShearCell     759 .4 609 .4 423 .4 269 .4
CalibrationHeap 29.725 1
CalibrationDrum          40.859 1
# Additional parameters for executables
-normalStress 1000 800 600 400 200