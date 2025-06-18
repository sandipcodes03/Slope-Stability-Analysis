% Interactive Slope Stability Analysis using Simplified Bishop's Method
clear; clc;

% Prompting user for inputs
H = input('Enter the height of the slope (in meters): ');
beta = input('Enter the slope angle (in degrees): ');
r = input('Enter the radius of the slip surface (in meters): ');
nSlices = input('Enter the number of slices for analysis: ');

% Number of soil layers
nLayers = input('Enter the number of soil layers: ');

% Initialize layers matrix
layers = zeros(nLayers, 3);  % [Cohesion (c), Friction Angle (phi), Unit Weight (gamma)]
for i = 1:nLayers
    fprintf('For soil layer %d:\n', i);
    layers(i, 1) = input('  Enter cohesion (c) in kPa: ');
    layers(i, 2) = input('  Enter friction angle (phi) in degrees: ');
    layers(i, 3) = input('  Enter unit weight (gamma) in kN/m^3: ');
end

% Prompt for pore water pressure function details
porePressureOption = input('Do you want to consider pore water pressure? (1 for Yes, 0 for No): ');

if porePressureOption == 1
    porePressureType = input('  Choose pore water pressure distribution (1 for linear, 2 for constant): ');
    if porePressureType == 1
        maxPorePressure = input('  Enter the maximum pore water pressure at the base (in kPa): ');
        porePressureFunc = @(y) maxPorePressure * (1 - y / H);  % Linear distribution
    elseif porePressureType == 2
        constantPorePressure = input('  Enter constant pore water pressure (in kPa): ');
        porePressureFunc = @(y) constantPorePressure;  % Constant value
    end
else
    porePressureFunc = @(y) 0;  % No pore water pressure
end

% Coordinates for the center of the circular slip surface
xCenter = 0;  % Horizontal distance of the center of circle
yCenter = H - r;  % Vertical distance of the center of circle

% Slice calculation (for a non-homogeneous slope)
sliceWidth = 2 * r * sind(beta) / nSlices;
thetaStart = atand(H / (2 * r));
thetaEnd = beta;
theta = linspace(thetaStart, thetaEnd, nSlices);

% Initialize summation variables
sumW = 0;
sumS = 0;
sumM = 0;

% Loop through each slice for Factor of Safety calculation
for i = 1:nSlices
    % Calculate the coordinates for the base of each slice
    xi = r * sind(theta(i));
    yi = r * cosd(theta(i));
    
    % Determine which layer the slice is in
    if yi > H/2
        layer = 1;  % Top layer
    else
        layer = min(nLayers, 2);  % Adjust based on user input
    end
    
    % Get soil properties from the current layer
    c = layers(layer, 1);
    phi = layers(layer, 2);
    gamma = layers(layer, 3);
    
    % Weight of the slice
    Wi = gamma * sliceWidth * (yCenter - yi);
    
    % Pore water pressure at the base of the slice
    u = porePressureFunc(yi);  % Pore water pressure function
    
    % Effective normal force on the slice
    N = Wi * cosd(theta(i)) - u * sliceWidth;
    
    % Shear strength of the soil at the base of the slice (Mohr-Coulomb)
    S = c + N * tand(phi);
    
    % Resisting force
    R = S * sliceWidth;
    
    % Summing for Factor of Safety calculation
    sumW = sumW + Wi * sind(theta(i));  % Driving force
    sumS = sumS + S;  % Shear strength
    sumM = sumM + R;  % Resisting moment
end

% Factor of Safety calculation (Simplified Bishop's Method)
FoS = sumS / sumW;

% Display the result
fprintf('The Factor of Safety (FoS) is: %.3f\n', FoS);
