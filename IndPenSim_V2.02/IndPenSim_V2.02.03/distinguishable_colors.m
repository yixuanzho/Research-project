function colors = distinguishable_colors(n_colors)
% DISTINGUISHABLE_COLORS - Create a set of distinct colors for plotting
%
% This is a simplified version of the distinguishable_colors function that
% returns a predefined set of distinct colors. The original function used
% perceptual distance calculations to generate maximally distinct colors.
%
% Input:
%   n_colors - Number of colors to generate
%
% Output:
%   colors - n_colors x 3 matrix of RGB values in the range [0,1]

% Define a set of base colors that are perceptually distinct
base_colors = [
    0.0000, 0.4470, 0.7410;  % Blue
    0.8500, 0.3250, 0.0980;  % Red
    0.9290, 0.6940, 0.1250;  % Yellow
    0.4940, 0.1840, 0.5560;  % Purple
    0.4660, 0.6740, 0.1880;  % Green
    0.3010, 0.7450, 0.9330;  % Light blue
    0.6350, 0.0780, 0.1840;  % Dark red
    0.0000, 0.0000, 1.0000;  % Bright blue
    0.0000, 0.5000, 0.0000;  % Dark green
    1.0000, 0.0000, 0.0000;  % Bright red
    0.7500, 0.7500, 0.0000;  % Olive
    0.7500, 0.0000, 0.7500;  % Magenta
    0.0000, 0.7500, 0.7500;  % Cyan
    0.7500, 0.7500, 0.7500;  % Gray
    0.1500, 0.1500, 0.1500;  % Dark gray
    0.0000, 1.0000, 0.0000;  % Bright green
    1.0000, 0.0000, 1.0000;  % Bright magenta
    0.0000, 1.0000, 1.0000;  % Bright cyan
    0.5000, 0.5000, 0.5000;  % Medium gray
    1.0000, 1.0000, 0.0000;  % Bright yellow
];

% Number of base colors
n_base = size(base_colors, 1);

% If we need more colors than we have in our base set
if n_colors > n_base
    % Repeat the base colors with slight variations
    colors = zeros(n_colors, 3);
    
    for i = 1:n_colors
        % Get base color (cycling through the base set)
        base_idx = mod(i-1, n_base) + 1;
        base_color = base_colors(base_idx, :);
        
        % Add small random variation for repeated colors
        if i > n_base
            % Add random variation (but ensure color stays in [0,1] range)
            variation = 0.1 * (rand(1,3) - 0.5);
            new_color = base_color + variation;
            new_color = min(max(new_color, 0), 1);  % Clip to [0,1]
            colors(i, :) = new_color;
        else
            colors(i, :) = base_color;
        end
    end
else
    % Just return the first n_colors from the base set
    colors = base_colors(1:n_colors, :);
end
end