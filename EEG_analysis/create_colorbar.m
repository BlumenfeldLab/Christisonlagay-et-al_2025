
function mycmap = create_colorbar(outerbounds, innerbounds)

% define different blocks
b1 = innerbounds(1) - outerbounds(1);
b2 = innerbounds(2) - innerbounds(1);
b3 = outerbounds(2) - innerbounds(2);
blocks = [b1, b2, b3];
blocks = blocks / min(blocks);

% find nearest integer ratio between a minumum of 5 and maximum of 30 bins (arbitrary)
span = (20 : 30)';
block_span = span * blocks;
res = sum(rem(block_span, 1), 2);
[~, min_res] = min(res);
block = round(block_span(min_res, :));

% define colors
darkblue = [0, 0, .5];
blue = [0, 0, 1];
turquoise = [0.4, 0.9, 0.9];
gray = [.9059, .9059, .9059];
yellow = [1, 1, 0];
red = [1, 0, 0];
darkred = [.5, 0, 0];

% construct dark blue to blue section
block1_1 = floor(block(1) / 3);
block1_2 = block1_1;
block1_3 = block(1) - block1_1 - block1_2;
b1_1 = [];
for i = 1 : 3
    b1_1(:, i) = linspace(darkblue(i), blue(i), block1_1)';
end
b1_2 = [];
for i = 1 : 3
    b1_2(:, i) = linspace(blue(i), turquoise(i), block1_2)';
end
b1_3 = [];
for i = 1 : 3
    b1_3(:, i) = linspace(turquoise(i), gray(i), block1_3)';
end
b1 = [b1_1; b1_2; b1_3];

% construct gray section
b2 = repmat(gray, block(2), 1);

% construct yellow to red to dark red section
block3_1 = floor(block(3) / 3);
block3_2 = block3_1;
block3_3 = block(3) - block3_1 - block3_2;
b3_1 = [];
for i = 1 : 3
    b3_1(:, i) = linspace(gray(i), yellow(i), block3_1)';
end
b3_2 = [];
for i = 1 : 3
    b3_2(:, i) = linspace(yellow(i), red(i), block3_2)';
end
b3_3 = [];
for i = 1 : 3
    b3_3(:, i) = linspace(red(i), darkred(i), block3_3)';
end
b3 = [b3_1; b3_2; b3_3];

% compile cmap
mycmap = [b1(1 : end - 1, :); b2(1 : end - 1, :); b3];
