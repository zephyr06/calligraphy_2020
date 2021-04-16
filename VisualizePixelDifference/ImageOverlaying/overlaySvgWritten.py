from PIL import Image
import cv2

one_name = "Wo-initial.png"
two_name = "Wo-opt-dynamic.png"
three_name = "Wo-opt-dynamic1.png"
four_name = "Wo-opt-dynamic2.png"

# set up filenames for reading files later
original_filename = '../Data/svgtopng/wo_svg2png.png'
one_filename = '../Data/written_results_scanned_using_scanner/' + one_name
two_filename = '../Data/written_results_scanned_using_scanner/' + two_name
three_filename = '../Data/written_results_scanned_using_scanner/' + three_name
four_filename = '../Data/written_results_scanned_using_scanner/' + four_name

# read images
original = cv2.imread(original_filename)
one = cv2.imread(one_filename)
two = cv2.imread(two_filename)
three = cv2.imread(three_filename)
four = cv2.imread(four_filename)

# resize images
dim = (640, 640)
original = cv2.resize(original, dim, interpolation=cv2.INTER_NEAREST)
one = cv2.resize(one, dim, interpolation=cv2.INTER_NEAREST)
two = cv2.resize(two, dim, interpolation=cv2.INTER_NEAREST)
three = cv2.resize(three, dim, interpolation=cv2.INTER_NEAREST)
four = cv2.resize(four, dim, interpolation=cv2.INTER_NEAREST)

# set the written image to RED for visibility
Conv_hsv_Gray = cv2.cvtColor(one, cv2.COLOR_BGR2GRAY)
ret, mask = cv2.threshold(Conv_hsv_Gray, 0, 255,cv2.THRESH_BINARY_INV |cv2.THRESH_OTSU)
one[mask == 255] = [0, 0, 255]
Conv_hsv_Gray = cv2.cvtColor(two, cv2.COLOR_BGR2GRAY)
ret, mask = cv2.threshold(Conv_hsv_Gray, 0, 255,cv2.THRESH_BINARY_INV |cv2.THRESH_OTSU)
two[mask == 255] = [0, 0, 255]
Conv_hsv_Gray = cv2.cvtColor(three, cv2.COLOR_BGR2GRAY)
ret, mask = cv2.threshold(Conv_hsv_Gray, 0, 255,cv2.THRESH_BINARY_INV |cv2.THRESH_OTSU)
three[mask == 255] = [0, 0, 255]
Conv_hsv_Gray = cv2.cvtColor(four, cv2.COLOR_BGR2GRAY)
ret, mask = cv2.threshold(Conv_hsv_Gray, 0, 255,cv2.THRESH_BINARY_INV |cv2.THRESH_OTSU)
four[mask == 255] = [0, 0, 255]

# save the resized images for calculating pixel difference later
# cv2.imwrite('./blended_WO_scanner/blended_WO_scanner_original.png',           original)
# cv2.imwrite('./blended_WO_scanner/blended_WO_scanner_simple_brush_init.png',  simple_brush_init)
# cv2.imwrite('./blended_WO_scanner/blended_WO_scanner_simple_brush_opt.png',   simple_brush_opt)
# cv2.imwrite('./blended_WO_scanner/blended_WO_scanner_dynamic_brush_init.png', dynamic_brush_init)
# cv2.imwrite('./blended_WO_scanner/blended_WO_scanner_dynamic_brush_opt.png',  dynamic_brush_opt)

print('original: ', original.shape)
print('one: ', one.shape)
print('two: ', two.shape)
print('three: ', three.shape)
print('four: ', four.shape)

# Using cv2.imwrite() method 
# Saving the image
# cv2.imwrite('resized_' + filename, resized) 

# blend images with original
alpha = 0.4
beta = 1 - alpha
one = cv2.addWeighted(original, alpha, one, beta, 0.0)
cv2.imwrite('./blended_WO_scanner/' + one_name, one)

two = cv2.addWeighted(original, alpha, two, beta, 0.0)
cv2.imwrite('./blended_WO_scanner/' + two_name, two)

three = cv2.addWeighted(original, alpha, three, beta, 0.0)
cv2.imwrite('./blended_WO_scanner/' + three_name, three)

four = cv2.addWeighted(original, alpha, four, beta, 0.0)
cv2.imwrite('./blended_WO_scanner/' + four_name, four)