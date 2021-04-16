from PIL import Image
import cv2

# im = Image.open("resized_result_SI_scanned.png")
# im.save("resized_result_SI_scanned.png")

# set up filenames for reading files later
original_filename = '../Data/si_simulation/original_SI_from_SVG.png'
simple_brush_init_filename = '../Data/si_simulation/si_init_simpleBrush.png'
simple_brush_opt_filename = '../Data/si_simulation/si_opt_simpleBrush.png' 
dynamic_brush_init_filename = '../Data/si_simulation/simu_init_dynamicBrush.png'
dynamic_brush_opt_filename = '../Data/si_simulation/simu_opt_dynamicBrush.png'

# read images
original = cv2.imread(original_filename)
simple_brush_init = cv2.imread(simple_brush_init_filename)
simple_brush_opt = cv2.imread(simple_brush_opt_filename)
dynamic_brush_init = cv2.imread(dynamic_brush_init_filename)
dynamic_brush_opt = cv2.imread(dynamic_brush_opt_filename)

# resize images
dim = (128, 128)
simple_brush_init = cv2.resize(simple_brush_init, dim, interpolation=cv2.INTER_NEAREST)
simple_brush_opt = cv2.resize(simple_brush_opt, dim, interpolation=cv2.INTER_NEAREST)
dynamic_brush_init = cv2.resize(dynamic_brush_init, dim, interpolation=cv2.INTER_NEAREST)
dynamic_brush_opt = cv2.resize(dynamic_brush_opt, dim, interpolation=cv2.INTER_NEAREST)

# save the resized images for calculating pixel difference later
cv2.imwrite('./resized_SI/resized_SI_original.png',           original)
cv2.imwrite('./resized_SI/resized_SI_simple_brush_init.png',  simple_brush_init)
cv2.imwrite('./resized_SI/resized_SI_simple_brush_opt.png',   simple_brush_opt)
cv2.imwrite('./resized_SI/resized_SI_dynamic_brush_init.png', dynamic_brush_init)
cv2.imwrite('./resized_SI/resized_SI_dynamic_brush_opt.png',  dynamic_brush_opt)

print('original: ', original.shape)
print('simple_brush_init: ', simple_brush_init.shape)
print('simple_brush_opt: ', simple_brush_opt.shape)
print('dynamic_brush_init: ', dynamic_brush_init.shape)
print('dynamic_brush_opt: ', dynamic_brush_opt.shape)

# Using cv2.imwrite() method 
# Saving the image
# cv2.imwrite('resized_' + filename, resized) 

# blend images with original
alpha = 0.5
simple_brush_init = cv2.addWeighted(original, alpha, simple_brush_init, alpha, 0.0)
cv2.imwrite('./blended_SI/blended_SI_simple_brush_init.png', simple_brush_init)

simple_brush_opt = cv2.addWeighted(original, alpha, simple_brush_opt, alpha, 0.0)
cv2.imwrite('./blended_SI/blended_SI_simple_brush_opt.png', simple_brush_opt)

dynamic_brush_init = cv2.addWeighted(original, alpha, dynamic_brush_init, alpha, 0.0)
cv2.imwrite('./blended_SI/blended_SI_dynamic_brush_init.png', dynamic_brush_init)

dynamic_brush_opt = cv2.addWeighted(original, alpha, dynamic_brush_opt, alpha, 0.0)
cv2.imwrite('./blended_SI/blended_SI_dynamic_brush_opt.png', dynamic_brush_opt)
