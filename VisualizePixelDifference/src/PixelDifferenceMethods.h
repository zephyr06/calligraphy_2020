// the configured settings and options for calculating pixel difference

#define PixelDifference_VERSION_MAJOR 1
#define PixelDifference_VERSION_MINOR 2

const int processImages(std::string originalImagePath, std::string resultImagePath, std::string resultFilename, int dimensions);
const Eigen::MatrixXd getImage(std::string path, int dimensions);
const Eigen::MatrixXd subtractMatrices(Eigen::MatrixXd original, Eigen::MatrixXd result);
void blendMatrices(Eigen::MatrixXd original, Eigen::MatrixXd result, std::string resultFilename);
const int calculatePixelDifference(Eigen::MatrixXd diffMatrix);
