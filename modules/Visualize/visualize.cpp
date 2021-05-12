// Copyright 2021 Oganyan Robert

#include <omp.h>
#include <vector>
#include <iostream>
#include <map>
#include <random>

#include <opencv2/core.hpp>
#include <opencv2/imgcodecs.hpp>
#include <opencv2/imgproc.hpp>


#include "../../modules/task_1/oganyan_r_mark_components/mark_components.h"
#include "../../modules/task_2/oganyan_r_mark_components_omp/mark_components_omp.h"


void convert_to_zeroone(cv::Mat *img) {
    for (int i = 0; i < (*img).rows; ++i) {
        for (int j = 0; j < (*img).cols; ++j) {
            if ((*img).at<uchar>(i, j) < 255) {
                (*img).at<uchar>(i, j) = 1;
            }
            if ((*img).at<uchar>(i, j) == 255) {
                (*img).at<uchar>(i, j) = 0;
            }
        }
    }
}


void convert_to_rdm_color(std::vector<cv::Vec3b> *img) {
    std::random_device rd;
    std::default_random_engine generator(rd());
    std::uniform_int_distribution<unsigned int> distribution(0, 255);
    std::map<uint, cv::Vec3b> colors;
    for (uint64_t i = 0; i < (*img).size(); ++i) {
        int cur_num = (*img)[i][0];
        if (cur_num > 1) {
            if (colors.find(cur_num) == colors.end()) {
                cv::Vec3b result;
                for (int j = 0; j < 3; ++j) {
                    result[j] = distribution(generator);
                }
                colors[cur_num] = result;
                (*img)[i] = result;
            } else {
                (*img)[i] = colors[cur_num];
            }
        } else {
            (*img)[i] = {255, 255, 255};
        }
    }
    return;
}

int main(int argc, char *argv[]) {
    std::cout << "Printing arvg info: " << "\n";
    for (int i = 0; i < argc; ++i) {
        std::cout << argv[i] << "\n";
    }

    if (argc < 3) {
        std::cout << "Wrong input" << "\n";
        return 0;
    }

    cv::Mat img = imread(argv[1], cv::IMREAD_GRAYSCALE);

    if (img.empty()) {
        std::cout << "Wrong data" << "\n";
        return 0;
    }
    convert_to_zeroone(&img);
    std::vector<int> array(img.datastart, img.dataend);
    std::pair<std::vector<int>, int> new_img = MarkComponentsPar(&array, img.rows, img.cols);
    std::vector<cv::Vec3b> result_image(new_img.first.begin(), new_img.first.end());
    convert_to_rdm_color(&result_image);
    cv::Mat final_img(img.rows, img.cols, CV_8UC3, result_image.data());
    cv::imwrite(argv[2], final_img, {cv::IMWRITE_JPEG_QUALITY});

    return 0;
}





