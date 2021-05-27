//// Copyright 2021 Oganyan Robert

#include <vector>
#include <iostream>
#include <map>
#include <random>
#include <iomanip>

#include <opencv2/core.hpp>
#include <opencv2/imgcodecs.hpp>
#include <opencv2/imgproc.hpp>
#include <omp.h>
//#include <tbb/parallel_for.h>
//#include <tbb/blocked_range2d.h>
//#include <tbb/task_scheduler_init.h>
#include <thread>
#include <queue>



#include "../../modules/task_1/oganyan_r_mark_components/mark_components.h"
#include "../../modules/task_2/oganyan_r_mark_components_omp/mark_components_tbb.h"
#include "../../modules/task_3/oganyan_r_mark_components_tbb/mark_components_tbb.h"
#include "../../modules//task_4/oganyan_r_mark_components_std/mark_components_std.h"

//static const std::vector<std::pair<int, int>> directions{
//        {-1, 0},
//        {0,  -1},
//        {0,  1},
//        {1,  0},
//};
//
//void bfs(std::vector<uint16_t>* img, std::pair<uint16_t, uint16_t> start,
//         uint16_t* number, uint16_t width, uint16_t height) {
//    if ((*img)[start.first * width + start.second] != 1) {
//        return;
//    }
//    std::queue<std::pair<uint16_t, uint16_t>> q;
//    q.push({start});
//    (*img)[start.first * width + start.second] = ++(*number);
//    while (!q.empty()) {
//        auto cur{q.front()};
//        q.pop();
//        for (auto &neighbor : directions) {
//            if (cur.first + neighbor.first >= height || cur.first + neighbor.first < 0
//                || cur.second + neighbor.second >= width || cur.second + neighbor.second < 0) {
//                continue;
//            }
//            if ((*img)[(cur.first + neighbor.first) * width + cur.second + neighbor.second] == 1) {
//                q.push({(cur.first + neighbor.first), cur.second + neighbor.second});
//                (*img)[(cur.first + neighbor.first) * width + cur.second + neighbor.second] = (*number);
//            }
//        }
//    }
//    return;
//}
//
//
//std::pair<std::vector<uint16_t>, uint16_t> MarkComponents(std::vector<uint16_t> img, uint16_t height, uint16_t width) {
//    if (img.size() == 0) {
//        throw std::invalid_argument("Size of the image cant be negative");
//    }
//    if (static_cast<int>(img.size()) != width * height) {
//        throw std::invalid_argument("Size of the image is not correct");
//    }
//    uint16_t count_comp{1};
//    for (uint16_t i{0}; i < height; ++i) {
//        for (uint16_t j{0}; j < width; ++j) {
//            bfs(&img, {i, j}, &count_comp, width, height);
//        }
//    }
//    return {img, count_comp - 1};
//}
//
//template <class T>
//class Disjoint_Set_Union {
//private:
//    std::vector<T> rank;
//    std::vector<T> parent;
//
//public:
//    explicit Disjoint_Set_Union(int size) : rank(size), parent(size) {
//    }
//
//    void make_set(int vertex) {
//        parent[vertex] = vertex;
//        rank[vertex] = 0;
//    }
//
//
//    void Init() {
//        for (std::size_t vertex = 0; vertex < rank.size(); ++vertex) {
//            make_set(vertex);
//        }
//    }
//
//
//    void Init(int size) {
//        rank.resize(size);
//        parent.resize(size);
//        for (int vertex = 0; vertex < size; ++vertex) {
//            make_set(vertex);
//        }
//    }
//
//    int find_set(int vertex) {
//        //  Bug here (parent[vertex] = [parent[paren[vertex]])
//        if (vertex == parent[vertex]) {
//            return vertex;
//        }
//
//        return parent[vertex] = find_set(parent[vertex]);
//    }
//
//    int find_set(int vertex, int last) {
//        if (vertex == parent[vertex]) {
//            return vertex;
//        }
//
//        if (last == parent[vertex]) {
//            return std::min(last, vertex);
//        }
//
//        return parent[vertex] = find_set(parent[vertex], vertex);
//    }
//
//    void union_sets(int fi_union, int se_union) {
//        fi_union = find_set(fi_union, fi_union);
//        se_union = find_set(se_union, se_union);
//        if (fi_union != se_union) {
//            if (rank[fi_union] < rank[se_union])
//                std::swap(fi_union, se_union);
//            parent[se_union] = fi_union;
//            if (rank[fi_union] == rank[se_union])
//                ++rank[fi_union];
//        }
//    }
//
//    const std::vector<T> get_rank() {
//        return this->rank;
//    }
//
//    const std::vector<T> get_parent() {
//        return this->parent;
//    }
//};
//
//std::pair<std::vector<int>, int> MarkComponentsPar(std::vector<int> *img, int height, int width) {
//    if ((*img).size() == 0) {
//        throw std::invalid_argument("Size of the image cant be negative");
//    }
//    if (static_cast<int>((*img).size()) != width * height) {
//        throw std::invalid_argument("Size of the image is not correct");
//    }
//    auto img_new = *img;
//    int count_comp{0};
//    Disjoint_Set_Union<int> dsu(height * width);
//    dsu.Init();
//    omp_set_num_threads(8);
//
//#pragma omp parallel default(none) shared(img_new, width, height, count_comp, dsu, directions)
//    {
//#pragma omp for schedule(static)
//        for (int i = 0; i < height; ++i) {
//            for (int j = 0; j < width; ++j) {
//                if (img_new[i * width + j]) {
//                    for (auto &neighbor : directions) {
//                        if (i + neighbor.first >= height || i + neighbor.first < 0
//                            || j + neighbor.second >= width || j + neighbor.second < 0) {
//                            continue;
//                        }
//                        int cur = (i + neighbor.first) * width + j + neighbor.second;
//                        if (img_new[cur] == 0) {
//                            continue;
//                        }
//                        dsu.union_sets(i * width + j, cur);
//                    }
//                }
//            }
//        }
//    }
//
//
//    for (int i = 0; i < height; ++i) {
//        for (int j = 0; j < width; ++j) {
//            int cur = i * width + j;
//            if (img_new[cur]) {
//                img_new[cur] = dsu.find_set(cur, cur) + 1;
//            }
//            if (img_new[cur] == cur + 1) {
//                count_comp++;
//            }
//        }
//    }
//
//    return {img_new, count_comp};
//}

//
//std::pair<std::vector<int>, int> MarkComponentsParTbb(std::vector<int> *img, int height, int width, int proc_num) {
//    if ((*img).size() == 0) {
//        throw std::invalid_argument("Size of the image cant be negative");
//    }
//    if (static_cast<int>((*img).size()) != width * height) {
//        throw std::invalid_argument("Size of the image is not correct");
//    }
//    auto img_new = *img;
//    int count_comp{0};
//    Disjoint_Set_Union<int> dsu(height * width);
//    dsu.Init();
//
//    tbb::task_scheduler_init init(proc_num);
//    tbb::parallel_for(tbb::blocked_range<int>(0, height, height/proc_num + 1),
//                      [&img_new, &dsu, width, height](tbb::blocked_range<int> block) {
//                          for (auto i = block.begin(); i != block.end(); ++i) {
//                              for (int j = 0; j < width; ++j) {
//                                  if (img_new[i * width + j]) {
//                                      for (auto &neighbor : directions) {
//                                          if (i + neighbor.first >= height || i + neighbor.first < 0
//                                              || j + neighbor.second >= width || j + neighbor.second < 0) {
//                                              continue;
//                                          }
//                                          int cur = (i + neighbor.first) * width + j + neighbor.second;
//                                          if (img_new[cur] == 0) {
//                                              continue;
//                                          }
//                                          dsu.union_sets(i * width + j, cur);
//                                      }
//                                  }
//                              }
//                          }
//                      });
//
//    //init.terminate();
//    for (int i = 0; i < height; ++i) {
//        for (int j = 0; j < width; ++j) {
//            int cur = i * width + j;
//            if (img_new[cur]) {
//                img_new[cur] = dsu.find_set(cur, cur) + 1;
//            }
//            if (img_new[cur] == cur + 1) {
//                count_comp++;
//            }
//        }
//    }
//
//    return {img_new, count_comp};
//}

//
//std::pair<std::vector<int>, int> MarkComponentsParStd(std::vector<int> *img, int height, int width, int num_proc) {
//    if ((*img).size() == 0) {
//        throw std::invalid_argument("Size of the image cant be negative");
//    }
//    if (static_cast<int>((*img).size()) != width * height) {
//        throw std::invalid_argument("Size of the image is not correct");
//    }
//    auto img_new = *img;
//    int count_comp{0};
//    Disjoint_Set_Union<int> dsu(height * width);
//    dsu.Init();
//
//    std::vector<std::thread> threads;
//    threads.reserve(num_proc);
//
//
//    int div = height / num_proc;
//    int mod = height % num_proc;
//
//    int last = 0;
//
//
//    for (int proc = 0; proc < num_proc; ++proc) {
//        threads.emplace_back([&img_new, width, height, &dsu, proc, div, mod, last]() {
//            int balance = (proc < mod) ? 1 : 0;
//            for (int i = last; i < last + div + balance; ++i) {
//                for (int j = 0; j < width; ++j) {
//                    if (img_new[i * width + j]) {
//                        for (auto &neighbor : directions) {
//                            if (i + neighbor.first >= height || i + neighbor.first < 0
//                                || j + neighbor.second >= width || j + neighbor.second < 0) {
//                                continue;
//                            }
//                            int cur = (i + neighbor.first) * width + j + neighbor.second;
//                            if (img_new[cur] == 0) {
//                                continue;
//                            }
//                            dsu.union_sets(i * width + j, cur);
//                        }
//                    }
//                }
//            }
//        });
//        last += div + ((proc < mod) ? 1 : 0);
//    }
//
//    for (auto &thread : threads) {
//        thread.join();
//    }
//
//
//    for (int i = 0; i < height; ++i) {
//        for (int j = 0; j < width; ++j) {
//            int cur = i * width + j;
//            if (img_new[cur]) {
//                img_new[cur] = dsu.find_set(cur, cur) + 1;
//            }
//            if (img_new[cur] == cur + 1) {
//                count_comp++;
//            }
//        }
//    }
//
//    return {img_new, count_comp};
//}

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
    // Sequential

    std::vector<uint16_t> array_seq(img.datastart, img.dataend);
    auto start = std::chrono::high_resolution_clock::now();
    std::pair<std::vector<uint16_t>, uint16_t> new_img_seq = MarkComponents(array_seq, img.rows, img.cols);
    auto end = std::chrono::high_resolution_clock::now();
    std::chrono::duration<double> diff = end - start;
    std::cout << "Sequential: " << std::setw(9) << diff.count() << "\n";
    std::vector<cv::Vec3b> result_image_seq(new_img_seq.first.begin(), new_img_seq.first.end());
    convert_to_rdm_color(&result_image_seq);
    cv::Mat final_img_seq(img.rows, img.cols, CV_8UC3, result_image_seq.data());
    cv::imwrite("/home/ogrob/pp_2021_spring_engineers/modules/Visualize/results/seq.png", final_img_seq,
                {cv::IMWRITE_JPEG_QUALITY});

    // OMP

    std::vector<int> array_omp(img.datastart, img.dataend);
    start = std::chrono::high_resolution_clock::now();
    std::pair<std::vector<int>, int> new_img_omp = MarkComponentsPar(&array_omp, img.rows, img.cols);
    end = std::chrono::high_resolution_clock::now();
    diff = end - start;
    std::cout << "OMP: " << std::setw(9) << diff.count() << "\n";
    std::vector<cv::Vec3b> result_image_omp(new_img_omp.first.begin(), new_img_omp.first.end());
    convert_to_rdm_color(&result_image_omp);
    cv::Mat final_img_omp(img.rows, img.cols, CV_8UC3, result_image_omp.data());
    cv::imwrite("/home/ogrob/pp_2021_spring_engineers/modules/Visualize/results/omp.png", final_img_omp,
                {cv::IMWRITE_JPEG_QUALITY});

    // TBB

    std::vector<int> array_tbb(img.datastart, img.dataend);
    start = std::chrono::high_resolution_clock::now();
    std::pair<std::vector<int>, int> new_img_tbb = MarkComponentsParTbb(&array_tbb, img.rows, img.cols, 8);
    end = std::chrono::high_resolution_clock::now();
    diff = end - start;
    std::cout << "OMP: " << std::setw(9) << diff.count() << "\n";
    std::vector<cv::Vec3b> result_image_tbb(new_img_tbb.first.begin(), new_img_tbb.first.end());
    convert_to_rdm_color(&result_image_tbb);
    cv::Mat final_img_tbb(img.rows, img.cols, CV_8UC3, result_image_tbb.data());
    cv::imwrite("/home/ogr/pp_2021_spring_engineers/modules/Visualize/results/tbb.png", final_img_tbb,
                {cv::IMWRITE_JPEG_QUALITY});

    // STD

    std::vector<int> array_std(img.datastart, img.dataend);
    start = std::chrono::high_resolution_clock::now();
    std::pair<std::vector<int>, int> new_img_std = MarkComponentsParStd(&array_omp, img.rows, img.cols, 4);
    end = std::chrono::high_resolution_clock::now();
    diff = end - start;
    std::cout << "OMP: " << std::setw(9) << diff.count() << "\n";
    std::vector<cv::Vec3b> result_image_std(new_img_std.first.begin(), new_img_std.first.end());
    convert_to_rdm_color(&result_image_std);
    cv::Mat final_img_std(img.rows, img.cols, CV_8UC3, result_image_std.data());
    cv::imwrite("/home/ogrob/pp_2021_spring_engineers/modules/Visualize/results/std.png", final_img_std,
                {cv::IMWRITE_JPEG_QUALITY});


    return 0;
}
