//
// Created by vladu on 03.08.2021.
//

#ifndef MPI_PROJECT_FUNCTIONS_H
#define MPI_PROJECT_FUNCTIONS_H

#include <boost/multi_array.hpp>
#include <iostream>
#include <fstream>
#include <vector>
#include "cmath"
#include <boost/numeric/ublas/matrix.hpp>
//#include "ImageMagick-7/Magick++.h"

typedef boost::numeric::ublas::matrix<double> array_2D;
//typedef boost::multi_array<double,2> array_2D;
//typedef array_2D::index index;

class bad_index:std::exception{
};

int save_state_file(array_2D temp_grid, const std::string& path);

//int draw_image(twoD_array<double> temp_grid, const std::string& path, double max_temp);

//Magick::Color get_color_from_temp(double temperature, double max);

class Mpi_config{
public:
    double conductivity, height, width, delta_x, delta_y, delta_t;
    long long save_rate, cycle_duration, max_temp;
    int processes;
    long grid_height, grid_width;
    array_2D starting_condition;
    Mpi_config(const std::string& config="../config.dat", const std::string& starting_cond="../base_cond.txt");
};
#endif //MPI_PROJECT_FUNCTIONS_H
