#include <iostream>
#include <iomanip>
#include <cmath>

#include "../include/progress_bar.h"

ProgressBar::ProgressBar(int total_length, double total_points): total_length{total_length}, total_points{total_points} {
    points = 0.0;
    this -> print_progress();
}

void ProgressBar::update_progress(double new_points) {
    points += new_points;
}

void ProgressBar::print_progress() {
    std::cout << "\r[" << std::string(std::floor(total_length * (points / total_points)), 'X') <<   //printing filled part
                     std::string(std::floor(total_length * (1 - points / total_points)), '-') <<    //printing empty part
              "] " << std::setprecision(3) << 100 * (points / total_points) << "%" << std::flush;
}

void ProgressBar::update_and_print_progress(double new_points) {
    this -> update_progress(new_points);
    this -> print_progress();
}