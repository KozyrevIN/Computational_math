#include <iostream>
#include <iomanip>
#include <cmath>

#ifndef progress_bar
#define progress_bar
    #include "../include/progress_bar.h"
#endif

ProgressBar::ProgressBar(int total_length, double total_points): total_length(total_length), total_points(total_points) {
    updatable = true;
    points = 0.0;
    this -> print_progress();
}

ProgressBar::ProgressBar(double total_points, std::string text): total_length(64), total_points(total_points), text(text) {
    updatable = true;
    points = 0.0;
    this -> print_progress();
}

void ProgressBar::update_progress(double new_points) {
    points += new_points;
}

void ProgressBar::print_progress() {
    if (updatable) {
        std::cout << "\r" << text << " [" << std::string(round(total_length * (points / total_points)), '#') <<   //printing filled part
                              std::string(total_length - round(total_length * (points / total_points)), '-') <<    //printing empty part
                       "] " << round(100 * (points / total_points)) << "% " <<std::flush;
        if (round(100 * (points / total_points)) == 100) {
            std::cout << '\n';
            updatable = false;
        }
    }
}

void ProgressBar::update_and_print_progress(double new_points) {
    this -> update_progress(new_points);
    this -> print_progress();
}

void ProgressBar::set_100() {
    points = total_points;
    this -> print_progress();
}