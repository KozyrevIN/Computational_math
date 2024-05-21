#include <string>

class ProgressBar {
private:
    const int total_length;
    double total_points;
    double points;
    std::string text;
    bool updatable;

public:
    ProgressBar(int total_length, double total_points);
    ProgressBar(double total_points, std::string text);

    void update_progress(double points);
    void print_progress();
    void update_and_print_progress(double points);
    void set_100();
};
