
class ProgressBar {
private:
    const int total_length;
    double total_points;
    double points;
    bool updatable;

public:
    ProgressBar(int, double);

    void update_progress(double);
    void print_progress();
    void update_and_print_progress(double);
    void set_100();
};
