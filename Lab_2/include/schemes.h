#ifndef solver
#define solver
    #include "solver.h"
#endif

class LeftAngleSolver : public Solver {
    public:
    LeftAngleSolver(int, int, int, int, ConvectionDiffusionProblem);

    private:
    void make_step(int) override;
};

class RightAngleSolver : public Solver {
    public:
    RightAngleSolver(int, int, int, int, ConvectionDiffusionProblem);

    private:
    void make_step(int) override;
};

class ImplicitAngleSolver : public Solver {
    public:
    ImplicitAngleSolver(int, int, int, int, ConvectionDiffusionProblem);

    private:
    double d;
    void make_step(int) override;
};

class FourPointSolver : public Solver {
    public:
    FourPointSolver(int, int, int, int, ConvectionDiffusionProblem);

    private:
    double y_cur_left[2];
    double y_prev_left[2];
    double d;
    void make_step(int) override;
};

class LaxSolver : public Solver {
    public:
    LaxSolver(int, int, int, int, ConvectionDiffusionProblem);

    private:
    double d_1;
    double d_2;
    void make_step(int) override;
};

class LaxWendroffSolver : public Solver {
    public:
    LaxWendroffSolver(int, int, int, int, ConvectionDiffusionProblem);

    private:
    double d;
    void make_step(int) override;
};