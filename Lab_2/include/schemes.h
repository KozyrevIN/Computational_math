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
    void make_step(int) override;
};

