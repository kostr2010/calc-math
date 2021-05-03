#ifndef LIB_H_INCLUDED
#define LIB_H_INCLUDED

#include <array>
#include <cassert>
#include <cstddef>
#include <functional>
#include <iostream>
#include <utility>
#include <vector>

namespace diff_sheme {

typedef float point_t;

class Grid {
  public:
    Grid() = delete;
    explicit Grid(const Grid& grid) : x_i(grid.x_i), f_i(grid.f_i), step_(grid.step_) {
    }
    explicit Grid(const std::vector<point_t>& x, const std::vector<point_t>& f) {
        assert(x.size() == f.size());

        x_i = x;
        f_i = f;

        n_nodes_ = x.size();
        step_ = (x[n_nodes_ - 1] - x[0]) / n_nodes_;
    }
    explicit Grid(const std::vector<point_t>& x, const std::function<point_t(point_t)>& f) {
        x_i = x;

        n_nodes_ = x.size();
        step_ = (x[n_nodes_ - 1] - x[0]) / n_nodes_;

        f_i.reserve(n_nodes_);

        for (size_t i = 0; i < n_nodes_; i++) {
            f_i[i] = f(x[i]);
        }
    }
    explicit Grid(const std::pair<point_t, point_t>& segment, size_t n_nodes,
                  const std::function<point_t(point_t)>& f) {
        n_nodes_ = n_nodes;
        step_ = step_ = (segment.second - segment.first) / n_nodes;

        x_i.reserve(n_nodes_);
        f_i.reserve(n_nodes_);

        x_i[0] = segment.first;
        f_i[0] = f(segment.first);
        x_i[n_nodes_ - 1] = segment.second;
        f_i[n_nodes_ - 1] = f(segment.second);

        for (size_t i = 1; i < n_nodes_ - 1; i++) {
            x_i[i] = segment.first + i * step_;
            f_i[i] = f(x_i[i]);
        }
    }

    ~Grid() = default;

    inline point_t step() const {
        return step_;
    }

    inline size_t n_nodes() const {
        return n_nodes_;
    }

    point_t x(const size_t node) const {
        assert(node < n_nodes_);

        return x_i[node];
    }

    point_t f(const size_t node) const {
        assert(node < n_nodes_);

        return f_i[node];
    }

    std::ostream& dump(std::ostream& os) const {
        for (size_t i = 0; i < n_nodes_ - 1; i++) {
            os << x_i[i] << ", ";
        }
        os << x_i[n_nodes_ - 1] << '\n';

        for (size_t i = 0; i < n_nodes_ - 1; i++) {
            os << f_i[i] << ", ";
        }
        os << f_i[n_nodes_ - 1] << '\n';

        return os;
    }

  private:
    std::vector<point_t> x_i;
    std::vector<point_t> f_i;

    size_t n_nodes_;
    point_t step_;
};

class DiffEquation {
  public:
    DiffEquation() = delete;
    explicit DiffEquation(std::function<point_t(point_t)> k, std::function<point_t(point_t)> q,
                          std::function<point_t(point_t)> f,
                          const std::pair<point_t, point_t>& deltas,
                          const std::pair<point_t, point_t>& epsilons)
        : k_(k), q_(q), f_(f), deltas_(deltas), epsilons_(epsilons) {
    }

    ~DiffEquation() = default;

    std::vector<point_t> SolveFixedCoeff(const std::pair<point_t, point_t>& segment,
                                         const size_t n_nodes, const point_t at) const {
        const point_t k_fixed = k_(at);
        const point_t q_fixed = q_(at);
        const point_t f_fixed = f_(at);

        const point_t step = (segment.second - segment.first) / n_nodes;

        const point_t a_0 = k_fixed;
        const point_t b_0 = -k_fixed - deltas_.first * step;
        const point_t c_0 = 0.0;
        const point_t d_0 = -epsilons_.first * step;

        const point_t a_l = k_fixed;
        const point_t b_l = -2.0 * k_fixed - q_fixed * step * step;
        const point_t c_l = k_fixed;
        const point_t d_l = -f_fixed * step * step;

        const point_t a_L = 0;
        const point_t b_L = -k_fixed - deltas_.second * step;
        const point_t c_L = k_fixed;
        const point_t d_L = -epsilons_.second * step;

        std::vector<point_t> u = {};
        std::vector<point_t> alphas = {};
        std::vector<point_t> betas = {};

        u.reserve(n_nodes);
        alphas.reserve(n_nodes);
        betas.reserve(n_nodes);

        alphas[0] = -a_0 / b_0;
        betas[0] = d_0 / b_0;

        // forward (calculating coefficients)
        for (size_t i = 1; i < n_nodes - 1; i++) {
            const auto coeff = b_l + c_l * alphas[i - 1];
            alphas[i] = -a_l / coeff;
            betas[i] = (d_l - c_l * betas[i - 1]) / coeff;
        }

        u[n_nodes - 1] = (d_L - c_L * betas[n_nodes - 1]) / (b_L + c_L * alphas[n_nodes - 1]);

        // backward (calculating u[i])
        for (int i = n_nodes - 2; i >= 0; i--) {
            u[i] = alphas[i] * u[i + 1] + betas[i];
        }

        return u;
    }

    std::vector<point_t> Solve(const std::pair<point_t, point_t>& segment,
                               const size_t n_nodes) const {
        Grid f_i(segment, n_nodes, f_);
        Grid q_i(segment, n_nodes, q_);

        const auto step = f_i.step();
        const auto step_sq = step * step;

        // k_(l + 1/2) <=> k_(l+1)
        // k_(l - 1/2) <=> k_(l)
        // k_0 & k_L are absent
        Grid k_i({segment.first + step / 2, segment.second - step / 2}, n_nodes - 1, k_);
        const auto k_0 = k_(segment.first);
        const auto k_L = k_(segment.second);

        const auto a_0 = k_0;
        const auto b_0 = -k_0 - deltas_.first * step;
        const auto c_0 = 0.0;
        const auto d_0 = -epsilons_.first * step;

        const auto a_L = 0.0;
        const auto b_L = -k_L - deltas_.second * step;
        const auto c_L = k_L;
        const auto d_L = -epsilons_.second * step;

        const std::function<point_t(size_t)> a_l = [&k_i](size_t l) { return k_i.f(l + 1); };
        const std::function<point_t(size_t)> b_l = [&k_i, &q_i, step_sq](size_t l) {
            return q_i.f(l) * step_sq - k_i.f(l + 1) - k_i.f(l);
        };
        const std::function<point_t(size_t)> c_l = [&k_i](size_t l) { return k_i.f(l); };
        const std::function<point_t(size_t)> d_l = [&f_i, step_sq](size_t l) {
            return f_i.f(l) * step_sq;
        };

        // forward
        std::vector<point_t> alphas = {-a_0 / b_0};
        std::vector<point_t> betas = {d_0 / b_0};
        alphas.reserve(n_nodes);
        betas.reserve(n_nodes);

        // FIXME: REDO K_I TO MATCH OTHER GRIDS IN TERMS OF INDEX

        for (size_t l = 1; l < n_nodes - 1; l++) {
            const auto coeff = b_l(l) + c_l(l) * alphas[l - 1];
            alphas[l] = -a_l(l) / coeff;
            betas[l] = (d_l(l) - c_l(l) * betas[l - 1]) / coeff;
        }

        // backward
        std::vector<point_t> u = {};
        u.reserve(n_nodes);
        u[n_nodes - 1] = (d_L - c_L * betas[n_nodes - 2]) / (b_L + c_L * alphas[n_nodes - 2]);

        for (int i = n_nodes - 2; i >= 0; i--) {
            u[i] = alphas[i] * u[i + 1] + betas[i];
        }

        return u;
    }

  private:
    std::function<point_t(point_t)> k_;
    std::function<point_t(point_t)>& q_;
    std::function<point_t(point_t)>& f_;

    const std::pair<point_t, point_t>& deltas_;
    const std::pair<point_t, point_t>& epsilons_;
};

}; // namespace diff_sheme

#endif
