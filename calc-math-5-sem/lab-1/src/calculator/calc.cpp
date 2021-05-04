#include "derivatives-inl.h"

#include <cmath>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <map>
#include <string>

const double       POINT     = 0.5;
const unsigned int PLOT_FROM = 0;
const unsigned int PLOT_TO   = 21;

int main() {
    std::map<std::string, std::pair<diff::Func, diff::Func>> functions{
        {"sin^2(x)",
         {[](double x) { return pow(sin(x), 2); }, [](double x) { return 2 * sin(x) * cos(x); }}},
        {"cos(sin(x))",
         {[](double x) { return cos(sin(x)); }, [](double x) { return (-cos(x)) * sin(sin(x)); }}},
        {"e^sin(cos(x))",
         {[](double x) { return exp(sin(cos(x))); },
          [](double x) { return sin(x) * cos(cos(x)) * -1 * exp(sin(cos(x))); }}},
        {"ln(x + 3)", {[](double x) { return log(x + 3); }, [](double x) { return 1 / (x + 3); }}},
        {"sqrt(x + 3)",
         {[](double x) { return sqrt(x + 3); }, [](double x) { return 1 / (2 * sqrt(x + 3)); }}}};

    std::ofstream out{};

    for (const auto& pair : functions) {
        out.open("output/calculator/" + pair.first + ".csv");

        if (!out.is_open()) {
            std::cerr << "cannot open file output/calculator/" << pair.first << " to write!\n";
            continue;
        }

        out << "h, first order direct, first order inverse, first order, second order, third "
               "order\n";

        for (unsigned int i = PLOT_FROM; i < PLOT_TO; i++) {
            double step          = 1 / pow(2, i);
            double derivative_an = pair.second.second(POINT);

            out << step;

            out << ", " << std::scientific
                << fabs(diff::GetDerivativeFirstOrderDirect(pair.second.first, POINT, step) -
                        derivative_an);

            out << ", " << std::scientific
                << fabs(diff::GetDerivativeFirstOrderInverse(pair.second.first, POINT, step) -
                        derivative_an);

            out << ", " << std::scientific
                << fabs(diff::GetDerivativeFirstOrderSimmetrical(pair.second.first, POINT, step) -
                        derivative_an);

            out << ", " << std::scientific
                << fabs(diff::GetDerivativeSecondOrder(pair.second.first, POINT, step) -
                        derivative_an);

            out << ", " << std::scientific
                << fabs(diff::GetDerivativeThirdOrder(pair.second.first, POINT, step) -
                        derivative_an);

            out << "\n";
        }

        out.close();
    }

    return 0;
}