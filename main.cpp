#include <iostream>
#include "vect.h"
#include "particle.h"
#include "field.h"

int main() {
    int64_t n, ticks;
    double size, temp;
    std::cin >> n >> size >> ticks >> temp;
    field test(n, size, temp);
    //test.printInfo();
    test.makeTicks(ticks);
    //test.printInfo();
    return 0;
}