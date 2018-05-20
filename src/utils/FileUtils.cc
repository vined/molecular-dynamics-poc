#include <fstream>
#include <iostream>

#include "FileUtils.h"


std::ifstream openFile(std::string fileName) {
    std::ifstream file;
    file.open(fileName);

    if (file.is_open()) {
        return file;
    }

    std::cout << "Failed to open file: " << fileName << std::endl;
    throw 200;
}

