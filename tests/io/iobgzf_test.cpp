#include <gtest/gtest.h>
#include <iostream>
#include <string>

#include "io/iobgzf.h"


TEST(iobgzfTest, BasicTest) {

    std::string temp_filename = "test.gz";
    ngslib::BGZFile file(temp_filename, "w");  // "wb" "uw"
    
    // Strings
    file << "Hello";               // 使用 const char* 重载
    file << std::string("World");  // 使用 string 重载
    
    // Numbers
    file << 42;                // 使用算术类型模板
    file << 3.14 << std::endl; // 使用算术类型模板
    file << "Hello" << " " << 42 << "\n" << std::endl;
    file.close();  // 确保数据被写入文件
    
    // // Custom type with operator<<
    // struct Point { 
    //     int x, y;
    //     friend std::ostream& operator<<(std::ostream& os, const Point& p) {
    //         return os << "(" << p.x << "," << p.y << ")";
    //     }
    // };
    // Point p{1, 2};
    // file << p;                 // 使用通用类型模板

    // Reading by lines
    ngslib::BGZFile infile(temp_filename, "rb");  // "rb"
    std::string line;
    std::cout << "Output: " << std::endl;
    infile.readline(line);
    std::cout << "-- Read line: " << line << std::endl;
    while (infile.readline(line)) {
        std::cout << "Read line: " << line << std::endl;
    }

    std::remove(temp_filename.c_str());
}

