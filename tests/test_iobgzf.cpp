#include <iostream>
#include "io/iobgzf.cpp"

// Example usage
void example() {
    ngslib::BGZFile file("test.gz", "w");  // "wb" "uw"
    
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
    ngslib::BGZFile infile("test.gz", "rb");  // "rb"
    std::string line;
    std::cout << "Output: " << std::endl;
    infile.getline(line);
        std::cout << "-- Read line: " << line << std::endl;
    while (infile.getline(line)) {
        std::cout << "Read line: " << line << std::endl;
    }
}

int main(int argc, char* argv[]) {
    example();
    return 0;
}

