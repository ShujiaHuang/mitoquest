// Author: Shujia Huang
// Date: 2021-08-18
#ifndef __INCLUDE_NGSLIB_UTILS_H__
#define __INCLUDE_NGSLIB_UTILS_H__

#include <iostream>
#include <fstream>
#include <unistd.h>
#include <filesystem>  // Available since C++17
#include <algorithm>
#include <sstream>
#include <string>
#include <cstring>    // use the 'strlen' function
#include <vector>
#include <set>
#include <iterator>
#include <cstdint>    // uint32_t

namespace ngslib {

    // define data types for variant calling
    struct GenomeRegion {
        std::string chrom; // chromosome name
        uint32_t start;    // 0-based start position
        uint32_t end;      // 0-based end position (exclusive)

        GenomeRegion() : chrom(""), start(0), end(0) {};
        GenomeRegion(const std::string& rid, uint32_t s, uint32_t e) : chrom(rid), start(s), end(e) {
            if (start > end) {
                throw std::invalid_argument("[ERROR] start postion is larger than end position in "
                                            "GenomeRegion: " + chrom + ":" + std::to_string(start) + 
                                            "-" + std::to_string(end));
            }
        };

        std::string to_string() const {
            return chrom + ":" + std::to_string(start) + "-" + std::to_string(end);
        }
    };

    // Get a path that unambiguously identifies the location of a file.
    std::string abspath(const std::string file_path);

    // Get stem directory and return
    std::string dirname(const std::string file_path);
    
    // Get file name from a path
    std::string basename(const std::string file_path);

    // Get filename suffix
    std::string suffix_name(const std::string file_path);
    
    // remove the filename extension
    std::string remove_filename_extension(const std::string filename);

    // Template function can only be defined in C++ header file
    template<typename T>
    std::string tostring(const T &value) {
        std::ostringstream ss;
        ss << value;
        return ss.str();
    }

    /** Check if a file is readable and exists.
     * 
     * @param name Name of a file to test
     * @return a bool type for file is readable and exists or not.
     * 
     */
    bool is_readable(const char *name);
    inline bool is_readable(const std::string &name) { return is_readable(name.c_str()); }

    // check a fold path is empty or not
    bool path_exists_and_not_empty(std::string folder_path);

    /**
     * @brief Make a folder if it doesn't exist, handling concurrent race conditions.
     * @param path  The directry path
     * 
     */
    bool safe_mkdir(std::string folder_path);
    bool safe_remove(std::string file_path);

    /**
     * @brief Get the last modification file object
     * 
     * @param directory_path 
     * @return std::string 
     * 
     */
    std::string get_last_modification_file(std::string directory_path);

    /**
     * @brief Get the unique strings vector: sort by length and then by ASCII
     * 
     * @param strings 
     * @return std::vector<std::string> 
     */
    std::vector<std::string> get_unique_strings(const std::vector<std::string>& strings);

    /**
     * @brief Get the first column from a file, this is used for getting filename from input filelist.
     * 
     * @param fn 
     * @return std::vector<std::string> 
     * 
     */
    std::vector<std::string> get_firstcolumn_from_file(const std::string &fn);

    template<typename T>
    std::string join(const std::vector<T> &input, const std::string delim="\t") {
        if (input.empty()) 
            return "";

        std::string out_str = tostring(input[0]);
        for (size_t i(1); i < input.size(); ++i) {
            out_str += (delim + tostring(input[i]));
        }

        return out_str;
    }

    void split(const std::string &in_str, std::vector<std::string> &out, const char *delim, bool is_append=false);
    
    // 重载了除 vector<string> 之外的所有其他基础数据类型(string 的操作和它们不同)，包括：int，double，float
    template<typename T>
    void split(const std::string &in_str, std::vector<T> &out, const char *delim, bool is_append=false) {
        
        if (!is_append) { out.clear(); }

        // Takes only space separated C++ strings when using 'stringstream'  
        std::istringstream ss;

        size_t i(0), find_start_idx(0), delim_len(std::strlen(delim)), len;
        std::string tok;
        T d;

        while(i != std::string::npos) {

            ss.clear();  // clear the stringstream pipe first

            i = in_str.find(delim, find_start_idx);
            len = (i==std::string::npos) ? in_str.length() - find_start_idx : i - find_start_idx;
            tok = in_str.substr(find_start_idx, len); 

            if (!tok.empty()) {
                ss.str(tok);
                ss >> d;
                out.push_back(d);
            } else {
                out.push_back(0);  // Empty value set to be 0. [int, double, float]
            }

            find_start_idx = i + delim_len;  // skip 'delim' character and set to find next 'delim' string
        }

        return;
    }

    // Template function to slice a vector from range X to Y
    template <typename T>
    std::vector<T> vector_slicing(const std::vector<T> &v, size_t x, size_t y) {
    
        if (x > y) 
            throw std::invalid_argument(
                "[ERROR] input value error in 'std::vector<T> vector_slicing"
                "(const std::vector<T> &v, int x, int y)', start (x) larger "
                "than end(y)!");

        // Begin and End iterator
        auto first = x < v.size() ? v.begin()+x : v.end();
        auto last  = y < v.size() ? v.begin()+y : v.end();
    
        // Copy the element: [first, last)
        std::vector<T> new_v(first, last);

        // Return the results
        return new_v;
    }

    // Template function to discover and return all duplicated elements as a vector
    template <typename T>
    std::vector<T> find_duplicates(std::vector<T> input) { // do not pass reference or pointer to the 'input'

        std::sort(input.begin(), input.end());  // The order of 'input' will be changed here
        std::vector<T> duplicates;
        for (size_t i(1); i < input.size(); ++i) {
            if (input[i] == input[i-1]) {
                if (duplicates.empty() || (input[i] != duplicates.back())) {
                    duplicates.push_back(input[i]);
                }
            }
        }

        return duplicates;
    }
    
}  // namespace ngslib

#endif
