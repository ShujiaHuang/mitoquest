/**
 * @file bgzf.h
 * 
 * @brief C++ wrapper for the codes of htslib/bgzf.h
 * Functions that read and write block gzipped files.
 * 
 * @author Shujia Huang
 * @date 2021-08-20
 * 
 */
#ifndef __INCLUDE_NGSLIB_IOBGZF_H__
#define __INCLUDE_NGSLIB_IOBGZF_H__

#include <iostream>
#include <sstream>  // Add this header for std::ostringstream
#include <string>
#include <vector>

#include <htslib/bgzf.h>
#include <htslib/kstring.h>


namespace ngslib {

    class BGZFile {
    private:
        std::string _fname;  // input file name
        std::string _mode;   // read and write mode, Mode matching.
        BGZF *_bgzf;         // file handle
        static const size_t DEFAULT_BUFFER_SIZE = 4096;

        // Prevent copying
        BGZFile(const BGZFile &b) = delete;             // reject using copy constructor (C++11 style).
        BGZFile &operator=(const BGZFile &b) = delete;  // reject using copy/assignment operator (C++11 style).

    public:
        BGZFile(): _bgzf(nullptr) {}  // default constructor, do nothing
        explicit BGZFile(const std::string &fn, const std::string mode = "rb") {
            open(fn, mode);  /* Open the specified file for reading or writing. */
        }
        ~BGZFile() { close(); /* call to close the file.*/ }

        BGZFile& write(const std::string &data); // Write operations
        BGZFile& read(std::string &data, size_t size = DEFAULT_BUFFER_SIZE); // Read Operations
        bool getline(std::string &line, char delim = '\n');

        // Utility methods
        /**
         * @brief call `bgzf_open` function to open file.
         * 
         * @param fn 
         * 
         * @param mode  Mode
         * The mode argument can be any of 'r', 'rb', 'a', 'ab', 'w', 'wb', 'uw', 'x', or 'xb' depending
         * on whether the file will be read or written.  
         * 
         * The default is the mode of fileobj if discernible; otherwise, the default is 'rb'.
         * A mode of 'r' is equivalent to one of 'rb', and similarly for 'w' and 'wb', 'a' and 
         * 'ab', and 'x' and 'xb'.
         * 
         */
        void open(const std::string &fn, const std::string mode);
        void close();
        bool is_open() const { return _bgzf != nullptr; }
        bool eof() const { return bgzf_check_EOF(_bgzf); }

        /**
         * @brief Force writing of all buffered data to disk.
         * 
         * This function ensures all buffered data is written to the underlying file.
         * It's useful when:
         * - You need to ensure data integrity
         * - Real-time writing is required
         * - Multiple processes are reading/writing the same file
         * 
         * @throw std::runtime_error if flush operation fails
         */
        void flush() { 
            if (bgzf_flush(_bgzf) < 0) {
                throw std::runtime_error("Flush failed");
            }
        }

        // Make stream operators friends
        friend BGZFile& operator>>(BGZFile& file, std::string& data);
        friend BGZFile& operator<<(BGZFile& file, const std::string& data);
        friend BGZFile& operator<<(BGZFile& file, const char* data);

        // Template version for arithmetic types
        template<typename T>
        friend typename std::enable_if<std::is_arithmetic<T>::value, BGZFile&>::type
        operator<<(BGZFile& file, const T& data);
        
        // For types with stream operator defined
        template<typename T>
        friend typename std::enable_if<!std::is_arithmetic<T>::value && 
                                       !std::is_convertible<T, const char*>::value &&
                                       !std::is_same<T, std::string>::value, BGZFile&>::type
        operator<<(BGZFile& file, const T& data);

    }; // class BGZFile

    // Non-member operator overloads
    inline BGZFile& operator>>(BGZFile& file, std::string& data) {
        return file.read(data);
    }

    inline BGZFile& operator<<(BGZFile& file, const std::string& data) {
        return file.write(data);
    }

    inline BGZFile& operator<<(BGZFile& file, const char* data) {
        return file.write(data);
    }

    template<typename T>
    typename std::enable_if<std::is_arithmetic<T>::value, BGZFile&>::type
    operator<<(BGZFile& file, const T& data) {
        return file.write(std::to_string(data));
    }

    template<typename T>
    typename std::enable_if<!std::is_arithmetic<T>::value && 
                            !std::is_convertible<T, const char*>::value &&
                            !std::is_same<T, std::string>::value, BGZFile&>::type
    operator<<(BGZFile& file, const T& data) {
        std::ostringstream ss;
        ss << data;
        return file.write(ss.str());
    }

    // Non-member operators for stream manipulators
    inline BGZFile& operator<<(BGZFile& file, std::ostream& (*manip)(std::ostream&)) {
        if (manip == static_cast<std::ostream& (*)(std::ostream&)>(std::endl)) {
            file.write("\n");
        }
        return file;
    }
}  // namespace ngslib

#endif

