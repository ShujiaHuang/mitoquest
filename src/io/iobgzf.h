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
#include <htslib/tbx.h>
#include <htslib/kstring.h>


namespace ngslib {

    class BGZFile {
    private:
        static const size_t DEFAULT_BUFFER_SIZE = 4096;

        std::string _fname;  // input file name
        std::string _mode;   // read and write mode, Mode matching.
        BGZF *_bgzf;         // file handle

        // Private method to open the file
        void _open(const std::string &fn, const std::string mode);

        // Prevent copying
        BGZFile(const BGZFile &b) = delete;             // reject using copy constructor (C++11 style).
        BGZFile &operator=(const BGZFile &b) = delete;  // reject using copy/assignment operator (C++11 style).

    public:
        BGZFile(): _bgzf(nullptr) {}  // default constructor, do nothing

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
        explicit BGZFile(const std::string &fn, const std::string mode = "rb") {
            _open(fn, mode);  /* Open the specified file for reading or writing. */
        }

        // Add move constructor and move assignment
        BGZFile(BGZFile&& other) noexcept : 
            _bgzf(other._bgzf), 
            _fname(other._fname), 
            _mode(other._mode) 
        {
            other._bgzf = nullptr;
            other._fname.clear();
            other._mode.clear();
        }
        
        BGZFile& operator=(BGZFile&& other) noexcept {
            if (this != &other) {
                close();  // Close current file if open

                // Move all members
                _bgzf = other._bgzf;
                _fname = std::move(other._fname);
                _mode = std::move(other._mode);
                
                // Reset other object's state
                other._bgzf = nullptr;
                other._fname.clear();
                other._mode.clear();
            }
            return *this;
        }

        // Add static factory method for creating multiple BGZFile objects
        static std::vector<std::unique_ptr<BGZFile>> open_multiple(
            const std::vector<std::string>& filenames, 
            const char* mode = "r");

        // Destructor
        ~BGZFile() { close(); /* call to close the file.*/ }
    
        // Functions for I/O operations
        BGZFile& write(const std::string &data); // Write operations

        // Read up to _size_ bytes from the file storing into _data_.
        BGZFile& read_bytes(std::string &data, size_t size = DEFAULT_BUFFER_SIZE);
        bool read(std::string &data, char delim = '\n'); 
        bool readline(std::string &line) { return read(line, '\n'); } // Read a line from the file
        bool readline_with_index(tbx_t* tbx, hts_itr_t* itr, std::string& line); // Add method to read a line with index

        // Utility methods
        bool is_open() const { return _bgzf != nullptr; }
        bool eof() const { return bgzf_check_EOF(_bgzf); }
        void close();

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
        friend BGZFile& operator>>(BGZFile& file, std::string& data);       // read
        friend BGZFile& operator<<(BGZFile& file, const std::string& data); // write
        friend BGZFile& operator<<(BGZFile& file, const char* data);        // write

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
        return file.read_bytes(data);
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

