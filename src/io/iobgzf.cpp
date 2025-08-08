#include "iobgzf.h"

namespace ngslib {
    void BGZFile::_open(const std::string &fn, const std::string mode) {
        _fname = fn; 
        _mode = mode;

        _bgzf = bgzf_open(fn.c_str(), mode.c_str());
        if (!_bgzf) {
            throw std::runtime_error("[iobgzf.cpp::BGZFile:_open] Failed to open file: " + fn);
        }
    }

    void BGZFile::close() {
        if (_bgzf) {
            int is_cl = bgzf_close(_bgzf);  // `bgzf_close` return 0 on success and -1 on error
            if (is_cl < 0) std::cerr << "[iobgzf.cpp::~BGZFile] " + _fname + " fail close." << std::endl;

            _bgzf = nullptr;  // must set to be a nullptr after bgzf_close
            _fname.clear();
            _mode.clear();
        }
    }

    // Add static factory method for creating multiple BGZFile objects
    // CAUSION: The key word of 'static' doest not need to be placed here.
    std::vector<std::unique_ptr<BGZFile>> BGZFile::open_multiple(
        const std::vector<std::string>& filenames, 
        const char* mode) 
    {
        std::vector<std::unique_ptr<BGZFile>> files;
        files.reserve(filenames.size());
        for (const auto& filename : filenames) {
            files.push_back(std::make_unique<BGZFile>(filename, mode));
        }

        return files;
    }

    // Write operations
    BGZFile& BGZFile::write(const std::string &data) {
        if (!_bgzf || _mode.find('w') == std::string::npos) {
            throw std::runtime_error("[iobgzf.cpp::BGZFile:write] File not open for writing");
        }
        
        ssize_t bytes_written = bgzf_write(_bgzf, data.c_str(), data.length());
        if (bytes_written < 0 || static_cast<size_t>(bytes_written) != data.length()) {
            throw std::runtime_error("[iobgzf.cpp::BGZFile:write] fail to write data: " + data);
        }

        return *this;
    }

    // Read operations
    BGZFile& BGZFile::reado(std::string &line) {
        // 用这个函数虽然可以和运算符匹配但坏处是不知道文件被读完了没有，这很糟糕
        bool is_eof = read(line, '\n');
        return *this;
    }

    BGZFile& BGZFile::read_bytes(std::string &data, size_t size) {
        if (!_bgzf || _mode.find('r') == std::string::npos) {
            throw std::runtime_error("[iobgzf.cpp::BGZFile:read] File not open for reading");
        }

        std::vector<char> buffer(size);
        ssize_t bytes_read = bgzf_read(_bgzf, buffer.data(), size);
        if (bytes_read < 0) {
            throw std::runtime_error("Read failed");
        }
        data.assign(buffer.data(), bytes_read);
        return *this;
    }

    // Read data by _delim_. default is a complete line
    bool BGZFile::read(std::string &data, char delim) {
        if (!_bgzf || _mode.find('r') == std::string::npos) {
            throw std::runtime_error("File not open for reading");
        }

        kstring_t s; s.l = s.m = 0; s.s = nullptr; // initialize kstring
        int ret = bgzf_getline(_bgzf, delim, &s);
        if (ret >= 0) {
            data.assign(s.s, s.l);  // assign the read data to the string
            free(s.s);              // free the allocated memory
            return true;
        } else {
            data.clear();  // clear data if hit the end of file
        }
        
        if (s.s) free(s.s);
        return false;  // EOF or error
    }

    bool BGZFile::readline_with_index(tbx_t* tbx, hts_itr_t* itr, std::string& line) {
        kstring_t s; s.s = nullptr; s.l = s.m = 0; // must be refreshed in loop
        int ret = tbx_bgzf_itr_next(_bgzf, tbx, itr, &s);
        if (ret < 0) {
            line.clear();
            free(s.s);
            return false;
        }

        line.assign(s.s, s.l);
        free(s.s);
        return true;
    }
}  // namespace ngslib