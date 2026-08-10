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
#include <memory>  // std::unique_ptr
#include <stdexcept>
#include <cstring> // std::memcmp
#include <cstdint> // uint8_t, uint16_t, uint32_t
#include <array>   // std::array

#include <htslib/bgzf.h>
#include <htslib/tbx.h>
#include <htslib/kstring.h>


namespace ngslib {

    // -----------------------------------------------------------------------
    //  BGZF block-level utilities
    //
    //  Low-level helpers for working with raw BGZF blocks.  These are useful
    //  for block-level concatenation, integrity checking, and streaming
    //  operations that need to parse or copy individual BGZF blocks without
    //  going through the decompression/recompression path.
    // -----------------------------------------------------------------------

    /// Size (in bytes) of a BGZF block header.
    ///
    /// Every BGZF block begins with a fixed 18-byte header laid out as follows:
    ///
    ///   Offset  Size  Field
    ///   ──────  ────  ─────────────────────────────────────────────
    ///    0        3   gzip ID1/ID2/CM   (0x1f 0x8b 0x08)
    ///    3        1   Compression method (0x08 = deflate)
    ///    4        1   Flags (FLG)
    ///    5        4   MTIME (modification time)
    ///    9        1   XFL  (extra flags)
    ///   10        1   OS   (operating system)
    ///   11        2   XLEN = 6  (extra-field length)
    ///   13        2   Sub-field SI1:SI2 = 'B':'C'  (BGZF identifier)
    ///   15        2   Sub-field LEN = 2
    ///   17        2   BSIZE = total_block_size - 1  (little-endian uint16)
    ///   ──────  ────
    ///   Total:  18 bytes
    ///
    /// The BSIZE field at offset 16-17 encodes the total compressed block
    /// size minus one.  Use bgzf_block_bsize() to extract it.
    static constexpr int BGZF_HEADER_SIZE = 18;

    /// Size (in bytes) of the BGZF EOF block — a 28-byte sentinel that marks
    /// the end of a BGZF stream.  Intermediate EOF blocks must be stripped
    /// when concatenating BGZF files; bgzf_close() writes the final one.
    static constexpr int BGZF_EOF_SIZE = 28;

    /// The 28-byte BGZF EOF block (a minimal, empty BGZF block).
    static constexpr std::array<uint8_t, BGZF_EOF_SIZE> BGZF_EOF_BLOCK = {{
        0x1f, 0x8b, 0x08, 0x04, 0x00, 0x00, 0x00, 0x00,
        0x00, 0xff, 0x06, 0x00, 0x42, 0x43, 0x02, 0x00,
        0x1b, 0x00, 0x03, 0x00, 0x00, 0x00, 0x00, 0x00,
        0x00, 0x00, 0x00, 0x00
    }};

    /**
     * @brief Validate a BGZF block header (first BGZF_HEADER_SIZE bytes of a raw block).
     *
     * Checks the three gzip magic bytes (0x1f 0x8b 0x08) at the start of
     * the block header.
     *
     * @param header  Pointer to at least BGZF_HEADER_SIZE bytes of raw BGZF block data.
     * @return 0 if the header is valid, -2 otherwise.
     */
    inline int bgzf_block_is_valid(const uint8_t *header) {
        if (header[0] != 0x1f || header[1] != 0x8b || header[2] != 0x08)
            return -2;
        return 0;
    }

    /**
     * @brief Extract BSIZE from a BGZF block header.
     *
     * BSIZE is stored as a 16-bit little-endian integer at offset 16 in the
     * block header (see BGZF_HEADER_SIZE layout).  Per the BGZF specification,
     * BSIZE = (total_block_size - 1).  Callers should add 1 to obtain the
     * actual block size in bytes.
     *
     * @param header  Pointer to at least BGZF_HEADER_SIZE bytes of raw BGZF block data.
     * @return Raw BSIZE value (actual block size minus 1).
     */
    inline uint16_t bgzf_block_bsize(const uint8_t *header) {
        return static_cast<uint16_t>(header[16]) |
               (static_cast<uint16_t>(header[17]) << 8);
    }

    /**
     * @brief Check whether a raw BGZF block is the 28-byte EOF sentinel.
     *
     * @param data  Pointer to the raw block data.
     * @param len   Length of the data in bytes.
     * @return true if the block is exactly the BGZF EOF block.
     */
    inline bool is_bgzf_eof(const uint8_t *data, size_t len) {
        return len == static_cast<size_t>(BGZF_EOF_SIZE) &&
               std::memcmp(data, BGZF_EOF_BLOCK.data(), BGZF_EOF_SIZE) == 0;
    }

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
            _fname(other._fname), 
            _mode(other._mode), 
            _bgzf(other._bgzf) 
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
        BGZFile& reado(std::string &line);       // Read operations (read one line per call)

        // Read up to _size_ bytes from the file storing into _data_.
        BGZFile& read_bytes(std::string &data, size_t size = DEFAULT_BUFFER_SIZE);
        bool read(std::string &data, char delim = '\n'); 
        bool readline(std::string &line) { return read(line, '\n'); } // Read a line from the file
        bool readline_with_index(tbx_t* tbx, hts_itr_t* itr, std::string& line); // Add method to read a line with index

        // Utility methods
        bool is_open() const { return _bgzf != nullptr; }
        bool eof() const { return _bgzf && bgzf_check_EOF(_bgzf); }
        /// Check whether the underlying file is BGZF/gzip compressed.
        /// Returns false for plain-text (uncompressed) files.
        bool is_compressed() const { return _bgzf && _bgzf->is_compressed; }
        void close();

        // Raw binary I/O for binary batchfile format (.bbf/.bbi)
        BGZF* bgzf_handle() { return _bgzf; }  // Access underlying BGZF handle

        // Write raw binary data (not text)
        ssize_t write_raw(const void *data, size_t len) {
            return bgzf_write(_bgzf, data, len);
        }

        // Read raw binary data (not text)
        ssize_t read_raw(void *data, size_t len) {
            return bgzf_read(_bgzf, data, len);
        }

        // -------------------------------------------------------------------
        //  BGZF block-level I/O (no decompression / recompression)
        //
        //  These methods operate on raw BGZF blocks directly, bypassing the
        //  decompression/recompression pipeline.  They are used for
        //  block-level concatenation, integrity checking, and streaming.
        // -------------------------------------------------------------------

        /// Read and decompress one BGZF block into the internal buffer.
        /// Returns 0 on success, -1 on error.
        int read_block() {
            return bgzf_read_block(_bgzf);
        }

        /// Number of bytes of decompressed data in the current block.
        int block_length() const {
            return _bgzf->block_length;
        }

        /// Pointer to the internal decompression buffer (read-only).
        const char* uncompressed_data() const {
            return static_cast<const char*>(_bgzf->uncompressed_block);
        }

        /// Read raw BGZF block bytes directly from the underlying file
        /// (no decompression).  Returns bytes read, 0 on EOF, -1 on error.
        ssize_t raw_block_read(void *data, size_t len) {
            return bgzf_raw_read(_bgzf, data, len);
        }

        /// Write raw BGZF block bytes directly to the underlying file
        /// (no compression).  Returns bytes written, or -1 on error.
        ssize_t raw_block_write(const void *data, size_t len) {
            return bgzf_raw_write(_bgzf, data, len);
        }

        // Seek to a BGZF virtual offset (for .bbi index-based random access)
        int seek_virtual(int64_t virtual_offset) {
            return bgzf_seek(_bgzf, virtual_offset, SEEK_SET);
        }

        // Get current BGZF virtual offset (for building .bbi index)
        int64_t tell_virtual() {
            return bgzf_tell(_bgzf);
        }

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
        // 用这个函数虽然可以和运算符匹配但坏处是不知道文件被读完了没有，这很糟糕，要考虑能否改进
        return file.reado(data);
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

