// The C++ codes for VCF/BCF file reading
// Author: Shujia Huang
// Date: 2025-04-18 (Resumed 2025-04-27)

#ifndef __INCLUDE_NGSLIB_VCF_H__
#define __INCLUDE_NGSLIB_VCF_H__

#include <iostream>
#include <string>
#include <memory>
#include <stdexcept>

#include <htslib/vcf.h>
#include <htslib/hts.h>
#include <htslib/tbx.h> // For tbx_index_load with VCF

#include "vcf_header.h"
#include "vcf_record.h"
#include "utils.h"

namespace ngslib {

    /**
     * @brief A VCF/BCF file reader class.
     *
     * Wraps htslib functions for opening, reading, and iterating through
     * VCF or BCF files. Manages the file pointer, header, and iterators.
     */
    class VCFFile {
    private:
        std::string _fname;  // input file name
        std::string _mode;   // read mode
        int _io_status;      // I/O status code from the last read operation

        htsFile *_fp;        // htsFile file pointer for VCF/BCF
        hts_idx_t *_idx;     // VCF/BCF index pointer (usually tbi or csi)
        hts_itr_t *_itr;     // Iterator for querying specific regions
        VCFHeader _hdr;      // VCFHeader object managing the bcf_hdr_t

        // Private method to open the file and read the header
        void _open(const std::string &fn, const std::string mode);

        // Disable copy constructor and assignment operator
        VCFFile(const VCFFile &) = delete;
        VCFFile &operator=(const VCFFile &) = delete;

    public:
        /**
         * @brief Default constructor. Initializes an invalid reader.
         */
        VCFFile() : _fname(""), _mode(""), _io_status(-1), _fp(nullptr), _idx(nullptr), _itr(nullptr), _hdr() {}

        /**
         * @brief Constructor that opens a VCF/BCF file for reading.
         * @param fn Path to the VCF/BCF file. "-" for stdin.
         * @param mode Read mode (e.g., "r", "rb", "w", "wb", "wz"). Defaults to "r".
         * @throws std::runtime_error if the file cannot be opened or header read.
         */
        explicit VCFFile(const std::string &fn, const std::string mode = "r") :
            _fname(""), _mode(""), _io_status(-1), _fp(nullptr), _idx(nullptr), _itr(nullptr), _hdr()
        {
            _open(fn, mode);
        }

        /**
         * @brief Constructor that opens a VCF/BCF file for writing and writes the header.
         * @param fn Path to the output VCF/BCF file. "-" for stdout.
         * @param hdr The VCFHeader object containing the header to write.
         * @param mode Write mode (e.g., "w" for uncompressed VCF, "wb" for BCF, "wz" for bgzipped VCF). Defaults to "w".
         * @throws std::runtime_error if the file cannot be opened or the header written.
         */
        explicit VCFFile(const std::string &fn, VCFHeader &hdr, const std::string mode = "w") :
            _fname(""), _mode(""), _io_status(-1), _fp(nullptr), _idx(nullptr), _itr(nullptr), _hdr(hdr) // Copy/share the header
        {
            // Check if the provided header is valid before proceeding
            if (!hdr.is_valid()) {
                throw std::runtime_error("Provided VCFHeader is invalid.");
            }
            // _hdr is now managing the same underlying bcf_hdr_t as the input hdr
            _open(fn, mode);
        }

        /**
         * @brief Destructor. Closes file, frees index and iterator, header managed by VCFHeader.
         */
        ~VCFFile();

        /**
         * @brief Closes resources explicitly. Called by the destructor.
         */
        void close();

        /**
         * @brief Checks if the reader is associated with a valid open file.
         * @return True if the file pointer is not null, false otherwise.
         */
        bool is_open() const { return _fp != nullptr; }

        /**
         * @brief Gets the VCF header.
         * @return A reference to the VCFHeader object.
         */
        VCFHeader &header() { return _hdr; }

        /**
         * @brief Gets the VCF header (const version).
         * @return A const reference to the VCFHeader object.
         */
        const VCFHeader &header() const { return _hdr; }

        /**
         * @brief Loads the index file (.tbi or .csi) for the VCF/BCF file.
         * The index file must exist at the expected location (<filename>.tbi or <filename>.csi).
         * @throws std::runtime_error if the index cannot be loaded.
         */
        void index_load();

        /**
         * @brief Sets up an iterator to fetch records from a specific genomic region.
         * Requires the index to be loaded first via index_load().
         * @param region Genomic region string (e.g., "chr1:1000-2000", "chrM").
         * @return True if the iterator was successfully created, false otherwise (e.g., region invalid, index not loaded).
         */
        bool fetch(const std::string &region);

        /**
         * @brief Sets up an iterator to fetch records from a specific genomic region.
         * Requires the index to be loaded first via index_load().
         * @param seq_id The sequence name (chromosome).
         * @param beg The 0-based start position.
         * @param end The 0-based end position (exclusive, or HTS_POS_MAX for end of contig).
         * @return True if the iterator was successfully created, false otherwise.
         */
        bool fetch(const std::string &seq_id, hts_pos_t beg, hts_pos_t end);

        /**
         * @brief Reads the next VCF/BCF record from the file or iterator.
         * If an iterator is active (after fetch()), reads the next record within the region.
         * Otherwise, reads the next record sequentially from the file.
         * @param rec The VCFRecord object to populate with the read data.
         * @return >= 0 on success (typically 0), -1 on end of file/iterator, < -1 on error.
         *         The return value is stored and accessible via io_status().
         */
        int read(VCFRecord &rec);

        /**
         * @brief Alias for read(). Reads the next VCF/BCF record.
         * @param rec The VCFRecord object to populate.
         * @return Result of the read operation (>= 0 success, -1 EOF, < -1 error).
         */
        int next(VCFRecord &rec) { return read(rec); }

        /**
         * @brief Writes a VCF/BCF record to the output file.
         * The record must be associated with the same header used to open the writer.
         * @param rec The VCFRecord object to write.
         * @return 0 on success, negative value on error.
         *         The return value is stored and accessible via io_status().
         */
        int write(VCFRecord &rec);

        /**
         * @brief Gets the status of the last read operation.
         * @return >= 0 on success, -1 on end of file/iterator, < -1 on error.
         */
        int io_status() const { return _io_status; }

        /**
         * @brief Boolean conversion operator. Checks if the last read was successful.
         * @return True if io_status() >= 0, false otherwise.
         */
        operator bool() const { return _io_status >= 0; }

        /**
         * @brief Friend function for stream output (basic info).
         */
        friend std::ostream &operator<<(std::ostream &os, const VCFFile &r);

    }; // class VCFFile

} // namespace ngslib

#endif // __INCLUDE_NGSLIB_VCF_H__
