// The C++ codes for VCF/BCF file reading
// Author: Shujia Huang
// Date: 2025-04-18
#include "vcf.h"

namespace ngslib {
    // Private method to open file and read header
    void VCFFile::_open(const std::string &fn, const std::string mode) {
        _fname = fn;
        _mode = mode;

        // Ensure mode is for reading
        if (_mode.empty()) {
            throw std::runtime_error("Read mode must be Provided: " + _mode);
        }

        _fp = hts_open(fn.c_str(), mode.c_str());
        if (_fp == nullptr) {
            close(); // Clean up potentially partially initialized state
            throw std::runtime_error("Failed to open VCF/BCF file: " + fn);
        }

        if (_mode[0] == 'r') {
            // Use htslib to read the header
            // Note: bcf_hdr_read returns a pointer to a bcf_hdr_t structure
            // which is managed by htslib. We need to transfer ownership to our VCFHeader.
            bcf_hdr_t *hdr_raw = bcf_hdr_read(_fp);
            if (hdr_raw == nullptr) {
                close(); // Clean up
                throw std::runtime_error("Failed to read VCF/BCF header from: " + fn);
            }
            // Transfer ownership to the VCFHeader object
            _hdr = VCFHeader(hdr_raw); // Use the private constructor taking bcf_hdr_t*
            if (!_hdr.is_valid()) {
                close(); // Clean up
                throw std::runtime_error("Invalid VCF/BCF header read from: " + fn);
            }

        } else if (_mode[0] == 'w') {
            // For writing, we need to ensure the header is valid
            if (!_hdr.is_valid()) {
                close(); // Clean up
                throw std::runtime_error("Cannot write VCF/BCF: Header is invalid.");
            }
            // Write the header immediately after opening
            if (bcf_hdr_write(_fp, _hdr.hts_header()) != 0) {
                close(); // Clean up
                throw std::runtime_error("Failed to write VCF/BCF header to: " + fn);
            }

        } else {
            close(); // Clean up
            throw std::runtime_error("Invalid mode for VCF/BCF file: " + _mode);
        }

        _io_status = 0; // Initial status is OK after successful open
    }

    // Destructor
    VCFFile::~VCFFile() {
        close();
    }

    // Explicit cleanup
    void VCFFile::close() {
        if (_itr) {
            bcf_itr_destroy(_itr);
            _itr = nullptr;
        }
        if (_idx) {
            hts_idx_destroy(_idx); // Use generic hts_idx_destroy
            _idx = nullptr;
        }
        // _hdr manages its own memory via shared_ptr

        if (_fp) {
            hts_close(_fp);
            _fp = nullptr;
        }

        _io_status = -1; // Mark as closed/invalid state
        _fname = "";
        _mode  = "";
    }

    // Load index
    void VCFFile::index_load() {
        if (!is_open()) {
            throw std::runtime_error("Cannot load index: VCF/BCF file is not open.");
        }
        if (_idx) { // Already loaded
            return;
        }

        // hts_idx_load2 tries both .tbi and .csi
        _idx = hts_idx_load2(_fname.c_str(), nullptr); // Pass NULL for default index filename
        if (_idx == nullptr) {
            // Try legacy tbx_index_load for VCF specifically (might find .tbi)
            // Note: bcf_index_load is a wrapper around hts_idx_load nowadays
            // _idx = tbx_index_load(_fname.c_str()); // Deprecated approach

            // If hts_idx_load2 failed, it's unlikely tbx_index_load would succeed,
            // unless the index has a very non-standard name not handled by hts_idx_load2.
            // Stick with hts_idx_load2 result.
             throw std::runtime_error("Failed to load VCF/BCF index for: " + _fname + ". Ensure index file (.tbi or .csi) exists.");
        }
    }

    // Fetch by region string
    bool VCFFile::fetch(const std::string &region) {
        if (!is_open() || !_idx) {
            _io_status = -2; // Indicate error state (e.g., index not loaded)
            return false;
        }
        if (_itr) { // Destroy previous iterator if exists
            bcf_itr_destroy(_itr);
            _itr = nullptr;
        }

        // Use bcf_itr_querys for string region parsing
        // Pass writable header pointer as the underlying hts_itr_querys expects void*
        _itr = bcf_itr_querys(_idx, _hdr.hts_header_writable(), region.c_str());
        if (_itr == nullptr) {
            // Query failed (e.g., invalid region string, contig not found)
            _io_status = -3; // Indicate query failure
            return false;
        }
        
        _io_status = 0; // Reset status for reading from iterator
        return true;
    }

    // Fetch by coordinates
    bool VCFFile::fetch(const std::string &seq_id, hts_pos_t beg, hts_pos_t end) {
         if (!is_open() || !_idx) {
            _io_status = -2;
            return false;
        }
        if (_itr) {
            bcf_itr_destroy(_itr);
            _itr = nullptr;
        }

        // Get the sequence ID (rid) from the header
        int rid = bcf_hdr_name2id(_hdr.hts_header(), seq_id.c_str());
        if (rid < 0) {
             _io_status = -4; // Contig not found in header
             return false;
        }

        // Use bcf_itr_queryi for coordinate-based query
        _itr = bcf_itr_queryi(_idx, rid, beg, end);
        if (_itr == nullptr) {
            // Query failed (e.g., invalid coordinates)
            _io_status = -3;
            return false;
        }
        
        _io_status = 0;
        return true;
    }

    // Read next record
    int VCFFile::read(VCFRecord &rec) {
        if (!is_open()) {
            _io_status = -1; // Not open
            return _io_status;
        }

        // Ensure the target record object has an allocated bcf1_t
        // If the passed VCFRecord is default-constructed (invalid), allocate one.
        if (!rec.is_valid()) {
             bcf1_t *b = bcf_init();
             if (!b) {
                 _io_status = -5; // Allocation failure
                 throw std::runtime_error("Failed to allocate memory for VCF record.");
             }
             // Assign the newly allocated record to the VCFRecord object
             rec = VCFRecord(b); // This transfers ownership to the shared_ptr in rec
        } else {
            // If the record is already valid, clear it for reuse.
            // This assumes the user might pass the same VCFRecord object repeatedly.
            rec.clear();
        }

        if (_itr) {
            // Read using iterator
            _io_status = bcf_itr_next(_fp, _itr, rec.hts_record_writable());
        } else {
            // Read sequentially
            // Signature: int bcf_read(htsFile *fp, const bcf_hdr_t *hdr, bcf1_t *rec)
            _io_status = bcf_read(_fp, _hdr.hts_header(), rec.hts_record_writable());
        }

        // bcf_read returns >= 0 on success, -1 on EOF, <-1 on error.
        // bcf_itr_next returns >= 0 on success, -1 on EOF, <-1 on error.
        return _io_status;
    }

    // Write record
    int VCFFile::write(VCFRecord &rec) {
        if (!is_open()) {
            _io_status = -1; // Not open
            return _io_status;
        }
        if (!rec.is_valid()) {
            _io_status = -2; // Invalid record passed
            return _io_status;
        }

        // Write the record using bcf_write
        // Signature: bcf_write(htsFile *fp, const bcf_hdr_t *hdr, bcf1_t *rec)
        _io_status = bcf_write(_fp, _hdr.hts_header(), rec.hts_record_writable());

        // bcf_write returns 0 on success, negative on error.
        return _io_status;
    }

    // Stream output operator
    std::ostream &operator<<(std::ostream &os, const VCFFile &r) {
        os << "VCFFile(file='" << r._fname << "', mode='" << r._mode
           << "', is_open=" << (r.is_open() ? "true" : "false")
           << ", status="   << r._io_status << ")";
        return os;
    }

} // namespace ngslib