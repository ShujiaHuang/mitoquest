#include "bam.h"


namespace ngslib {

    void Bam::_open(const std::string &fn, const std::string mode, const std::string ref_fn) {
        if ((mode[0] == 'r') && (!is_readable(fn))) {
            throw std::runtime_error("[bam.cpp::Bam:_open] file not found - " + fn);
        }

        // Open a SAM/BAM/CRAM file and return samFile file pointer.
        if (is_cram(fn)) {

            if (ref_fn.empty()) {
                throw std::runtime_error("[bam.cpp::Bam:_open] Reference file is required for CRAM file: " + fn);
            }

            _fp = sam_open_format(fn.c_str(), "r", NULL);
            if (_fp) {
                if (hts_set_fai_filename(_fp, ref_fn.c_str()) != 0) {
                    sam_close(_fp);
                    throw std::runtime_error("[bam.cpp::Bam:_open] Fail to set reference file: " + ref_fn);
                }
            }
            
            _reference_path = ref_fn;

        } else {
            _fp = sam_open(fn.c_str(), mode.c_str()); 
        }

        if (!_fp) {
            throw std::runtime_error("[bam.cpp::Bam:_open] file open failure.");
        }

        _fname = fn;
        _mode = mode;
        _itr = NULL;
        _idx = NULL;

        // Check index
        this->index_load();
        _io_status = 0;  // reset, 0 means everything is fine.

        return;
    }

    // copy constructor
    Bam::Bam(const Bam &b) : _itr(NULL), _idx(NULL) {
        _open(b._fname, b._mode, b._reference_path);  // reopen the bamfile
    }

    Bam &Bam::operator=(const Bam &b) {
        destroy();
        _open(b._fname, b._mode, b._reference_path);  // reopen the bamfile
        return *this;
    }

    void Bam::destroy() {
        // sam_close function is an alias name of hts_close.
        if (_fp) sam_close(_fp);
        if (_idx) hts_idx_destroy(_idx);
        if (_itr) sam_itr_destroy(_itr);

        _io_status = -1;
    }

    samFile *Bam::fp() const {
        return _fp;
    }

    BamHeader &Bam::header() {
        if (!_hdr) {
            _hdr = BamHeader(_fp);  // call copy constructor of BamHeader.
        }
        return _hdr;
    }

    hts_idx_t *Bam::idx() {
        if (!_idx) {
            this->index_load();
        }
        return _idx;
    }

    void Bam::index_load() {

        if (_idx)
            return;

        _idx = sam_index_load(_fp, _fname.c_str());
        if (!_idx) {
            this->destroy();
            throw std::runtime_error("[bam.cpp::Bam:index_load] Unable to load index file, please "
                                     "ensure index has been created (.bai/.crai) for " + _fname);
        }
    }

    // fetch 这个函数在使用多线程的时候会不会发生问题？尝试多区间处理方式？
    // 特别是类成员参数, _fp/_idx 在并行处理时是否存在问题? (htslib/thread_pool.h 参考一下)
    // 最好不要在一份文件中做多线程，而是以文件为单位跑多线程，从而在根上避免？ (本 BaseVar 就是这样做)
    // Create a SAM/BAM/CRAM iterator for one region.
    bool Bam::fetch(const std::string &region) {

        if (!_idx) index_load();           // May not be thread safety?
        if (!_hdr) _hdr = BamHeader(_fp);  // If NULL, set BAM header to _hdr.

        // Reset a iterator, An iterator on success; NULL on failure
        if (_itr) sam_itr_destroy(_itr);
        _itr = sam_itr_querys(_idx, _hdr.h(), region.c_str());

        if (!_itr) {
            throw std::runtime_error("[bam.cpp::Bam:fetch] Fail to fetch "
                                     "the alignment data in : " + region);
        }
        return _itr != NULL;
    }

    bool Bam::fetch(const std::string &seq_name, hts_pos_t beg, hts_pos_t end) {

        if (!_idx) index_load();  // May not be thread safety?
        if (!_hdr) _hdr = BamHeader(_fp);  // If NULL, set BAM header to _hdr.

        // Reset a iterator, An iterator on success; NULL on failure
        if (_itr) sam_itr_destroy(_itr);
        _itr = sam_itr_queryi(_idx, _hdr.name2id(seq_name), beg, end);

        if (!_itr) {
            std::string region = seq_name + ":" + tostring(beg) + "-" + tostring(end);
            throw std::runtime_error("[bam.cpp::Bam:fetch] Fail to fetch "
                                     "the alignment data in: " + region);
        }

        return _itr != NULL;
    }

    // Get alignement read into BamRecord
    // 如何解决线程中关于共享变量的问题? (用智能指针？)
    int Bam::read(BamRecord &br) {

        // If NULL, initial the BAM header by _fp.
        if (!_hdr.h()) _hdr = BamHeader(_fp);

        if (!_itr) {
            _io_status = br.load_read(_fp, _hdr.h());
        } else {
            _io_status = br.next_read(_fp, _itr);
        }

        // Destroy BamRecord and set br to be NULL if fail to read data
        if (_io_status < 0) 
            br.destroy();

        return _io_status;
    }

    std::ostream &operator<<(std::ostream &os, const Bam &b) {

        if (b) {
            os << b._fname;
        }
        return os;
    }

}  // namespace ngslib

