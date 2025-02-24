#include "mt_variant_caller.h"

// MtVariantCaller implementation
MtVariantCaller::MtVariantCaller(const Config& config) : _config(config) {
    // load fasta
    reference = _config.reference_file;

    _get_calling_interval();
    print_calling_interval();

    // keep the order of '_samples_id' as the same as input bamfiles
    _get_sample_id_from_bam();
}

MtVariantCaller::~MtVariantCaller() {
    // 析构函数的实现，如果不需要特殊操作，可以为空
}

bool MtVariantCaller::run() {

    clock_t cpu_start_time = clock();
    time_t real_start_time = time(0);

    bool is_success = _caller_process();

    // Time information
    time_t now = time(0);
    std::string ct(ctime(&now));
    ct.pop_back();
    std::cout << "[INFO] " + ct + ". Done for Variant calling, "
              << difftime(now, real_start_time) << " (CPU time: "
              << (double)(clock() - cpu_start_time) / CLOCKS_PER_SEC << ") seconds elapsed."
              << std::endl;

    return is_success;
}

void MtVariantCaller::_get_sample_id_from_bam() {
    time_t real_start_time = time(0);

    // Loading sample ID in BAM/CRMA files from RG tag.
    if (_config.filename_has_samplename)
        std::cout << "[INFO] BaseVar'll load samples id from filename directly, becuase you set "
                     "--filename-has-samplename.\n";

    std::string samplename, filename;
    size_t si;
    for (size_t i(0); i < _config.bam_files.size(); ++i) {

        if ((i+1) % 1000 == 0)
            std::cout << "[INFO] loading "   << i+1 << "/" << _config.bam_files.size()
                      << " alignment files." << std::endl;

        if (_config.filename_has_samplename) {
            filename = ngslib::remove_filename_extension(ngslib::basename(_config.bam_files[i]));
            si = filename.find('.');
            samplename = si > 0 && si != std::string::npos ? filename.substr(0, si) : filename;
        } else {
            // Get sampleID from BAM header, a bit time-consuming.
            ngslib::BamHeader bh(_config.bam_files[i]);
            samplename = bh.get_sample_name();
        }

        if (!samplename.empty()) {
            _samples_id.push_back(samplename);
        } else {
            throw std::invalid_argument("[MtVariantCaller::_load_sample_id_from_bam] " +
                                        _config.bam_files[i] + " sample ID not found.\n");
        }
    }

    // check samples duplication
    std::vector<std::string> duplicate_samples = ngslib::find_duplicates(_samples_id);
    if (!duplicate_samples.empty()) {
        std::cout << "[WARNING] Find " << duplicate_samples.size() << " duplicated samples within "
                  << "the input bamfiles: " + ngslib::join(duplicate_samples, ",") + "\n"
                  << std::endl;
    }

    // Time information
    time_t now = time(0);
    std::string ct(ctime(&now));
    ct.pop_back();  // rm the trailing '\n' put by `asctime`
    std::cout << "[INFO] " + ct + ". Done for loading all samples' id from alignment files, "
              << difftime(now, real_start_time) << " seconds elapsed.\n"
              << std::endl;

    return;
}

void MtVariantCaller::_get_calling_interval() {

    std::vector<GenomeRegion> regions;
    if (!_config.calling_regions.empty()) {
        std::vector<std::string> rg_v;
        ngslib::split(_config.calling_regions, rg_v, ",");

        for (size_t i(0); i < rg_v.size(); ++i) {
            if (rg_v[i].length() == 0) continue; // ignore empty string
            regions.push_back(_make_genome_region(rg_v[i]));
        }
    } else {
        // Call the whole genome
        int n = reference.nseq();
        for (size_t i(0); i < n; ++i) {
            std::string ref_id = reference.iseq_name(i);
            regions.push_back(GenomeRegion(ref_id, 1, reference.seq_length(ref_id)));
        }
    }

    _calling_intervals.clear();
    // split region into small pieces by chunk_size
    for (size_t i(0); i < regions.size(); ++i) {
        if (regions[i].end < regions[i].start) {
            throw std::invalid_argument("[ERROR] start postion is larger than end position in "
                                        "-r/--regions " + regions[i].ref_id + ":" +
                                        std::to_string(regions[i].start) + "-" + std::to_string(regions[i].end));
        }

        uint32_t total_length = regions[i].end - regions[i].start + 1;
        for (uint32_t j(0); j < total_length; j += _config.chunk_size) {
            uint32_t start = regions[i].start + j;
            uint32_t end = std::min(regions[i].end, start + _config.chunk_size - 1);
            _calling_intervals.push_back(GenomeRegion(regions[i].ref_id, start, end));
        }
    }
    return;
}

MtVariantCaller::GenomeRegion MtVariantCaller::_make_genome_region(std::string gregion) {

    // Genome Region, 1-based
    std::string ref_id;
    uint32_t reg_start, reg_end;

    std::vector<std::string> gr;
    ngslib::split(gregion, gr, ":");     // get reference id
    ref_id = gr[0];

    if (gr.size() == 2) {                // 'start-end' or start
        std::vector<uint32_t> gs;
        ngslib::split(gr[1], gs, "-");   // get position coordinate
        reg_start = gs[0];
        reg_end = (gs.size() == 2) ? gs[1] : reference.seq_length(ref_id);
    } else {
        reg_start = 1;
        reg_end = reference.seq_length(ref_id);  // the whole ``ref_id`` length
    }

    if (reg_start > reg_end) {
        throw std::invalid_argument("[ERROR] start postion is larger than end position in "
                                    "-r/--regions " + gregion);
    }

    return GenomeRegion(ref_id, reg_start, reg_end);  // 1-based
}

void MtVariantCaller::print_calling_interval() {

    std::string ref_id;
    uint32_t reg_start, reg_end;
    std::cout << "---- Calling Intervals ----\n";
    for (size_t i(0); i < _calling_intervals.size(); ++i) {
        std::cout << i+1 << " - " << _calling_intervals[i].ref_id << ":"
                  << _calling_intervals[i].start << "-" << _calling_intervals[i].end << "\n";
    }
    std::cout << "\n";
    return;
}

bool MtVariantCaller::_caller_process() {

    // Get filepath and stem name
    std::string _bname = ngslib::basename(_config.output_file);
    size_t si = _bname.find(".vcf");
    std::string stem_bn = (si > 0 && si != std::string::npos) ? _bname.substr(0, si) : _bname;

    std::string outdir = ngslib::dirname(_config.output_file);
    std::string cache_outdir = outdir + "/cache_" + stem_bn;
    ngslib::safe_mkdir(cache_outdir);

    // Variant calling
    // 以区间为单位进行变异检测, 每个区间里先按照样本调用多线程，然后合并样本，多线程遍历位点并行处理
    bool is_success = true;
    std::vector<std::string> sub_vcf_files;
    for (size_t i(0); i < _calling_intervals.size(); ++i) {



        // Call variants in parallel
        std::string regstr = _calling_intervals[i].ref_id + "_" +
                             std::to_string(_calling_intervals[i].start) + "_" +
                             std::to_string(_calling_intervals[i].end);
        std::string sub_vcf_fn = cache_outdir + "/" + stem_bn + "." + regstr + ".vcf.gz";
        sub_vcf_files.push_back(sub_vcf_fn);

    }

    // Merge VCF
    std::string header = vcf_header_define(_config.reference_file, _samples_id);
    merge_file_by_line(sub_vcf_files, _config.output_file, header, true);

    const tbx_conf_t bf_tbx_conf = {1, 1, 2, 0, '#', 0};  // {preset, seq col, beg col, end col, header-char, skip-line}
    if ((ngslib::suffix_name(_config.output_file) == ".gz") &&           // create index
        tbx_index_build(_config.output_file.c_str(), 0, &bf_tbx_conf))   // file suffix is ".tbi"
    {
        throw std::runtime_error("tbx_index_build failed: Is the file bgzip-compressed? "
                                 "Check this file: " + _config.output_file + "\n");
    }

    return is_success;
}

