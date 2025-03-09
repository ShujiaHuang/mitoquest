#include "mt_variant_caller.h"

// MtVariantCaller implementation
MtVariantCaller::MtVariantCaller(int argc, char* argv[]) {

    Config config;
    // Set default values
    config.min_mapq                = 20;
    // config.min_baseq               = 20;
    config.heteroplasmy_threshold  = 0.01;
    config.thread_count            = 1;
    config.chunk_size              = 1000;
    config.pairs_map_only          = false;
    config.proper_pairs_only       = false;
    config.filename_has_samplename = false;

    static const struct option MT_CMDLINE_LOPTS[] = {
        {"reference",          required_argument, 0, 'R'},
        {"output",             required_argument, 0, 'o'},

        {"bam-list",           optional_argument, 0, 'b'},
        {"regions",            optional_argument, 0, 'r'},
        {"min-MQ",             optional_argument, 0, 'q'},
        // {"min-BQ",             optional_argument, 0, 'Q'},
        {"threads",            optional_argument, 0, 't'},
        {"threshold",          optional_argument, 0, 'j'},
        {"chunk",              optional_argument, 0, 'c'},

        {"pairs-map-only",          no_argument, 0, 'p'}, // 小写 p
        {"proper-pairs-only",       no_argument, 0, 'P'},
        {"filename-has-samplename", no_argument, 0, '1'},
        {"help",                    no_argument, 0, 'h'},

        // must set this value, to get the correct value from getopt_long
        {0, 0, 0, 0}
    };

    int opt;
    std::vector<std::string> bam_filelist;
    while ((opt = getopt_long(argc, argv, "f:b:o:r:q:c:j:t:pPh", MT_CMDLINE_LOPTS, NULL)) != -1) {
        switch (opt) {
            case 'f': config.reference_file = optarg;                    break;
            case 'o': config.output_file    = optarg;                    break;
            case 'b': 
                bam_filelist = ngslib::get_firstcolumn_from_file(optarg);
                config.bam_files.insert(
                    config.bam_files.end(), 
                    bam_filelist.begin(), 
                    bam_filelist.end()
                );
                break;

            case 'r': config.calling_regions         = optarg;            break;
            case 'q': config.min_mapq                = std::atoi(optarg); break;
            // case 'Q': config.min_baseq               = std::atoi(optarg); break;
            case 'p': config.pairs_map_only          = true;              break;
            case 'P': config.proper_pairs_only       = true;              break;
            case '1': config.filename_has_samplename = true;              break;
            case 'c': config.chunk_size              = std::atoi(optarg); break;
            case 'j': config.heteroplasmy_threshold  = std::atof(optarg); break;
            case 't': config.thread_count            = std::atoi(optarg); break;
            case 'h': usage(config); exit(0);
            
            default:
                std::cerr << "Unknown argument: " << opt << std::endl;
                exit(1);
        }
    }

    // Collect BAM/CRAM files
    while (optind < argc) {
        config.bam_files.push_back(argv[optind++]);
    }

    /* Make sure we set valid arguments */
    if (config.reference_file.empty() || config.output_file.empty()) {
        std::cerr << "Error: Missing required arguments. \n\n";
        usage(config);
        exit(1);
    }

    if (config.bam_files.empty()) {
        std::cerr << "Error: Missing required BAM/CRAM files.\n\n";
        usage(config);
        exit(1);
    }

    if (config.min_mapq < 0) {
        std::cerr << "Error: Quality score must be non-negative\n";
        exit(1);
    }

    if (config.heteroplasmy_threshold <= 0.0 || config.heteroplasmy_threshold > 1.0) {
        std::cerr << "Error: Heteroplasmy threshold must be between 0 and 1\n";
        exit(1);
    }

    if (config.thread_count < 1) {
        std::cerr << "Error: Thread count must be at least 1\n";
        exit(1);
    }

    if (config.chunk_size < 100) {
        std::cerr << "Error: Chunk size must be at least 100\n";
        exit(1);
    }

    // Output the commandline options "f:b:o:r:q:c:j:t:pPh"
    std::cout <<
        "[INFO] arguments: "
        "mitoquest caller -f " + config.reference_file + " "
        "-t " << config.thread_count           << " "
        "-q " << config.min_mapq               << " "
        "-j " << config.heteroplasmy_threshold << " "
        "-c " << config.chunk_size             << (config.filename_has_samplename ? " "
        "--filename-has-samplename" : "")      << (config.pairs_map_only ?          " "
        "--pairs-map-only"          : "")      << (config.proper_pairs_only ?       " "
        "--proper-pairs-only"       : "")      << " "
        "-o " + config.output_file + " "       << config.bam_files[0] + " [ ... " 
        << config.bam_files.size() << " bamfiles in total]. \n" << std::endl;

    // set parameters
    _config = config;

    // load fasta
    reference = _config.reference_file;
    _get_calling_interval(); // _print_calling_interval();

    // keep the order of '_samples_id' as the same as input bamfiles
    _get_sample_id_from_bam();
}

MtVariantCaller::~MtVariantCaller() {
    // 析构函数的实现，如果不需要特殊操作，则为空
}

void MtVariantCaller::usage(const Config &config) {
    std::cout << MITOQUEST_DESCRIPTION            << "\n"
              << "Version: " << MITOQUEST_VERSION << "\n\n"
            //   << "Author: "  << MITOQUEST_AUTHOR << " <" << MITOQUEST_AUTHOR_EMAIL << ">\n\n"

              << "Usage: mitoquest caller [options] -f ref.fa -o output.vcf.gz in1.bam [in2.bam ...]\n\n"
              << "Required options:\n"
              << "  -f, --reference FILE       Reference FASTA file\n"
              << "  -o, --output FILE          Output VCF file\n\n"

              << "Optional options:\n"
              << "  -b, --bam-list FILE        list of input BAM/CRAM filenames, one per line.\n"
              << "  -r, --regions REG[,...]    Comma separated list of regions in which to process (default: entire genome).\n"
              << "                             REG format: chr:start-end (e.g.: chrM or chrM:1-1000,chrM:8000-8200)\n"
              << "  -q, --min-MQ INT           skip alignments with mapQ smaller than INT (default: " << config.min_mapq  << ")\n"
            //   << "  -Q, --min-BQ INT           skip bases with base quality smaller than INT (default: " << config.min_baseq << ")\n"
              << "  -p, --pairs-map-only       Only use the paired reads which mapped to the some chromosome.\n"
              << "  -P, --proper-pairs-only    Only use properly paired reads.\n"
              << "  --filename-has-samplename  If the name of bamfile is something like 'SampleID.xxxx.bam', set this\n"
              << "                             argrument could save a lot of time during get the sample id from BAMfile.\n"
              << "  -j, --threshold FLOAT      Heteroplasmy threshold (default: " << config.heteroplasmy_threshold << ")\n"
              << "  -c, --chunk INT            Chunk size for parallel processing (default: " << config.chunk_size << ")\n"
              << "  -t, --threads INT          Number of threads (default: " << config.thread_count << ")\n"
              << "  -h, --help                 Print this help message.\n\n";
}

bool MtVariantCaller::run() {
    return _caller_process();
}

void MtVariantCaller::_get_sample_id_from_bam() {
    time_t real_start_time = time(0);

    // Loading sample ID in BAM/CRMA files from RG tag.
    if (_config.filename_has_samplename)
        std::cout << "[INFO] load samples id from filename directly, becuase you set "
                     "--filename-has-samplename\n";

    _samples_id.clear();

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
            ngslib::BamHeader bh(_config.bam_files[i], _config.reference_file);
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
    std::cout << "[INFO] " + ct + ". Done for loading all " + std::to_string(_samples_id.size())
              << " samples' id from alignment files, " << difftime(now, real_start_time)
              << " seconds elapsed.\n" << std::endl;

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
        uint32_t total_length = regions[i].end - regions[i].start + 1;
        for (uint32_t j(0); j < total_length; j += _config.chunk_size) {
            uint32_t start = regions[i].start + j;
            uint32_t end = std::min(regions[i].end, start + _config.chunk_size - 1);
            _calling_intervals.push_back(GenomeRegion(regions[i].ref_id, start, end));
        }
    }

    return;
}

GenomeRegion MtVariantCaller::_make_genome_region(std::string gregion) {

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

void MtVariantCaller::_print_calling_interval() {

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

    std::string outdir = ngslib::dirname(ngslib::abspath(_config.output_file));
    std::string cache_outdir = outdir + "/cache_" + stem_bn;
    ngslib::safe_mkdir(cache_outdir);

    // 以区间为单位进行变异检测, 每个区间里先按照样本调用多线程，然后合并样本，多线程遍历位点并行处理
    bool is_success = true;
    std::vector<std::string> sub_vcf_files;

    for (auto &gr: _calling_intervals) {
        // 按样本进行 pileup，在函数里按样本进行多线程，记录每个样本在每个位点上的 pileup 信息
        std::vector<PosVariantMap> samples_pileup_v;
        samples_pileup_v.reserve(_samples_id.size() + 1);  // reserve the memory before push_back

        //////////////////////////////////////////////
        bool is_empty = _fetch_base_in_region(gr, samples_pileup_v);
        if (is_empty) {
            std::cerr << "[WARNING] No reads in region: " << gr.ref_id << ":"
                      << gr.start << "-" << gr.end << "\n";
            continue;
        }

        //////////////////////////////////////////////
        // Call variants in parallel
        std::string regstr = gr.ref_id + "_" + std::to_string(gr.start) + "_" + std::to_string(gr.end);
        std::string sub_vcf_fn = cache_outdir + "/" + stem_bn + "." + regstr + ".vcf.gz";
        sub_vcf_files.push_back(sub_vcf_fn);

        is_empty = _variant_discovery(gr, samples_pileup_v, sub_vcf_fn);
        if (is_empty) {
            std::cout << "[INFO] No variants in region: " << gr.ref_id << ":"
                      << gr.start << "-" << gr.end << "\n";
        }
    }

    // Merge VCF
    std::string header = vcf_header_define(_config.reference_file, _samples_id);
    merge_file_by_line(sub_vcf_files, _config.output_file, header, IS_DELETE_CACHE);

    const tbx_conf_t bf_tbx_conf = {1, 1, 2, 0, '#', 0};  // {preset, seq col, beg col, end col, header-char, skip-line}
    if ((ngslib::suffix_name(_config.output_file) == ".gz") &&          // create index
        tbx_index_build(_config.output_file.c_str(), 0, &bf_tbx_conf))  // file suffix is ".tbi"
    {
        throw std::runtime_error("tbx_index_build failed: Is the file bgzip-compressed? "
                                 "Check this file: " + _config.output_file + "\n");
    }

    if (IS_DELETE_CACHE) {
        for (auto fn: sub_vcf_files) {
            ngslib::safe_remove(fn);
        }

        ngslib::safe_remove(cache_outdir);
    }

    return is_success;
}

bool MtVariantCaller::_fetch_base_in_region(const GenomeRegion gr, std::vector<PosVariantMap> &samples_pileup_v) {

    ThreadPool thread_pool(this->_config.thread_count);  // set multiple-thread

    std::vector<std::future<PosVariantMap>> pileup_results;
    pileup_results.reserve(this->_config.bam_files.size());

    std::string fa_seq = this->reference[gr.ref_id];     // use the whole sequence of ``ref_id`` for simply
    // Loop all alignment files
    for(size_t i(0); i < this->_config.bam_files.size(); ++i) { // The same order as this->_samples_id
        pileup_results.emplace_back(
            thread_pool.submit(call_pileup_in_sample, 
                               this->_config.bam_files[i], 
                               std::cref(fa_seq), 
                               gr,
                               std::cref(this->_config))
        );
    }

    samples_pileup_v.clear(); // clear the vector before push_back
    bool is_empty = true;
    for (auto && p: pileup_results) { // Run and make sure all processes are finished
        // 一般来说，只有当 valid() 返回 true 的时候才调用 get() 去获取结果，这也是 C++ 文档推荐的操作。
        if (p.valid()) {
            // get() 调用会改变其共享状态，不再可用，也就是说 get() 只能被调用一次，多次调用会触发异常。
            PosVariantMap pm = p.get();
            samples_pileup_v.push_back(pm);

            if (!pm.empty()) is_empty = false;
        }
    }

    if (samples_pileup_v.size() != this->_samples_id.size())
        throw std::runtime_error("[_fetch_base_in_region] 'samples_pileup_v.size()' "
                                 "should be the same as '_config.bam_files.size()'");
    return is_empty;  // no cover reads in 'GenomeRegion' if empty.
}

bool MtVariantCaller::_variant_discovery(const GenomeRegion gr, const std::vector<PosVariantMap> &samples_pileup_v, const std::string out_vcf_fn) {
    // 1. integrate the variant information of all samples in the region
    // 2. call the variant by the integrated information
    // 3. output the variant information to the VCF file
    ThreadPool thread_pool(this->_config.thread_count);  // set multiple-thread
    std::vector<std::future<VCFRecord>> results;

    for (uint32_t pos(gr.start); pos < gr.end + 1; ++pos) {
        std::vector<VariantInfo> vvi;
        vvi.reserve(samples_pileup_v.size()); // reserve the memory before push_back

        // get the variant information of all samples in the position
        PosVariantMap::const_iterator smp_pos_it;
        bool is_empty = true;
        for (size_t i(0); i < samples_pileup_v.size(); ++i) {

            smp_pos_it = samples_pileup_v[i].find(pos);
            if (smp_pos_it != samples_pileup_v[i].end()) {
                vvi.push_back(smp_pos_it->second);
                if (is_empty) is_empty = false;
            } else {
                VariantInfo vi(gr.ref_id, pos, 0, 0); // empty VariantInfo
                vvi.push_back(vi);
            }
        }

        // ignore the position which no reads cover of all samples
        if (!is_empty) { 
            // performance multi-thread here
            results.emplace_back(thread_pool.submit(call_variant_in_pos, vvi));  // return VCFRecord
        }
    }

    BGZF *out_vcf_fp = bgzf_open(out_vcf_fn.c_str(), "w");
    if (out_vcf_fp == NULL) {
        throw std::runtime_error("[ERROR] Cannot open file: " + out_vcf_fn);
    }

    bool is_empty = true;
    for (auto && p: results) { // Run and make sure all processes are finished
        if (p.valid()) {
            VCFRecord vcf_record = p.get();

            // write to a file
            if (vcf_record.is_valid()) {
                if (is_empty) is_empty = false;

                std::string vcf_str = vcf_record.to_string() + "\n";
                if (bgzf_write(out_vcf_fp, vcf_str.c_str(), vcf_str.length()) != vcf_str.length()) {
                    throw std::runtime_error("[ERROR] fail to write data");
                }
            }
        }
    }
    
    int is_cl = bgzf_close(out_vcf_fp);
    if (is_cl < 0) throw std::runtime_error("[ERROR] " + out_vcf_fn + " fail close.");
    
    return is_empty;  // no variant in the region if empty.
}

// Seek the base information of each sample in the region.
PosVariantMap call_pileup_in_sample(const std::string sample_bam_fn, 
                                    const std::string &fa_seq,
                                    const GenomeRegion gr,
                                    const MtVariantCaller::Config &config) 
{
    // The expend size of region, 100bp is enough.
    static const uint32_t REG_EXTEND_SIZE = 100;
    uint32_t extend_start     = gr.start > REG_EXTEND_SIZE ? gr.start - REG_EXTEND_SIZE : 1;
    uint32_t extend_end       = gr.end + REG_EXTEND_SIZE;
    std::string extend_regstr = gr.ref_id + ":" + std::to_string(extend_start) + "-" + std::to_string(extend_end);

    // 位点信息存入该变量, 且由于是按区间读取比对数据，key 值无需再包含 ref_id，因为已经不言自明
    PosMap sample_posinfo_map; // key: position, value: alignment information
    ngslib::Bam bf(sample_bam_fn, "r", config.reference_file); // open bamfile in reading mode (one sample per bamfile)
    if (bf.fetch(extend_regstr)) { // Set 'bf' only fetch alignment reads in 'exp_regstr'.
        hts_pos_t map_ref_start, map_ref_end;  // hts_pos_t is uint64_t
        std::vector<ngslib::BamRecord> sample_target_reads; 
        ngslib::BamRecord al;  // alignment read
        while (bf.next(al) >= 0) {  // -1 => hit the end of alignement file.
            if (al.mapq() < config.min_mapq || al.is_duplicate() || al.is_qc_fail() ||
                (al.is_paired() && config.proper_pairs_only && !al.is_proper_pair()))
            {
                continue;
            }

            if (al.is_paired() && config.pairs_map_only) {
                std::string tid_name(al.tid_name(bf.header()));
                std::string mate_tid_name(al.mate_tid_name(bf.header()));

                // only use the paired reads which mapped to the same chromosome
                if (tid_name != mate_tid_name) continue;
            }

            map_ref_start = al.map_ref_start_pos() + 1;  // al.map_ref_start_pos() is 0-based, convert to 1-based
            map_ref_end   = al.map_ref_end_pos();        // al.map_ref_end_pos() is 1-based

            // Only fetch reads which in [reg_start, reg_end]
            if (gr.start > map_ref_end) continue;
            if (gr.end < map_ref_start) break;

            sample_target_reads.push_back(al);  // record all the proper reads of sample
        }

        if (sample_target_reads.size() > 0) {
            // get alignment information of [i] sample.
            seek_position(fa_seq, sample_target_reads, gr, sample_posinfo_map);
        }
    }

    // 信息抽提：先计算并返回每个样本在该区间里每一个位点的潜在突变信息（类似 pileup or gvcf），若不如此处理，内存吃不消
    // 这样做的好处还可为多样本 joint-calling 打下基础.
    PosVariantMap sample_pileup_map;
    for (auto &pos_align_info: sample_posinfo_map) {
        VariantInfo vi = basetype_caller_unit(pos_align_info.second, config.heteroplasmy_threshold);

        // key: position, value: variant information, 由于是按区间抽提，key 值无需再包含 ref_id，因为已经不言自明 
        sample_pileup_map.insert({pos_align_info.first, vi});
    }

    // return the variant information for all the ref_pos of the sample in the region, no matter its a variant or not.
    return sample_pileup_map;
}

void seek_position(const std::string &fa_seq,   // must be the whole chromosome sequence
                   const std::vector<ngslib::BamRecord> &sample_map_reads,  // record the alignment reads of sample
                   const GenomeRegion gr,
                   PosMap &sample_posinfo_map)  // key: position, value: alignment information
{
    if (!sample_posinfo_map.empty())
        throw std::runtime_error("[seek_position] 'sample_posinfo_map' must be empty.");

    // A vector of: (cigar_op, read position, reference position, read base, read_qual, reference base)
    std::vector<ngslib::ReadAlignedPair> aligned_pairs;
    for (auto &al: sample_map_reads) {  // loop mapping reads
        AlignBase ab;
        ab.map_strand = al.map_strand();  // '*', '-' or '+'
        ab.mapq = al.mapq();

        // std::string ref_bases;
        uint32_t map_ref_pos;
        aligned_pairs = al.get_aligned_pairs(fa_seq);
        for (size_t i(0); i < aligned_pairs.size(); ++i) {
            map_ref_pos = aligned_pairs[i].ref_pos + 1;  // ref_pos is 0-based, convert to 1-based;

            if (gr.start > map_ref_pos) continue;
            if (gr.end < map_ref_pos) break;

            // 'BAM_XXX' are macros for CIGAR, which defined in 'htslib/sam.h'
            if (aligned_pairs[i].op == BAM_CMATCH ||  /* CIGAR: M */ 
                aligned_pairs[i].op == BAM_CEQUAL ||  /* CIGAR: = */
                aligned_pairs[i].op == BAM_CDIFF)     /* CIGAR: X */
            {
                // SNV. Only one character
                ab.ref_base  = aligned_pairs[i].ref_base[0];
                ab.read_base = aligned_pairs[i].read_base[0];
                ab.base_qual = aligned_pairs[i].read_qual[0];
            } else if (aligned_pairs[i].op == BAM_CINS) {  /* CIGAR: I */
                // Insertion
                if (!aligned_pairs[i].ref_base.empty()) {
                    std::cerr << al << "\n";
                    throw std::runtime_error("[ERROR] We got reference base in insertion region.");
                }

                // roll back one position to the leftmost to insertion break point.
                --map_ref_pos;
                ab.ref_base  = fa_seq[aligned_pairs[i].ref_pos-1]; // break point's ref base
                ab.read_base = fa_seq[aligned_pairs[i].ref_pos-1] + aligned_pairs[i].read_base;

                // Need to convert the read_base to upper case if it is a insertion in case of the ref is lower case.
                std::transform(ab.read_base.begin(), ab.read_base.end(), ab.read_base.begin(), ::toupper);

                // mean quality of the whole insertion sequence
                double total_score = 0;
                for (size_t i = 0; i < aligned_pairs[i].read_base.size(); ++i)
                    total_score += aligned_pairs[i].read_qual[i];
                ab.base_qual = uint8_t(total_score / aligned_pairs[i].read_base.size());

            } else if (aligned_pairs[i].op == BAM_CDEL) {  /* CIGAR: D */
                // Deletion.
                if (!aligned_pairs[i].read_base.empty()) {
                    std::cerr << al << "\n";
                    throw std::runtime_error("[ERROR] We got read bases in deletion region.");
                }

                // roll back one position to the leftmost to deletion break point.
                --map_ref_pos;
                ab.ref_base  = fa_seq[aligned_pairs[i].ref_pos-1] + aligned_pairs[i].ref_base;
                ab.read_base = fa_seq[aligned_pairs[i].ref_pos-1]; // break point's ref base

                // Need to convert the read_base to upper case if it is a deletion in case of the ref is lower case.
                std::transform(ab.read_base.begin(), ab.read_base.end(), ab.read_base.begin(), ::toupper);

                // set to be mean quality of the whole read if deletion
                ab.base_qual = uint8_t(al.mean_qqual()) + 33; // 33 is the offset of base QUAL;

            } else { 
                // Skip other kinds of CIGAR symbals.
                continue;
            }
            // qpos is 0-based, conver to 1-based to set the rank of base on read.
            ab.rpr = aligned_pairs[i].qpos + 1;

            // 以 map_ref_pos 为 key，将所有的 read_bases 信息存入 map 中，多个突变共享同个 ref_pos，
            if (sample_posinfo_map.find(map_ref_pos) == sample_posinfo_map.end()) {
                // First level. If the position is not in the map, insert it.
                AlignInfo pos_align_info(gr.ref_id, map_ref_pos);
                sample_posinfo_map.insert({map_ref_pos, pos_align_info});
            }
            
            // collected all the base informations of the read into the map.
            sample_posinfo_map[map_ref_pos].align_bases.push_back(ab);
        }
    }

    return;
}

VariantInfo basetype_caller_unit(const AlignInfo &pos_align_info, double min_af) {

    BaseType::BatchInfo smp_bi(pos_align_info.ref_id, pos_align_info.ref_pos);
    for (auto &ab: pos_align_info.align_bases) {
        smp_bi.ref_bases.push_back(ab.ref_base);
        smp_bi.mapqs.push_back(ab.mapq);
        smp_bi.map_strands.push_back(ab.map_strand);

        smp_bi.base_pos_ranks.push_back(ab.rpr);
        smp_bi.align_bases.push_back(ab.read_base);
        smp_bi.align_base_quals.push_back(ab.base_qual);
    }

    BaseType bt(&smp_bi, min_af);
    bt.lrt();  // use likelihood ratio test to detect candidate variant

    return get_pileup(bt, &smp_bi);
}

VariantInfo get_pileup(const BaseType &bt, const BaseType::BatchInfo *smp_bi) {

    VariantInfo vi(bt.get_ref_id(), bt.get_ref_pos(), bt.get_total_depth(), bt.get_var_qual());
    for (size_t i(0); i < bt.get_active_bases().size(); ++i) {
        std::string b = bt.get_active_bases()[i];
        std::string ref_base = bt.get_bases2ref().at(b);

        vi.ref_bases.push_back(ref_base);
        vi.alt_bases.push_back(b);
        vi.depths.push_back(bt.get_base_depth(b));
        vi.freqs.push_back(bt.get_lrt_af(b));

        std::string upper_ref_base(ref_base);
        std::transform(upper_ref_base.begin(), upper_ref_base.end(), upper_ref_base.begin(), ::toupper);

        if (b == upper_ref_base) {
            vi.var_types.push_back("REF");
        } else if (b.size() == 1 && ref_base.size() == 1) {
            vi.var_types.push_back("SNV");
        } else if (b.size() > ref_base.size()) {
            vi.var_types.push_back("INS");
        } else if (b.size() < ref_base.size()) {
            vi.var_types.push_back("DEL");
        } else {
            vi.var_types.push_back("MNV");
        }

        std::pair<double, double> ci = calculate_confidence_interval(bt.get_base_depth(b), bt.get_total_depth());
        vi.ci.push_back(ci);
        vi.strand_bias.push_back(strand_bias(upper_ref_base, b, smp_bi->align_bases, smp_bi->map_strands));
    }

// std::cout << "*****: " << smp_bi->ref_id << "\t" << smp_bi->ref_pos << "\t" 
//           << ngslib::join(smp_bi->ref_bases, ",") << "\t" << ngslib::join(smp_bi->align_bases, ",") << "\t: " 
//           << ngslib::join(bt.get_active_bases(), "-") << "\t" 
//           << ngslib::join(vi.ref_bases, ",") << "\t" << ngslib::join(vi.alt_bases) << "\t" 
//           << ngslib::join(vi.depths, ",") << "\t" << ngslib::join(vi.freqs, ",") << "\t" << vi.total_depth << "-" << bt.get_total_depth() << "\n";
    return vi;
}

VCFRecord call_variant_in_pos(std::vector<VariantInfo> vi) {
    // 1. call the variant by the integrated information
    // 2. output the variant information to the VCF file
    if (vi.empty()) {
        return VCFRecord(); // Return empty record if no variants
    }

    VCFRecord vcf_record;
    const auto& first_var = vi[0];  // Use first variant info as reference
    
    // Set basic VCF fields
    vcf_record.chrom = first_var.ref_id;
    vcf_record.pos = first_var.ref_pos;
    
    // First pass: collect all REF sequences and find the longest one
    std::string shared_ref;  // The final REF to use in VCF
    for (size_t i = 0; i < vi.size(); i++) {
        for (const auto& ref : vi[i].ref_bases) {
            if (ref.length() > shared_ref.length()) {
                shared_ref = ref;
            }
        }
    }

    // Set the shared REF
    vcf_record.ref = shared_ref; // Set the original shared REF sequence
    std::transform(shared_ref.begin(), shared_ref.end(), shared_ref.begin(), ::toupper); // Convert to upper case

    // Second pass: collect and normalize ALT sequences
    std::set<std::string> unique_alts;
    for (size_t i = 0; i < vi.size(); i++) { // Loop all samples in the position, and collect the ALT information
        for (size_t j = 0; j < vi[i].alt_bases.size(); j++) {
            std::string alt = vi[i].alt_bases[j];
            std::string ref = vi[i].ref_bases[j];
            std::transform(ref.begin(), ref.end(), ref.begin(), ::toupper); // Convert to upper case
            
            // Normalize ALT sequence
            if (ref != shared_ref && shared_ref.length() > ref.length()) {
                alt += shared_ref.substr(ref.length());
                vi[i].alt_bases[j] = alt; // Update the ALT sequence
            }
            unique_alts.insert(alt);
        }
    }
    unique_alts.erase(shared_ref); // Remove REF from ALTs

    if (unique_alts.empty()) {
        return VCFRecord(); // Return empty record if no variants
    }

    // Set ALT field: Unique and sorted ALT sequences by length and then by ASCII
    vcf_record.alt = ngslib::get_unique_strings(std::vector<std::string>(unique_alts.begin(), unique_alts.end()));

    std::map<std::string, int> ac; // only record ref and non-ref allele counts (AC)
    double an = 0;                 // Total allele number in called GT

    std::vector<std::string> ref_alt_order = {shared_ref};
    for (size_t i = 0; i < vcf_record.alt.size(); i++) {
        ref_alt_order.push_back(vcf_record.alt[i]);
        ac[vcf_record.alt[i]] = 0; // Initialize AC
    }

    // Set QUAL field: The biggest QUAL value of all samples
    vcf_record.qual = 0;
    
    // Process sample information
    vcf_record.format = "GT:GQ:DP:AD:HF:CI:SB:FS:SOR:VT";
    for (size_t sample_idx = 0; sample_idx < vi.size(); sample_idx++) {
        const auto& sample_var = vi[sample_idx];

        // record the biggest qual value
        if (sample_var.qual > vcf_record.qual && sample_var.qual != 10000) {
            vcf_record.qual = sample_var.qual;
        }

        // Re-order ALTs present in this sample according to 'vcf_record.alt'
        std::vector<size_t> gt_indices; // Genotype indices
        std::vector<std::string> sample_alts;
        std::vector<int>    allele_depths;
        std::vector<double> allele_freqs;

        std::vector<std::string> ci_strings;
        std::vector<std::string> sb_strings;
        std::vector<std::string> fs_strings;
        std::vector<std::string> sor_strings;
        std::vector<std::string> var_types;

        // Collect and format sample information 
        for (size_t gti = 0; gti < ref_alt_order.size(); gti++) {  // gti == 0 represents the REF GT for mtDNA
            const auto& alt = ref_alt_order[gti];

            for (size_t j = 0; j < sample_var.alt_bases.size(); j++) {
                if (sample_var.alt_bases[j] == alt) {

                    an += 1;
                    ac[alt] += 1;

                    sample_alts.push_back(alt);
                    gt_indices.push_back(gti);
                    allele_depths.push_back(sample_var.depths[j]);
                    // allele_freqs.push_back(sample_var.freqs[j]); // 这里不要用 lrt 计算出来的 allele frequency，因为可能不知为何会有负数（极少情况下）
                    // use the allele frequency calculated by allele_depth/total_depth
                    allele_freqs.push_back(double(sample_var.depths[j])/double(sample_var.total_depth)); 
                    ci_strings.push_back(format_double(sample_var.ci[j].first) + "," + format_double(sample_var.ci[j].second));
                    sb_strings.push_back(std::to_string(sample_var.strand_bias[j].ref_fwd) + "," + 
                                         std::to_string(sample_var.strand_bias[j].ref_rev) + "," + 
                                         std::to_string(sample_var.strand_bias[j].alt_fwd) + "," + 
                                         std::to_string(sample_var.strand_bias[j].alt_rev));
                    fs_strings.push_back(sample_var.strand_bias[j].fs != 10000 ? 
                                         format_double(sample_var.strand_bias[j].fs) : "10000"); // it's a phred-scale score
                    sor_strings.push_back(sample_var.strand_bias[j].sor != 10000 ? 
                                          format_double(sample_var.strand_bias[j].sor) : "10000");
                    var_types.push_back(sample_var.var_types[j]);
                }
            }
        }

        bool alt_found = sample_alts.size() > 0;
        std::string gt = alt_found ? ngslib::join(gt_indices, "/") : ".";       // set generate genotype (GT)
        std::string sample_info = gt + ":" +                                    // GT, genotype
                                  std::to_string(int(sample_var.qual)) + ":" +  // GQ, genotype quality
                                  std::to_string(sample_var.total_depth);       // DP, total depth
        if (alt_found) {
            sample_info += ":";
            sample_info += ngslib::join(allele_depths, ",") + ":" +             // AD, allele depth
                           ngslib::join(allele_freqs, ",")  + ":" +             // HF, allele frequency
                           ngslib::join(ci_strings, ";")    + ":" +             // CI, confidence interval
                           ngslib::join(sb_strings, ";")    + ":" +             // SB, strand bias
                           ngslib::join(fs_strings, ",")    + ":" +             // FS, Fisher strand bias
                           ngslib::join(sor_strings, ",")   + ":" +             // SOR, Strand odds ratio
                           ngslib::join(var_types, ",");                        // VT, Variant type
        }
        vcf_record.samples.push_back(sample_info);
    }
    
    // Set INFO field
    std::vector<int> ac_v;
    std::vector<std::string> af;
    for (const auto& alt : vcf_record.alt) { // Only record the counts (AC) and frequencies (AF) of non-ref allele in INFO field
        ac_v.push_back(ac[alt]);
        af.push_back(format_double(ac[alt] / an, 4)); 
    }
    vcf_record.info = "AF=" + ngslib::join(af, ",") + ";AC=" + ngslib::join(ac_v, ",") + ";AN=" + std::to_string(int(an));

    return vcf_record;
}